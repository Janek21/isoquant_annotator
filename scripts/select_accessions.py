#!/usr/bin/env python
# Select the best runs for a species from the LyRic master TSV.
# Filters to transcriptomic bulk, size-caps via ENA, tags developmental stage,
# and picks the highest-read-count runs per platform × stage.
# ENA URL resolution happens at download time in srr_dw.sh.

import argparse
import urllib.request
import pandas as pd

SPECIES_KEYWORDS = {
    "Plasmodium":      ["sporozoite", "merozoite", "ring", "trophozoite", "schizont", "gametocyte", "ookinete"],
    "Eimeria":         ["oocyst", "sporozoite", "merozoite", "microgamete", "macrogamete"],
    "Toxoplasma":      ["tachyzoite", "bradyzoite", "sporozoite", "oocyst"],
    "Sarcocystis":     ["sarcocyst", "bradyzoite", "merozoite", "sporocyst"],
    "Cryptosporidium": ["oocyst", "sporozoite", "merozoite", "gamont"],
    "Paramecium":      ["vegetative", "conjugation", "autogamy"],
    "Tetrahymena":     ["vegetative", "starvation", "conjugation"],
}


def get_species_keywords(species):
    for genus, words in SPECIES_KEYWORDS.items():
        if genus.lower() in str(species).lower():
            return words
    return []


def keyword_extraction(description, keyword_list):
    found = [w.capitalize() for w in keyword_list if w in description.lower()]
    if found:
        return ", ".join(found)
    if "hpi" in description or "hrs post infection" in description:
        return "Timecourse_HPI"
    if "mixed" in description:
        return "Mixed stages"
    return "Unspecified"


def get_size_gb(run_id, timeout=20):
    """Return total fastq size in GB from ENA, or None on failure."""
    api = ("https://www.ebi.ac.uk/ena/portal/api/filereport"
           f"?accession={run_id}&result=read_run&fields=fastq_bytes&format=tsv")
    
    try:
        with urllib.request.urlopen(api, timeout=timeout) as resp:
            lines = resp.read().decode().splitlines()
    except Exception as e:
        print(f"Warning: ENA size lookup failed for {run_id} ({e}). Skipping.")
        return None
        
    # Check if we got actual data beyond the header row
    if len(lines) < 2:
        return None
        
    for row in lines[1:]:
        if row.strip():
            cols = row.split("\t")
            # cols[0] is the run_accession, cols[1] is fastq_bytes
            if len(cols) > 1 and cols[1]:
                total = sum(int(b) for b in cols[1].split(";") if b.isdigit())
                real_size=total / (1024 ** 3)
                print(f"    Size for {run_id}: {round(real_size, 3)} GB")
                return real_size
                
    return None


def select(df, top_reads, max_gb):
    # bulk transcriptomic only
    df = df[df["Source"].str.contains("TRANSCRIPTOMIC", case=False, na=False)].reset_index(drop=True)

    # SRA_id column is "ExperimentID:RunID" – keep only the run ID
    df[["Exp_id", "SRA_id"]] = df["SRA_id"].str.split(":", n=1, expand=True)

    # size filter (one lightweight ENA call per candidate)
    print(f"Checking sizes for {len(df)} candidate runs (cap: {max_gb} GB) ...")
    sizes = df["SRA_id"].apply(get_size_gb)
    df = df[sizes.notna() & (sizes <= max_gb)].reset_index(drop=True)
    print(f"Runs remaining after size filter: {len(df)}")
    if df.empty:
        return df

    keywords = get_species_keywords(df["Species"].iloc[0])
    df["Tissue_stage"] = df["Description"].apply(lambda d: keyword_extraction(d, keywords))

    # top-N highest-read-count runs per platform × developmental stage
    df = df.sort_values(by=["Platform", "Tissue_stage", "Read_count"], ascending=[True, True, False])
    best = df.groupby(["Platform", "Tissue_stage"]).head(top_reads).reset_index(drop=True)
    return best


def main():
    p = argparse.ArgumentParser(description="Select runs for IsoQuant from the LyRic master TSV.")
    p.add_argument("-i", "--input",    required=True, help="Species-filtered master TSV (no header).")
    p.add_argument("-o", "--output",   required=True, help="Selected runs table (srr_select.tsv).")
    p.add_argument("-s", "--srr_list", required=True, help="Plain accession list for the downloader (srr_list.tsv).")
    p.add_argument("-t", "--topReads", type=int,   default=2,    help="Top runs per platform × stage (default 2).")
    p.add_argument("-m", "--max_size", type=float, default=9e15, help="Max run size in GB (default: no limit).")
    args = p.parse_args()

    colnames = ["SRA_id", "Description", "TaxonID", "Lineage", "Species","Source", "Strategy", "Platform", "Read_count", "Date"]
    data = pd.read_csv(args.input, sep="\t", header=None, names=colnames)

    best = select(data, args.topReads, args.max_size)
    if best.empty:
        print("No runs selected.")
        return

    out_cols = ["SRA_id", "Description", "Tissue_stage", "TaxonID", "Lineage", "Species", "Source", "Strategy", "Platform", "Read_count", "Date"]
    best[out_cols].to_csv(args.output, sep="\t", index=False, header=False)
    print(f"Saved {len(best)} runs to {args.output}.")

    # plain accession list – one per line, consumed by srr_dw.sh
    best["SRA_id"].to_csv(args.srr_list, index=False, header=False)
    print(f"Saved accession list to {args.srr_list}.")


if __name__ == "__main__":
    main()
