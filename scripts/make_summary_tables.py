#!/usr/bin/env python3
"""Build summary tables from the per-species files written by evaluation.sh.

Reads the shared summary/ tree:
  summary/counts/<species>_<taxonID>_metrics.tsv          one row of counts + derived metrics
  summary/busco_lineage/<species>_<taxonID>_Lbusco.json    taxon-driven lineage BUSCO
  summary/busco_eukaryote/<species>_<taxonID>_Ebusco.json  Eukaryota BUSCO

Each _metrics.tsv carries (header + one data row):
  species, gene_count, transcript_count, genome_size_bp,
  coding_transcripts, transcriptome_transcripts,
  gene_density_per_mb, transcript_density_per_mb, isoforms_per_gene, coding_fraction

Writes three TSVs into the summary directory:
  1. counts_summary.tsv    every metrics row, one species per line
  2. busco_summary.tsv     lineage_busco, eukaryote_busco, lineage_used (Completeness only)
  3. general_summary.tsv   the first two tables merged on the species column

Run manually once all evaluation jobs have finished:
  python3 scripts/make_summary_tables.py [summary_dir]
"""
import glob
import json
import os
import re
import sys

import pandas as pd

METRIC_COLUMNS = [
    "species",
    "gene_count",
    "transcript_count",
    "genome_size_bp",
    "coding_transcripts",
    "transcriptome_transcripts",
    "gene_density_per_mb",
    "transcript_density_per_mb",
    "isoforms_per_gene",
    "coding_fraction",
]


def busco_fields(path):
    """Return (completeness_C, lineage_name) from a BUSCO short_summary JSON.

    Only the Completeness score is kept from the one-line summary, e.g.
    'C:60.8%[S:40.2%,D:20.6%],F:9.1%,M:30.2%,n:1591' -> '60.8'.
    """
    try:
        with open(path) as fh:
            data = json.load(fh)
    except (FileNotFoundError, json.JSONDecodeError):
        return None, None
    results = data.get("results", {})
    match = re.search(r"C:([\d.]+)%", results.get("one_line_summary", "") or "")
    if match:
        completeness = match.group(1)
    else:
        complete = results.get("Complete percentage")
        completeness = str(complete) if complete is not None else None
    lineage = data.get("lineage_dataset", {}).get("name")
    return completeness, lineage


def species_keys(directory, suffix):
    """Species identifiers (<species>_<taxonID>) from files ending in suffix."""
    keys = set()
    for path in glob.glob(os.path.join(directory, "*" + suffix)):
        keys.add(os.path.basename(path)[: -len(suffix)])
    return keys


def read_metrics(counts_dir):
    """Concatenate every per-species _metrics.tsv into one DataFrame."""
    frames = []
    for path in sorted(glob.glob(os.path.join(counts_dir, "*_metrics.tsv"))):
        try:
            frames.append(pd.read_csv(path, sep="\t", dtype=str))
        except (pd.errors.EmptyDataError, FileNotFoundError):
            continue
    if not frames:
        return pd.DataFrame(columns=METRIC_COLUMNS)
    return pd.concat(frames, ignore_index=True).reindex(columns=METRIC_COLUMNS)


def main():
    summary_dir = sys.argv[1] if len(sys.argv) > 1 else "summary"
    counts_dir = os.path.join(summary_dir, "counts")
    lineage_dir = os.path.join(summary_dir, "busco_lineage")
    euk_dir = os.path.join(summary_dir, "busco_eukaryote")

    metrics_df = read_metrics(counts_dir)

    busco_species = (species_keys(lineage_dir, "_Lbusco.json")
                     | species_keys(euk_dir, "_Ebusco.json"))
    busco_records = []
    for sp in sorted(busco_species):
        lineage_busco, lineage_used = busco_fields(
            os.path.join(lineage_dir, f"{sp}_Lbusco.json"))
        eukaryote_busco, _ = busco_fields(
            os.path.join(euk_dir, f"{sp}_Ebusco.json"))
        busco_records.append({
            "species": sp,
            "lineage_busco": lineage_busco,
            "eukaryote_busco": eukaryote_busco,
            "lineage_used": lineage_used,
        })
    busco_df = pd.DataFrame(
        busco_records,
        columns=["species", "lineage_busco", "eukaryote_busco", "lineage_used"])

    if metrics_df.empty and busco_df.empty:
        print(f"No result files found under {summary_dir}/. Nothing to do.")
        return

    # third table: merge the two tables on the species column (outer keeps
    # metrics-only and BUSCO-only species)
    general_df = metrics_df.merge(busco_df, on="species", how="outer")

    counts_path = os.path.join(summary_dir, "counts_summary.tsv")
    busco_path = os.path.join(summary_dir, "busco_summary.tsv")
    general_path = os.path.join(summary_dir, "general_summary.tsv")

    metrics_df.sort_values("species").to_csv(
        counts_path, sep="\t", index=False, na_rep="NA")
    busco_df.sort_values("species").to_csv(
        busco_path, sep="\t", index=False, na_rep="NA")
    general_df.sort_values("species").to_csv(
        general_path, sep="\t", index=False, na_rep="NA")

    n = general_df["species"].nunique()
    print(f"Wrote 3 tables for {n} species to {summary_dir}/:")
    print("  counts_summary.tsv, busco_summary.tsv, general_summary.tsv")


if __name__ == "__main__":
    main()
