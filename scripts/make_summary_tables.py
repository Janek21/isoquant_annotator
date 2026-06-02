#!/usr/bin/env python3
"""Build summary tables from the per-species files written by evaluation.sh.

Reads the shared summary/ tree:
  summary/counts/<species>_<taxonID>_gc.txt   gene-model count
  summary/counts/<species>_<taxonID>_tc.txt   transcript-model count
  summary/busco_lineage/<species>_<taxonID>_Lbusco.json    taxon-driven lineage BUSCO
  summary/busco_eukaryote/<species>_<taxonID>_Ebusco.json  Eukaryota BUSCO

Writes three TSVs into the summary directory:
  1. counts_summary.tsv   gene_count, transcripts_count
  2. busco_summary.tsv     lineage_busco, eukaryote_busco, lineage_used  (Completeness % only)
  3. general_summary.tsv   the first two tables merged on the species column

Run manually once all evaluation jobs have finished:
  python3 scripts/make_summary.py [summary_dir]
"""
import csv
import glob
import json
import os
import re
import sys


def read_count(path):
    """Return the integer string stored in a count file, or 'NA' if absent."""
    try:
        with open(path) as fh:
            return fh.read().strip() or "NA"
    except FileNotFoundError:
        return "NA"


def busco_fields(path):
    """Return (completeness_C, lineage_name) from a BUSCO short_summary JSON.

    Only the Completeness score is kept from the one-line summary, e.g.
    'C:60.8%[S:40.2%,D:20.6%],F:9.1%,M:30.2%,n:1591' -> '60.8%'.
    """
    try:
        with open(path) as fh:
            data = json.load(fh)
    except (FileNotFoundError, json.JSONDecodeError):
        return "NA", "NA"
    results = data.get("results", {})
    match = re.search(r"C:([\d.]+%)", results.get("one_line_summary", "") or "")
    if match:
        completeness = match.group(1)
    else:
        complete = results.get("Complete percentage")
        completeness = f"{complete}%" if complete is not None else "NA"
    lineage = data.get("lineage_dataset", {}).get("name", "NA")
    return completeness, str(lineage).strip()


def species_keys(directory, suffix):
    """Species identifiers (<species>_<taxonID>) from files ending in suffix."""
    keys = set()
    for path in glob.glob(os.path.join(directory, "*" + suffix)):
        keys.add(os.path.basename(path)[: -len(suffix)])
    return keys


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def read_tsv(path):
    """Return (header, {species: row}) for a TSV keyed on its first column."""
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        rows = {row[0]: row for row in reader if row}
    return header, rows


def main():
    summary_dir = sys.argv[1] if len(sys.argv) > 1 else "summary"
    counts_dir = os.path.join(summary_dir, "counts")
    lineage_dir = os.path.join(summary_dir, "busco_lineage")
    euk_dir = os.path.join(summary_dir, "busco_eukaryote")

    # union of species across every result type so nothing is dropped
    species = species_keys(counts_dir, "_gc.txt")
    species |= species_keys(lineage_dir, "_Lbusco.json")
    species |= species_keys(euk_dir, "_Ebusco.json")
    species = sorted(species)

    if not species:
        print(f"No result files found under {summary_dir}/. Nothing to do.")
        return

    counts_rows, busco_rows = [], []
    for sp in species:
        gene_count = read_count(os.path.join(counts_dir, f"{sp}_gc.txt"))
        transcripts_count = read_count(os.path.join(counts_dir, f"{sp}_tc.txt"))
        lineage_busco, lineage_used = busco_fields(
            os.path.join(lineage_dir, f"{sp}_Lbusco.json"))
        eukaryote_busco, _ = busco_fields(
            os.path.join(euk_dir, f"{sp}_Ebusco.json"))

        counts_rows.append([sp, gene_count, transcripts_count])
        busco_rows.append([sp, lineage_busco, eukaryote_busco, lineage_used])

    counts_path = os.path.join(summary_dir, "counts_summary.tsv")
    busco_path = os.path.join(summary_dir, "busco_summary.tsv")
    general_path = os.path.join(summary_dir, "general_summary.tsv")

    write_tsv(counts_path,
              ["species", "gene_count", "transcripts_count"], counts_rows)
    write_tsv(busco_path,
              ["species", "lineage_busco", "eukaryote_busco", "lineage_used"], busco_rows)

    # third table: merge the two tables just written, on the species column
    counts_header, counts_data = read_tsv(counts_path)
    busco_header, busco_data = read_tsv(busco_path)
    general_header = counts_header + busco_header[1:]   # drop the duplicate 'species'
    general_rows = []
    for sp in sorted(set(counts_data) | set(busco_data)):
        crow = counts_data.get(sp, [sp] + ["NA"] * (len(counts_header) - 1))
        brow = busco_data.get(sp, [sp] + ["NA"] * (len(busco_header) - 1))
        general_rows.append(crow + brow[1:])
    write_tsv(general_path, general_header, general_rows)

    print(f"Wrote 3 tables for {len(species)} species to {summary_dir}/:")
    print("  counts_summary.tsv, busco_summary.tsv, general_summary.tsv")


if __name__ == "__main__":
    main()
