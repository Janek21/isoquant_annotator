#!/usr/bin/env bash
# Aggregate the per-species gene/transcript count files written by evaluation.sh
# (summary/counts/<species>_<taxonID>_gc.txt and _tc.txt) into a single table.
# Run manually once all evaluation jobs have finished.
# Usage: bash scripts/make_counts_summary.sh [summary_dir]
set -euo pipefail

summary_dir="${1:-summary}"
counts_dir="$summary_dir/counts"
out_tsv="$summary_dir/counts_summary.tsv"

if [ ! -d "$counts_dir" ]; then
	echo "ERROR: counts directory not found: $counts_dir" >&2
	exit 1
fi

printf 'species\tgene_count\ttranscripts_count\n' > "$out_tsv"

shopt -s nullglob
n=0
for gc in "$counts_dir"/*_gc.txt; do
	species=$(basename "$gc" _gc.txt)          # <species_name>_<taxonID>
	tc="$counts_dir/${species}_tc.txt"
	gene_count=$(<"$gc")
	transcript_count="NA"
	[ -f "$tc" ] && transcript_count=$(<"$tc")
	printf '%s\t%s\t%s\n' "$species" "$gene_count" "$transcript_count" >> "$out_tsv"
	n=$((n + 1))
done

echo "Wrote $n species row(s) to $out_tsv"
