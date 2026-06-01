#!/usr/bin/env bash
#SBATCH --job-name=iq_eval
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=90
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# Merges the per-platform IsoQuant outputs (pacbio + nanopore) and evaluates
# the merged annotation with BUSCO.
# Usage: sbatch evaluation.sh <species_name> [busco_db]
set -euo pipefail

species_name="$1"
busco_db="${2:-/no_backup/rg/references/busco_downloads}"

# location of get_busco_db.py (from the LyRic repo) and an NCBI Entrez email
ncbi_email="${NCBI_EMAIL:-your_email@example.com}"

sp=$(echo "$species_name" | cut -f2 -d"_")
genome=$(ls "$species_name"/data/fasta/*.fa* 2>/dev/null | grep -vE '\.fai$' | head -1)
out="$species_name/output/eval"
mkdir -p "$out" logs

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate buscomania

echo "════════════════════════════════════════"
echo " IsoQuant evaluation: $species_name"
echo " genome : $genome"
echo "════════════════════════════════════════"

# ── 1. gather per-platform GTFs and merge (pass through if only one) ──
gtfs=()
for plat in pacbio nanopore; do
	g=$(find "$species_name/output/$plat" -name "*transcript_models.gtf" 2>/dev/null | head -1 || true)
	if [ -n "$g" ]; then
		echo "[1/5] Found $plat model: $g"
		gtfs+=("$g")
	fi
done

if [ "${#gtfs[@]}" -eq 0 ]; then
	echo "ERROR: no transcript_models.gtf under output/pacbio or output/nanopore"
	exit 1
elif [ "${#gtfs[@]}" -eq 1 ]; then
	echo "[1/5] Single platform; no merge needed."
	merged="${gtfs[0]}"
else
	echo "[1/5] Merging ${#gtfs[@]} platforms with AGAT ..."
	merged="$out/merged_${sp}.gff"
	merge_args=()
	for g in "${gtfs[@]}"; do merge_args+=(--gff "$g"); done
	agat_sp_merge_annotations.pl "${merge_args[@]}" --out "$merged"
fi
echo "    Annotation to evaluate: $merged"

# ── 2. longest isoform per gene ───────────────────────────────────
echo "[2/5] Extracting longest isoforms ..."
agat_sp_keep_longest_isoform.pl --gff "$merged" --out "$out/longest_${sp}.gtf"

# ── 3. transcriptome with gffread ─────────────────────────────────
echo "[3/5] Building transcriptome ..."
# gffread cannot read gzipped FASTA; decompress to a temp file if needed
if [[ "$genome" == *.gz ]]; then
	genome_plain="$out/genome_tmp.fa"
	echo "    Decompressing genome for gffread ..."
	pigz -dcp "${SLURM_CPUS_PER_TASK:-8}" "$genome" > "$genome_plain"
else
	genome_plain="$genome"
fi

gffread "$out/longest_${sp}.gtf" -g "$genome_plain" -w "$out/transcripts_${sp}.fa"

if [[ "$genome" == *.gz ]]; then
	rm -f "$genome_plain"
fi

# ── 4. ORF prediction with TransDecoder ───────────────────────────
echo "[4/5] Predicting ORFs (TransDecoder) ..."
TD2.LongOrfs -t "$out/transcripts_${sp}.fa" -O "$out/td_work"
TD2.Predict -t "$out/transcripts_${sp}.fa" -O "$out/td_work"
mv "./transcripts_${sp}.fa.TD2.pep" "$out/prot_${sp}.fa"
mv ./*.fa.TD2.* "$out/td_work" 2>/dev/null || true

# ── 5. BUSCO (taxon-driven lineage) ───────────────────────────────
# most frequent TaxonID among the selected runs
taxonID=$(cut -f4 "$species_name/srr_select.tsv" | sort | uniq -c | sort -nr | awk '{print $2}' | head -n1)
echo "[5/5] TaxonID: $taxonID"

busco_lineage=$(python3 "scripts/get_busco_db.py" \
	-e "$ncbi_email" \
	-t "$taxonID" \
	-b "$busco_db/file_versions.tsv" \
	-v odb12)
echo "      BUSCO lineage: $busco_lineage"

if [ -z "$busco_lineage" ]; then
	echo "ERROR: could not resolve a BUSCO lineage for TaxonID $taxonID"
	exit 1
fi

busco -i "$out/prot_${sp}.fa" \
	-o "busco_${sp}" \
	--out_path "$out" \
	-m protein \
	-l "$busco_lineage" \
	--download_path "$busco_db" \
	-c "${SLURM_CPUS_PER_TASK:-8}" \
	-f

# ── 6. collect the BUSCO JSON summary into a shared folder ─────────
busco_summary_dir="busco_summary"
mkdir -p "$busco_summary_dir"
busco_json="$out/busco_${sp}/short_summary.specific.${busco_lineage}.busco_${sp}.json"
busco_json_dest="$busco_summary_dir/${species_name}_${taxonID}_busco.json"
mv "$busco_json" "$busco_json_dest"
ln "$busco_json_dest" "$busco_json"   #keep it accessible at the original BUSCO output location too
echo "[6/6] BUSCO JSON summary collected into $busco_summary_dir/"

# count gene and transcript models in the prediction (col3 feature type;
# IsoQuant GTF uses "transcript", AGAT GFF uses "mRNA" — match both)
gene_count=$(cut -f3 "$merged" | grep -cxF "gene" || true)
transcript_count=$(cut -f3 "$merged" | grep -cxE 'transcript|mRNA' || true)
echo "$gene_count" > "$busco_summary_dir/${species_name}_${taxonID}_gc.txt"
echo "$transcript_count" > "$busco_summary_dir/${species_name}_${taxonID}_tc.txt"
echo "      Gene models: $gene_count | Transcript models: $transcript_count"

rm -rf agat_log_*
echo "Done. Merged annotation: $merged"
echo "BUSCO results in: $out/busco_${sp}/"
echo "BUSCO JSON summaries in: $busco_summary_dir/"
