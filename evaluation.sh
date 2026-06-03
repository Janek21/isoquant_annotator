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

#NCBI Entrez email
ncbi_email="${NCBI_EMAIL:-nqvsisnkflvflitqoy@kjkpc.net}"

sp=$(echo "$species_name" | cut -f2 -d"_")
genome=$(ls "$species_name"/data/fasta/*.fa* 2>/dev/null | grep -vE '\.fai$' | head -1)
out="$species_name/output/eval"
mkdir -p "$out" logs

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate buscomania

#AGAT config so parallel jobs don't collide --no-log so AGAT never creates agat_log_*
agat_cfg="$out/agat_${species_name}_${SLURM_ARRAY_TASK_ID:-$$}.yaml"
agat config --expose --no-log --output "$agat_cfg" >/dev/null 2>&1
trap 'rm -f "$agat_cfg"' EXIT

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
	agat_sp_merge_annotations.pl "${merge_args[@]}" --config "$agat_cfg" --out "$merged"
fi
echo "    Annotation to evaluate: $merged"

# ── 2. longest isoform per gene ───────────────────────────────────
echo "[2/5] Extracting longest isoforms ..."
agat_sp_keep_longest_isoform.pl --gff "$merged" --config "$agat_cfg" --out "$out/longest_${sp}.gtf"

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

# most frequent TaxonID among the selected runs (drives genetic code + BUSCO lineage)
taxonID=$(cut -f4 "$species_name/srr_select.tsv" | sort | uniq -c | sort -nr | awk '{print $2}' | head -n1)
echo "    TaxonID: $taxonID"

#resolve the correct NCBI translation table for this taxon.
#AGAT to use the right table
gcode=$(python3 "scripts/get_genetic_code.py" -e "$ncbi_email" -k "${NCBI_API_KEY:-}" -t "$taxonID" 2>/dev/null)
if ! [[ "$gcode" =~ ^[0-9]+$ ]]; then
	echo "    Could not resolve genetic code for taxon $taxonID; defaulting to table 1."
	gcode=1
fi
echo "    Translation table for $taxonID: $gcode"

# ── 4. ORF prediction with TransDecoder ───────────────────────────
echo "[4/5] Predicting ORFs (TransDecoder, genetic code $gcode) ..."
#ensure files are generated in particular folders(no naming clash)
td_work="$out/td_work"
mkdir -p "$td_work"
transcripts_abs="$(realpath "$out/transcripts_${sp}.fa")"   #transcriptome built in step 3

(cd "$td_work" && #move to folder for TD2 execution ONLY
	#Find ORFs in transcripts
	TD2.LongOrfs -t "$transcripts_abs" -O . -G "$gcode"
	#Select most probable ORFs to create proteins
	TD2.Predict -t "$transcripts_abs" -O . -G "$gcode" #-O is output of ORFs
)
#move TD2 prot files to correct folders
mv "$td_work/transcripts_${sp}.fa.TD2.pep" "$out/prot_${sp}.fa"

# ── 5. BUSCO (taxon-driven lineage + Eukaryota) ───────────────────
echo "[5/5] BUSCO (TaxonID: $taxonID)"

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

euk_lineage="eukaryota_odb12"

# taxon-driven lineage run
busco -i "$out/prot_${sp}.fa" \
	-o "busco_${sp}" \
	--out_path "$out" \
	-m protein \
	-l "$busco_lineage" \
	--download_path "$busco_db" \
	-c "${SLURM_CPUS_PER_TASK:-8}" \
	-f

# Eukaryota run
echo "      Eukaryota BUSCO lineage: $euk_lineage"
busco -i "$out/prot_${sp}.fa" \
	-o "busco_euk_${sp}" \
	--out_path "$out" \
	-m protein \
	-l "$euk_lineage" \
	--download_path "$busco_db" \
	-c "${SLURM_CPUS_PER_TASK:-8}" \
	-f

# ── 6. collect results into the shared summary/ tree ──────────────
summary_dir="summary"
busco_lineage_dir="$summary_dir/busco_lineage"
busco_euk_dir="$summary_dir/busco_eukaryote"
counts_dir="$summary_dir/counts"
mkdir -p "$busco_lineage_dir" "$busco_euk_dir" "$counts_dir"

# taxon-driven lineage BUSCO JSON
Lbusco_json="$out/busco_${sp}/short_summary.specific.${busco_lineage}.busco_${sp}.json"
Lbusco_dest="$busco_lineage_dir/${species_name}_${taxonID}_Lbusco.json"
mv "$Lbusco_json" "$Lbusco_dest"
ln "$Lbusco_dest" "$Lbusco_json"   #keep it accessible at the original BUSCO output location too

# Eukaryota BUSCO JSON
Ebusco_json="$out/busco_euk_${sp}/short_summary.specific.${euk_lineage}.busco_euk_${sp}.json"
Ebusco_dest="$busco_euk_dir/${species_name}_${taxonID}_Ebusco.json"
mv "$Ebusco_json" "$Ebusco_dest"
ln "$Ebusco_dest" "$Ebusco_json"   #keep it accessible at the original BUSCO output location too

echo "[6/6] BUSCO JSON summaries collected into $busco_lineage_dir/ and $busco_euk_dir/"

# count gene and transcript models in the prediction (col3 feature type;
# IsoQuant GTF uses "transcript", AGAT GFF uses "mRNA" — match both)
# the per-species files below are aggregated later by scripts/make_counts_summary.sh
gene_count=$(cut -f3 "$merged" | grep -cxF "gene" || true)
transcript_count=$(cut -f3 "$merged" | grep -cxE 'transcript|mRNA' || true)
echo "$gene_count" > "$counts_dir/${species_name}_${taxonID}_gc.txt"
echo "$transcript_count" > "$counts_dir/${species_name}_${taxonID}_tc.txt"
echo "      Gene models: $gene_count | Transcript models: $transcript_count"

echo "Done. Merged annotation: $merged"
echo "BUSCO results in: $out/busco_${sp}/ and $out/busco_euk_${sp}/"
echo "Summary outputs in: $summary_dir/ (busco_lineage/, busco_eukaryote/, counts/)"
echo "Build the summary tables with: python3 scripts/make_summary_tables.py"
