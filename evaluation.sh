#!/usr/bin/env bash
#SBATCH --job-name=iq_eval
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=90
#SBATCH --output=logs/eval/%x_%j.out
#SBATCH --error=logs/eval/%x_%j.err
#
# Merges the per-platform IsoQuant outputs (pacbio + nanopore) and evaluates
# the merged annotation with BUSCO.
# Usage: sbatch evaluation.sh <species_name> [busco_db]
set -euo pipefail

echo ">STARTING at $(date)"

species_name="$1"
busco_db="${2:-/no_backup/rg/references/busco_downloads}"

#NCBI Entrez email
ncbi_email="${NCBI_EMAIL:-nqvsisnkflvflitqoy@kjkpc.net}"

sp=$(echo "$species_name" | cut -f2 -d"_")
genome=$(ls "$species_name"/data/fasta/*.fa* 2>/dev/null | grep -vE '\.fai$')
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
		# AGAT derives its temp dir name (agat_tmp_<input-basename>) in the CWD.
		# Every IsoQuant model is OUT.transcript_models.gtf, so concurrent eval
		# jobs from the shared submit dir collide on agat_tmp_OUT.transcript_models.
		# Feed AGAT a uniquely named symlink so each job gets its own temp dir.
		link="$out/${species_name}_${plat}.transcript_models.gtf"
		ln -sf "$(realpath "$g")" "$link"
		gtfs+=("$link")
	fi
done

if [ "${#gtfs[@]}" -eq 0 ]; then
	echo "ERROR: no transcript_models.gtf under output/pacbio or output/nanopore"
	exit 1
elif [ "${#gtfs[@]}" -eq 1 ]; then
	echo "[1/5] Single platform; converting GTF to GFF3 ..."
	merged="$out/merged_${sp}.gff"
	agat_convert_sp_gxf2gxf.pl -g "${gtfs[0]}" --config "$agat_cfg" -o "$merged"
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

# IsoQuant's reference-free mode can emit a transcript model that runs past the
# contig end (e.g. kawagutii: a 1..244 model on the 232 bp contig
# VSDK01014088.1). gffread then aborts the whole transcriptome build with
# "GFaSeqGet: subsequence cannot be larger than ...". Collect the IDs of any
# out-of-bounds models and let gffread discard them (--nids), so it handles the
# gene/transcript/exon hierarchy itself. Valid species yield an empty list and
# the gffread call is unchanged.
oob_ids="$out/oob_${sp}.ids"
awk -v genome="$genome_plain" '
	BEGIN {
		# contig lengths straight from the genome FASTA (no samtools needed)
		while ((getline line < genome) > 0) {
			if (substr(line, 1, 1) == ">") {
				split(substr(line, 2), h, /[ \t]/); cur = h[1]; len[cur] = 0
			} else {
				len[cur] += length(line)
			}
		}
		close(genome)
	}
	/^#/ { next }
	($1 in len) && ($5 + 0 > len[$1]) && match($9, /transcript_id=[^;]+/) {
		print substr($9, RSTART + 14, RLENGTH - 14)
	}
' "$out/longest_${sp}.gtf" | sort -u > "$oob_ids"

if [ -s "$oob_ids" ]; then
	echo "    Dropping $(wc -l < "$oob_ids") out-of-bounds transcript model(s) before gffread"
	gffread --nids "$oob_ids" "$out/longest_${sp}.gtf" -g "$genome_plain" -w "$out/transcripts_${sp}.fa"
else
	gffread "$out/longest_${sp}.gtf" -g "$genome_plain" -w "$out/transcripts_${sp}.fa"
fi

rm -f "$oob_ids"
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
td_work="$(realpath "$out/td_work")"
mkdir -p "$td_work"
transcripts_abs="$(realpath "$out/transcripts_${sp}.fa")"   #transcriptome built in step 3

(cd "$td_work" && #move to folder for TD2 execution ONLY
	#Find ORFs in transcripts
	TD2.LongOrfs -t "$transcripts_abs" -O $td_work -G "$gcode"
	#Select most probable ORFs to create proteins
	TD2.Predict -t "$transcripts_abs" -O $td_work -G "$gcode" -v #-O is output of ORFs
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
pred_dir="$summary_dir/pred"
mkdir -p "$busco_lineage_dir" "$busco_euk_dir" "$counts_dir" "$pred_dir"

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

# predicted (merged) annotation relocated to summary/, hardlinked back to the species location
# (canonical inode lives in summary/pred so the species folder can be removed safely)
pred_dest="$pred_dir/${species_name}_${taxonID}_pred.gff"
rm -f "$pred_dest"                  #refresh on reruns
mv "$merged" "$pred_dest"          #relocate the prediction into the central summary tree
ln "$pred_dest" "$merged"          #link it back so the original species location stays valid

echo "[6/6] BUSCO JSON summaries collected into $busco_lineage_dir/ and $busco_euk_dir/"
echo "      Predicted annotation collected into $pred_dir/"

# count gene/transcript models with gffread. --keep-genes normalises every input
# into gene + transcript records: real gene features are preserved (so AGAT
# clustered loci are honoured rather than the per-transcript gene_id attribute,
# which --table @geneid would overcount), and one gene + one transcript is
# synthesised per id when the input lacks gene/transcript features. Counting the
# normalised col3 feature types is then uniform across annotation flavours.
read -r gene_count transcript_count < <(
	{ gffread "$merged" --keep-genes -o - 2>/dev/null || true; } | awk -F'\t' '
		/^#/ { next }
		$3 ~ /^([A-Za-z_]*gene)$/                { g++; next }
		$3 ~ /^(transcript|mRNA|[A-Za-z_]*RNA)$/ { t++ }
		END { print g + 0, t + 0 }'
)
echo "$gene_count" > "$counts_dir/${species_name}_${taxonID}_gc.txt"
echo "$transcript_count" > "$counts_dir/${species_name}_${taxonID}_tc.txt"
echo "      Gene models: $gene_count | Transcript models: $transcript_count"

# genome size = total assembly length (exact; sum of contig lengths, incl. N gaps)
fai="${genome}.fai"
if [ -s "$fai" ]; then
	genome_size=$(cut -f2 "$fai" | awk '{s+=$1} END{print s+0}')
elif [[ "$genome" == *.gz ]]; then
	genome_size=$(pigz -dcp "${SLURM_CPUS_PER_TASK:-8}" "$genome" \
		| awk '/^>/{next} {s+=length($0)} END{print s+0}')
else
	genome_size=$(awk '/^>/{next} {s+=length($0)} END{print s+0}' "$genome")
fi
echo "$genome_size" > "$counts_dir/${species_name}_${taxonID}_gs.txt"
echo "      Genome size: ${genome_size} bp"

echo "Done. Merged annotation: $merged"
echo "BUSCO results in: $out/busco_${sp}/ and $out/busco_euk_${sp}/"
echo "Summary outputs in: $summary_dir/ (busco_lineage/, busco_eukaryote/, counts/)"
echo "Build the summary tables with: python3 scripts/make_summary_tables.py"
echo ">ENDING at $(date)"
