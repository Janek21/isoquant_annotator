#!/bin/bash
# Prepare a species, download its reads from ENA, and launch IsoQuant
# automatically once the downloads finish (SLURM dependency).
# Usage: bash isoquant_prepare.sh <Genus_species> [master_tsv]

species_name="$1"
master_tsv="${2:-../data/longread_protists.tsv}"
busco_db="${3:-}"   # forwarded to evaluation.sh (empty = use its default)
scripts_dir="scripts"

sp=$(echo "$species_name" | cut -f2 -d"_")
sp_extra=$(echo "$species_name" | cut -f3 -d"_")
data="$species_name/data"
mkdir -p "$species_name" logs

#1. Copy assembly and annotation if its present
mkdir -p "$species_name/data/fasta"

#find genome and uncompress if needed
plain_genome=$(find ../data/species/"$species_name"*/GC* -type f -name "GC*.fna")
gzipped_genome=$(find ../data/species/"$species_name"*/GC* -type f -name "GC*.fna.gz")
if [ -n "$plain_genome" ]; then
	echo "Genome uncompressed"
	#copy regular genome to species data folder
	cp "$plain_genome" "$species_name/data/fasta/assembly_${sp}_genomic.fa"
else
	echo "Genome compressed"
	#uncompress genome to species data folder
	unpigz -c "$gzipped_genome" > "$species_name/data/fasta/assembly_${sp}_genomic.fa"
fi

#annotation
mkdir -p "$species_name/data/input"
function gtf_fix() { #gff to gtf conversion
	awk '$3 != "region"' "$species_name/data/input/clean_annotation.gtf" > "$species_name/data/input/clean_annotation.tmp.gtf"
	gffread "$species_name/data/input/clean_annotation.tmp.gtf" -T -o "$species_name/data/input/clean_annotation.fixed.gtf"
	echo "GFF to GTF successful"
}
	
#define the reference database argument for IsoQuant (overwriten if no annot)
genedb_arg="--genedb $species_name/data/input/clean_annotation.fixed.gtf"

plain_annotation=$(find ../data/species/"$species_name"*/GC* -type f \( -name "*GC*.gff" -o -name "*GC*.gff3" \))
gzipped_annotation=$(find ../data/species/"$species_name"*/GC* -type f \( -name "*GC*.gff.gz" -o -name "*GC*.gff3.gz" \))
if [ -n "$plain_annotation" ]; then
	echo "GFF uncompressed to GTF"
	#copy regular annotation to species data folder
	agat_convert_sp_gff2gtf.pl --gff "$plain_annotation" -o "$species_name/data/input/clean_annotation.gtf"
	gtf_fix
elif [ -n "$gzipped_annotation" ]; then
	echo "GFF compressed to GTF"
	#uncompress annotation to species data folder
	unpigz -c "$gzipped_annotation" | agat_convert_sp_gff2gtf.pl --gff - -o "$species_name/data/input/clean_annotation.gtf"
	gtf_fix
else #no annotation
	echo "No annotation file found. IsoQuant will run in de novo mode (without --genedb)."
	genedb_arg=""
fi

#2. select runs from the master TSV and resolve their ENA URLs
echo "Searching master TSV for $sp $sp_extra"
search_res=$(grep -i "$sp" "$master_tsv" | grep -i "$sp_extra")
if [ -z "$search_res" ]; then
	echo "No match for both terms. Falling back to: $sp"
	search_res=$(grep -i "$sp" "$master_tsv")
fi
echo "$search_res" > "$species_name/full_srr.tsv"

python3 "$scripts_dir/select_accessions.py" -i "$species_name/full_srr.tsv" -o "$species_name/srr_select.tsv" -s "$species_name/srr_list.tsv" -t 15 -m 8

srr_count=$(wc -l < "$species_name/srr_list.tsv")
if [ "$srr_count" -eq 0 ]; then
	echo "No runs selected. Stopping."
	exit 1
fi
echo "Selected $srr_count runs."

#3. submit the ENA download array (one task per run)
lyric_dir="$HOME/git/lyric_annotator/$species_name/data/fastq"

#check lyric existence
if [ -d "$lyric_dir" ];then
	mkdir -p "$species_name/data/fastq/"
	#if fastq exist already, copy them
	ln -v "$lyric_dir"/* "$species_name/data/fastq/"
	echo "$species_name data from LyRic."
	dl_jobid=""

else #if no data from lyric
	dl_jobid=$(sbatch --parsable \
		--job-name="ena_download_${sp}" \
		--array=0-$(( srr_count - 1 )) \
		--cpus-per-task=2 \
		--mem=8G \
		--time=180 \
		"$scripts_dir/srr_dw.sh" "$species_name")
	echo "Download array submitted: job $dl_jobid"

fi

#4. submit IsoQuant. If a download array was submitted, start after it has finished. With LyRic no dw, so jut submit straight
if [ -n "$dl_jobid" ]; then
	iq_dep=(--dependency=afterany:"$dl_jobid")
else
	iq_dep=()
fi
iq_jobid=$(sbatch --parsable \
	--job-name="isoquant_${sp}" \
	"${iq_dep[@]}" \
	--cpus-per-task=4 \
	--mem=16G \
	--time=260 \
	"$scripts_dir/isoquant_run.sh" "$species_name" "$genedb_arg")
echo "IsoQuant submitted: job $iq_jobid${dl_jobid:+ (starts after job $dl_jobid)}" #only 2nd part if it waits for dependency

#5. merge platforms + evaluate, once both IsoQuant array tasks finish
ev_jobid=$(sbatch --parsable \
	--job-name="iq_eval_${sp}" \
	--dependency=afterok:"$iq_jobid" \
	--cpus-per-task=4 \
	--mem=12G \
	--time=90 \
	"$scripts_dir/evaluation.sh" "$species_name" "$busco_db")
echo "Evaluation submitted: job $ev_jobid (starts after job $iq_jobid)"
