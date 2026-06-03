#!/bin/bash

#SBATCH --output=logs/run/%x.%A_%a.out
#SBATCH --error=logs/run/%x.%A_%a.err
#SBATCH --job-name=isoquant
#SBATCH --qos=normal
#SBATCH --time=260
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
##array: task 0 = PacBio, task 1 = Nanopore
#SBATCH --array=0-1

start_time=$(date +%s)
echo ">STARTING at $(date)"

source $(conda info --base)/etc/profile.d/conda.sh
conda activate isoquant

species_name="$1" #Cyanidioschyzon_merolae_strain_10D
genedb_arg="$2" #$2 = RefAnn or empty for ref-free mode
data="$species_name/data"
select_tsv="$species_name/srr_select.tsv"   #selected runs: col1=run id, col9=platform
echo "Processing species: $species_name (task $SLURM_ARRAY_TASK_ID)"

if [ -n "$genedb_arg" ]; then
	echo "Using reference annotation: $genedb_arg"
else
	echo "No annotation provided. Running IsoQuant in reference-free mode."
fi

#platform settings for this array task
if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
	TARGET_PLATFORM="pacbio"; ISOQUANT_TYPE="pacbio_ccs"; OUT_DIR="pacbio"
	RUNS=$(awk -F'\t' 'tolower($9) ~ /pacbio/ {print $1}' "$select_tsv")
else
	TARGET_PLATFORM="nanopore"; ISOQUANT_TYPE="nanopore"; OUT_DIR="nanopore"
	RUNS=$(awk -F'\t' 'tolower($9) ~ /nanopore/ || tolower($9) ~ /ont/ {print $1}' "$select_tsv")
fi

if [ -z "$RUNS" ]; then
	echo "No $TARGET_PLATFORM runs selected for $species_name. Nothing to do."
	exit 0
fi

#collect the fastq files for selected runs
FASTQ_ARGS=""
for run in $RUNS; do
	for fq in $data/fastq/*${run}*fastq.gz; do
		[ -f "$fq" ] && FASTQ_ARGS="$FASTQ_ARGS --fastq $fq"
	done
done

if [ -z "$FASTQ_ARGS" ]; then
	echo "ERROR: runs selected ($RUNS) but no matching fastq in $data/fastq/. Exiting."
	exit 1
fi

#reference genome (ignore .fai index)
reference=$(ls $data/fasta/*.fa* 2>/dev/null | grep -vE '\.fai$' | head -1)
if [ -z "$reference" ]; then
	echo "ERROR: no reference fasta in $data/fasta/. Exiting."
	exit 1
fi

echo "Running IsoQuant ($TARGET_PLATFORM) ..."
isoquant --threads "$SLURM_CPUS_PER_TASK" \
	--reference "$reference" \
	$genedb_arg \
	$FASTQ_ARGS \
	--data_type "$ISOQUANT_TYPE" \
	-o "$species_name/output/$OUT_DIR"

#peak mem
cgroup_dir=$(awk -F: '{print $NF}' /proc/self/cgroup)
if [ -f "/sys/fs/cgroup$cgroup_dir/memory.peak" ]; then
	peak_mem=$(cat "/sys/fs/cgroup$cgroup_dir/memory.peak")
	peak_mem_mb=$(awk "BEGIN {printf \"%.2f\", $peak_mem / 1048576}")
	echo ">Peak memory was $peak_mem_mb MegaBytes"
fi

elapsed_time=$(( $(date +%s) - start_time ))
echo "It takes $((elapsed_time / 60)) minutes"
echo ">ENDING at $(date)"
