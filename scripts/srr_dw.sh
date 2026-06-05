#!/bin/bash
# scripts/srr_dw.sh  –  SLURM array job: download one SRR/ERR per task from ENA
#
# Called by isoquant_prepare.sh as:
#   sbatch --array=0-N scripts/srr_dw.sh <species_name>
#
# Reads accession IDs from <species_name>/srr_list.tsv (one per line).
# Resolves the FASTQ URL from ENA's filereport API at download time.
# Saves to <species_name>/data/fastq/<accession>.fastq.gz

#SBATCH --output=logs/%x.%A_%a.out
#SBATCH --error=logs/%x.%A_%a.err
#SBATCH --job-name=srr_dw
#SBATCH --time=180
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2

start_time=$(date +%s)
echo ">STARTING at $(date)"

set -euo pipefail

SPECIES="$1"
FASTQ_DIR="${SPECIES}/data/fastq"
SRR_LIST="${SPECIES}/srr_list.tsv"

mkdir -p "$FASTQ_DIR"

#accession for this array task
ACCESSION=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SRR_LIST")

if [ -z "$ACCESSION" ]; then
	echo "ERROR: no accession at index $SLURM_ARRAY_TASK_ID in $SRR_LIST"
	exit 1
fi

OUT="${FASTQ_DIR}/${ACCESSION}.fastq.gz"

if [ -f "$OUT" ]; then
	echo "Already exists: $OUT"
	exit 0
fi

echo "Downloading $ACCESSION"

#ENA for the FTP URL at runtime (ENA returns paths without ftp:// prefix)
FTP_URL=$(curl -sf "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACCESSION}&result=read_run&fields=fastq_ftp" \
	| tail -n +2 | cut -f2 | tr ';' '\n' | head -1)

if [ -z "$FTP_URL" ]; then
	echo "ERROR: could not retrieve FTP URL for $ACCESSION from ENA"
	exit 1
fi

wget -q --tries=5 --timeout=120 -O "$OUT" "https://${FTP_URL}"
echo "Done: $OUT"

# Peak memory
cgroup_dir=$(awk -F: '{print $NF}' /proc/self/cgroup)
if [ -f "/sys/fs/cgroup$cgroup_dir/memory.peak" ]; then
	peak_mem=$(cat "/sys/fs/cgroup$cgroup_dir/memory.peak")
	peak_mem_mb=$(awk "BEGIN {printf \"%.2f\", $peak_mem / 1048576}")
	echo ">Peak memory was $peak_mem_mb MegaBytes"
fi

elapsed_time=$(( $(date +%s) - start_time ))
echo "It takes $((elapsed_time / 60)) minutes"
echo ">ENDING at $(date)"
