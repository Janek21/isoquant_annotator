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

species_name="$1"

out_dir="$species_name/data/fastq"
mkdir -p "$out_dir"

# accession for this array task (line SLURM_ARRAY_TASK_ID + 1)
accession=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$species_name/srr_list.tsv")

if [ -z "$accession" ]; then
	echo "ERROR: no accession at index $SLURM_ARRAY_TASK_ID in $species_name/srr_list.tsv"
	exit 1
fi

result_file="$out_dir/$accession.fastq.gz"

if [ -s "$result_file" ]; then #check if file has already been downloaded
	echo "Skipping $accession: $result_file already exists."
else
	echo "--- Processing $accession ---"

	#Resolve the FASTQ URL from ENA (returns paths without a scheme prefix)
	ftp_url=$(curl -sf "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${accession}&result=read_run&fields=fastq_ftp" | tail -n +2 | cut -f2 | tr ';' '\n' | head -1)

	if [ -z "$ftp_url" ]; then
		echo "Error: could not retrieve FTP URL for $accession from ENA."
		exit 1
	fi

	#Download to a temp file first (avoids leaving a partial file that the skip-check above would mistake for a finished download)
	tmp_file="$out_dir/$accession.fastq.gz.part"
	wget -q --tries=5 --waitretry=60 --random-wait --timeout=120 -O "$tmp_file" "https://${ftp_url}"

	#Validate gzip integrity so a truncated download fails loudly here instead of leaving a corrupt result file behind.
	zcat "$tmp_file" > /dev/null
	mv "$tmp_file" "$result_file"
	echo "Complete: new file is $result_file"
fi

# Record memory usage at the end
cgroup_dir=$(awk -F: '{print $NF}' /proc/self/cgroup)
# Check if the path exists to avoid errors on different cgroup versions
if [ -f "/sys/fs/cgroup$cgroup_dir/memory.peak" ]; then
	peak_mem=$(cat "/sys/fs/cgroup$cgroup_dir/memory.peak")
	peak_mem_mb=$(awk "BEGIN {printf \"%.2f\", $peak_mem / 1048576}")
	echo ">Peak memory was $peak_mem_mb MegaBytes"
fi

# Record end
elapsed_time=$(( $(date +%s) - start_time ))
echo "It takes $((elapsed_time / 60)) minutes"
echo ">ENDING at $(date)"
