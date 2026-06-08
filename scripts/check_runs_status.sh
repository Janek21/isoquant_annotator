#!/bin/bash
# Classify every species listed in dataspecie.txt by how far it got in the
# IsoQuant pipeline.
# Buckets (mutually exclusive):
#   NO FOLDER   - species in the list with no working dir (never started)
#   NOT RUN     - folder exists but IsoQuant produced no transcript models yet
#   NOT EVAL    - IsoQuant ran but the annotation was not merged/evaluated
#   DONE        - merged annotation exists, either as the canonical copy in
#                 summary/pred/ or the per-species output/eval/ copy. Detected
#                 even if the working dir was removed after migration, since
#                 evaluation.sh relocates the prediction to summary/pred/.
#
# Unlike LyRic, isoquant_execute.sh uses the species name verbatim as the
# working-dir name (no sanitization), so each raw line from dataspecie.txt is
# already the folder name (e.g. "Entamoeba_histolytica_HM-1:IMSS").
#
# Usage (from the repo root):
#   bash scripts/check_runs_status.sh [dataspecie.txt] [base_dir]
# Defaults:
#   dataspecie.txt = /home/jj/Desktop/Data_science/CRG/TFM2/projects/busco_references/dataspecie.txt
#   base_dir       = .  (where the species working dirs and summary/ live)

shopt -s nullglob

datafile="${1:-/home/jj/Desktop/Data_science/CRG/TFM2/projects/busco_references/dataspecie.txt}"
base_dir="${2:-.}"
summary_pred="$base_dir/summary/pred"

if [ ! -s "$datafile" ]; then
	echo "Species list not found or empty: $datafile" >&2
	exit 1
fi

#--- stage detection (isoquant layout) ---

has_run() {  #IsoQuant produced a non-empty transcript_models.gtf (pacbio and/or nanopore, nested)
	local m
	m=$(find "$1/output" -name '*transcript_models.gtf' -size +0c 2>/dev/null | head -1)
	[ -n "$m" ]
}

has_eval() {  #merged annotation exists: canonical summary/pred copy, or per-species output/eval copy
	local name="$1"
	local pred=("$summary_pred/${name}"_[0-9]*_pred.gff)
	[ -s "${pred[0]:-}" ] && return 0
	local merged=("$base_dir/$name"/output/eval/merged_*.gff)
	[ -s "${merged[0]:-}" ]
}

#--- classify ---

no_folder=()
not_run=()
not_eval=()
done_sp=()

while IFS= read -r raw || [ -n "$raw" ]; do
	raw="${raw%$'\r'}"                 # strip stray CR
	[ -z "${raw// }" ] && continue     # skip blank lines
	work_name="$raw"                   # isoquant uses the name verbatim as the dir
	sp_dir="$base_dir/$work_name"

	if has_eval "$work_name"; then     # DONE first: working dir may be gone post-migration
		done_sp+=("$raw")
	elif [ ! -d "$sp_dir" ]; then
		no_folder+=("$raw")
	elif has_run "$sp_dir"; then
		not_eval+=("$raw")
	else
		not_run+=("$raw")
	fi
done < "$datafile"

#--- report ---

print_group() {
	local title="$1"; shift
	printf '\n== %s (%d) ==\n' "$title" "$#"
	if [ "$#" -eq 0 ]; then
		echo "  (none)"
	else
		printf '  %s\n' "$@"
	fi
}

total=$(( ${#no_folder[@]} + ${#not_run[@]} + ${#not_eval[@]} + ${#done_sp[@]} ))
echo "Species list: $datafile"
echo "Base dir:     $base_dir"
echo "Total species: $total"

print_group "NO FOLDER (never started)"   "${no_folder[@]}"
print_group "FOLDER, ISOQUANT NOT RUN"    "${not_run[@]}"
print_group "RUN but NOT EVALUATED"       "${not_eval[@]}"
print_group "DONE (merged annotation)"    "${done_sp[@]}"
