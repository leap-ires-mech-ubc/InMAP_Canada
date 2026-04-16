#!/bin/bash
#Select first level from GEMMACHDATA 
inpath="/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/BASEGM_2015_017/"
outpath="/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/BASEGM_surface/"
output_file="/home/tfmrodge/projects/def-rscholes/tfmrodge/missing_or_corrupt_files.txt"
level_value="0.998749971"
#for f in "$inpath"/20190808*.nc "$inpath"/20190809*.nc; do /home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/BASEGM_2015_017/2019120800_007.nc
# for f in "$inpath"/2019080800_015.nc; do
# while IFS= read -r filename; do 
#     fullpath="${inpath}${filename}"
#     bn=$(basename "$filename")
#     out_file="${outpath}${bn%.nc}.nc"

#     if [[ -s "$fullpath" ]]; then
#         echo "Processing $bn → $out_file"
#         cdo sellevel,"$level_value" "$fullpath" "$out_file"
#     else
#         echo "Skipping $bn: file is empty or missing."
#     fi
# done < "$output_file"
# done
#Also check if anything else is missin
# Loop through all .nc files
# outpath="/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/BASEGM_surface/"
# for f in "$outpath"/*.nc; do
#   if [[ ! -s "$f" ]]; then
#     echo "$(basename "$f")" >> "$output_file"
#   fi
# done
#Also check if anything else is missin
set -euo pipefail

outpath="/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/BASEGM_surface/"
output_file="/home/tfmrodge/projects/def-rscholes/tfmrodge/missing_or_corrupt_files.txt"

# Variable to check for; override via env CHECK_VAR if needed (e.g., CHECK_VAR=air_density)
CHECK_VAR="${CHECK_VAR:-RHO}"
# Accept common aliases; extend if your files use different names
VAR_ALIASES="${VAR_ALIASES:-RHO rho RHO_SFC air_density}"

# Ensure cdo is available
module load cdo >/dev/null 2>&1 || true
command -v cdo >/dev/null 2>&1 || { echo "ERROR: cdo not found on PATH."; exit 1; }

mkdir -p "$(dirname "$output_file")"
: > "$output_file"  # truncate

log_bad() {
  # status, file, reason
  printf "%s\t%s\t%s\n" "$1" "$2" "$3" >> "$output_file"
}

# Scan all .nc files in outpath
while IFS= read -r -d '' f; do
  bn=$(basename "$f")

  # 1) Empty file?
  if [[ ! -s "$f" ]]; then
    log_bad "missing_or_empty" "$f" "size=0 or file absent"
    continue
  fi

  # 2) Can cdo open it?
  if ! cdo -s showname "$f" >/dev/null 2>&1; then
    # Example error: "Unsupported file structure"
    log_bad "corrupt_or_unreadable" "$f" "cdo showname failed"
    continue
  fi

  # 3) Does it contain the required variable (RHO or alias)?
  #    Use tokens split on spaces; accept group prefixes (e.g., 'group/RHO')
  show=$(cdo -s showname "$f" | tr ' ' '\n')
  has_var=0
  for v in $VAR_ALIASES; do
    if grep -Eq "(^|/)${v}$" <<<"$show"; then
      has_var=1; break
    fi
  done
  if (( has_var == 0 )); then
    log_bad "var_missing" "$f" "none of (${VAR_ALIASES}) found"
    continue
  fi

  # If we reach here, file looks good; do not log it (you asked for bad-only)
done < <(find "$outpath" -maxdepth 1 -type f -name '*.nc' -print0)

echo "Bad files report written to: $output_file"


