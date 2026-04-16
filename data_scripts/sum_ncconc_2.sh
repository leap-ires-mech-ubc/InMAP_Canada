#!/usr/bin/env bash
#SBATCH --time=23:59:00
#SBATCH --account=def-rscholes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=250G
#SBATCH --array=2
#SBATCH --job-name='annual_avg_filtered'
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
############################################
# Get annual average concentrations from GEMMACH outputs. Inputs should already be pre-processed
# to have all necessary variables
############################################

set -euo pipefail
trap 'status=$?; echo "ERROR at line ${LINENO}: ${BASH_COMMAND}" >&2; exit $status' ERR

############################################
# Helper: run CDO, print which file/op failed on error
############################################
cdo_or_die() {
  local op_desc="$1"   # e.g., "yearmean+mergetime"
  local infile="$2"    # e.g., path/glob or a label for stage
  shift 2
  if ! cdo "$@"; then
    echo "ERROR: CDO failed during ${op_desc} on: ${infile}" >&2
    exit 2
  fi
}

############################################
# Allow SLURM or CLI use
############################################
ARRAY_ID="${SLURM_ARRAY_TASK_ID:-${1:-}}"
if [[ -z "$ARRAY_ID" ]]; then
  echo "Usage: $0 <id> (or submit via sbatch --array)" >&2
  exit 1
fi

############################################
# Module loading (safe outside SLURM)
############################################
if command -v module >/dev/null 2>&1; then
  module load cdo
fi
if ! command -v cdo >/dev/null 2>&1; then
  echo "ERROR: cdo not found on PATH (load module or activate environment)" >&2
  exit 1
fi
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}

############################################
# Prefix lookup
############################################
PREFIX_MAP=/home/tfmrodge/projects/def-rscholes/tfmrodge/InMAP_Canada/data_scripts/config_gemmachouts.txt
prefix=$(awk -v id="$ARRAY_ID" '$1==id {print $2}' "$PREFIX_MAP")
if [[ -z "$prefix" ]]; then
  echo "ERROR: No prefix for id $ARRAY_ID in $PREFIX_MAP" >&2
  exit 1
fi

############################################
# Input/output paths
############################################
# Preprocessor should have already calculated the correct variables in each file
in_dir="/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/${prefix}/20260129_base_vars"
out_dir="/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/20260129_scenario_sums"
mkdir -p "$out_dir"

stamp=$(date +%Y%m%d)
native_out="${out_dir}/${stamp}_${prefix}_annual_mean.nc"
lcc_out="${out_dir}/${stamp}_${prefix}_annual_mean_LCC.nc"
LCC_GRID="/home/tfmrodge/projects/def-rscholes/tfmrodge/InMAP_Canada/data_scripts/gridfile_lcc.txt"

echo "Computing annual mean for ${prefix}"
echo "Input dir   : $in_dir"
echo "Output dir  : $out_dir"
echo "Threads     : $OMP_NUM_THREADS"

############################################
# Gather and sanity-check files
############################################
shopt -s nullglob
mapfile -t files < <(printf "%s\n" "$in_dir"/*.nc | sort)
if [[ ${#files[@]} -eq 0 ]]; then
  echo "ERROR: No .nc files found in $in_dir" >&2
  exit 1
fi
echo "Files found: ${#files[@]} (expect ~8760)"

############################################
# Merge + yearmean
############################################
# Use the validated array to avoid empty globs or too literal arguments.
echo "Merging & computing annual mean…"
# CDO expects: cdo [opts] yearmean -mergetime f1 f2 ... native_out
cdo_or_die "yearmean+mergetime" "$in_dir/*.nc" \
  -L -P "$OMP_NUM_THREADS" -O timmean -mergetime "${files[@]}" "$native_out"

############################################
# Reproject
############################################
echo "Reprojecting to LCC…"
cdo_or_die "remapbil(LCC)" "$native_out" \
  -L -P "$OMP_NUM_THREADS" -O remapbil,"$LCC_GRID" "$native_out" "$lcc_out"

############################################
# Optional: quick variable check (names only)
############################################
echo "Output variables:"
cdo -s showname "$native_out" | tr ' ' '\n' | sort -u

echo "Done!"
echo "Native annual mean: $native_out"
echo "LCC annual mean   : $lcc_out"