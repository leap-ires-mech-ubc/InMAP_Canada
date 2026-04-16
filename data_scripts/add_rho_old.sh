#!/bin/bash
#SBATCH --time=23:59:00
#SBATCH --account=def-agiang01 #def-rscholes
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=4000M # --mem=250G
#SBATCH --array=2-9
#SBATCH --job-name='add_rho_and_compute_base'
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -euo pipefail
trap 'status=$?; echo "ERROR at line ${LINENO}: ${BASH_COMMAND}" >&2; exit $status' ERR
if [[ "${DEBUG:-0}" == "1" ]]; then set -x; fi

# ===== Modules =====
module load nco cdo   # If GNU parallel is a module on your cluster: module load parallel

# ===== Environment =====
export OMP_NUM_THREADS=4 #${SLURM_CPUS_PER_TASK:-8}
export CDO_PCTL_NTHREADS=$OMP_NUM_THREADS
export HDF5_DISABLE_VERSION_CHECK=1
export HDF5_USE_FILE_LOCKING=FALSE


# ===== Inputs =====
PREFIX_MAP=/home/tfmrodge/projects/def-rscholes/tfmrodge/InMAP_Canada/data_scripts/config_gemmachouts.txt
ARRAY_ID=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is not set}
prefix=$(awk -v ArrayTaskID="$ARRAY_ID" '$1==ArrayTaskID {print $2}' "$PREFIX_MAP")
if [[ -z "${prefix}" ]]; then
  echo "ERROR: No prefix found for array id ${ARRAY_ID} in ${PREFIX_MAP}" >&2
  exit 1
fi

in_dir="/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/${prefix}"
rhopath="/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/BASEGM_surface"
out_dir="${in_dir}/base_vars"          # final outputs (BASE* only)
mkdir -p "$out_dir"
tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

# ===== Config =====
# Scenario vars you need to compute BASE* (only these will be pulled from the scenario file)
export SCEN_VARS="${SCEN_VARS:-AF,TNI1,TNO3,THN3,TAM1,TSU1,TOC1}"
# RHO var name if known; will auto-detect if not found
export RHO_VARNAME="${RHO_VARNAME:-RHO}"
# Derived variables to keep at the end
export VARS_KEEP="${VARS_KEEP:-BASEPM25,BASEPNO3,BASEPNH4,BASEPSO4,BASESOA,BASEPRIM25}"

# Expressions to compute (order matters)
EXPR_FORMULA="$(cat <<'EXPR'
BASEPM25=AF;
BASEPNO3=(TNI1+TNO3+THN3)*RHO;
BASEPNH4=TAM1*RHO;
BASEPSO4=TSU1*RHO;
BASESOA=TOC1*RHO;
BASEPRIM25=BASEPM25-BASESOA-BASEPSO4-BASEPNH4-BASEPNO3;
EXPR
)"


# ===== Special scenarios that require leap-day averaging =====
SPECIAL_SCENARIOS=("BAU_2020_105" "COVID_2020_003")

is_special_scenario() {
  local p="$1"
  for s in "${SPECIAL_SCENARIOS[@]}"; do
    if [[ "$p" == "$s" ]]; then return 0; fi
  done
  return 1
}

if is_special_scenario "$prefix"; then
  export IS_SPECIAL="1"
else
  export IS_SPECIAL="0"
fi

echo "Building BASE variables from ${prefix} using RHO from BASEGM_surface ..."

# ===== Function (per file) =====
merge_and_compute_base() {
  local f="$1"
  local bn; bn=$(basename "$f")
  local out_file="${out_dir}/${bn}"

  # Parse basename pattern: YYYYMMDD00_0HH.nc
  local ts="${bn%%_*}"           # e.g., 2020111700
  local suffix="${bn#${ts}}"     # e.g., _020.nc
  local yyyy="${ts:0:4}"
  local mm="${ts:4:2}"
  local dd="${ts:6:2}"

  # ---- Map to 2019 RHO (reuse suffix _0HH.nc exactly)
  local rho_bn="2019${mm}${dd}00${suffix}"
  local rho_file="${rhopath}/${rho_bn}"

  # ---- Special scenarios: synthesize leap day 2020-02-29 from 2019-02-28 and 2019-03-01
  if [[ "$IS_SPECIAL" == "1" && "${yyyy}${mm}${dd}" == "20200229" ]]; then
    local prev_rho="${rhopath}/2019022800${suffix}"
    local next_rho="${rhopath}/2019030100${suffix}"
    if [[ ! -s "$prev_rho" || ! -s "$next_rho" ]]; then
      echo "ERROR: Missing RHO inputs for leap-day average: prev=$prev_rho next=$next_rho" >&2
      return 0   # skip gracefully to not kill batch
    fi
    local tmp_rho_avg="${tmpdir}/RHO_2020022900${suffix}"
    cdo -L -O ensmean "$prev_rho" "$next_rho" "$tmp_rho_avg"
    rho_file="$tmp_rho_avg"
  fi

  # ---- Check RHO exists
  if [[ ! -s "$rho_file" ]]; then
    echo "Skipped (no RHO): $bn (expected $(basename "$rho_file"))"
    return 0
  fi

  # ---- Detect RHO variable name if needed
  local rho_var="$RHO_VARNAME"
  if ! (cdo -s showname "$rho_file" || true) | tr ' ' '\n' | grep -Fxq "$rho_var"; then
    for v in RHO rho RHO_SFC air_density; do
      if (cdo -s showname "$rho_file" || true) | tr ' ' '\n' | grep -Fxq "$v"; then
        rho_var="$v"; break
      fi
    done
  fi
  if ! (cdo -s showname "$rho_file" || true) | tr ' ' '\n' | grep -Fxq "$rho_var"; then
    echo "ERROR: Could not find RHO-like variable in $(basename "$rho_file")" >&2
    return 0
  fi

  # ---- Ensure required scenario vars are present
  local present
  present=$((cdo -s showname "$f" || true) | tr ' ' '\n' \
             | grep -F -x -e AF -e TNI1 -e TNO3 -e THN3 -e TAM1 -e TSU1 -e TOC1 || true)
  local missing_list=()
  for need in AF TNI1 TNO3 THN3 TAM1 TSU1 TOC1; do
    if ! grep -Fxq "$need" <<<"$present"; then
      missing_list+=("$need")
    fi
  done
  if (( ${#missing_list[@]} > 0 )); then
    echo "Skipped (missing scenario vars ${missing_list[*]}): $bn"
    return 0
  fi
  local scen_csv
  scen_csv=$(echo "$present" | paste -sd, -)

  # ---- Optional: harmonize RHO vertical dim name to target if singleton
  # Detect likely vertical dims
  local target_zdim rho_zdim
  target_zdim=$(ncdump -h "$f"        | grep -E 'dimensions:' -A999 \
                  | grep -Eo '^[[:space:]]*(lev|level|level_1|height|pressure)[[:space:]]*=' \
                  | sed -E 's/^[[:space:]]*//;s/[[:space:]]*=.*$//' | head -n1 || true)
  rho_zdim=$(   ncdump -h "$rho_file" | grep -E 'dimensions:' -A999 \
                  | grep -Eo '^[[:space:]]*(lev|level|level_1|height|pressure)[[:space:]]*=' \
                  | sed -E 's/^[[:space:]]*//;s/[[:space:]]*=.*$//' | head -n1 || true)

  local tmp_rho="$rho_file"
  if [[ -n "$target_zdim" && -n "$rho_zdim" && "$target_zdim" != "$rho_zdim" ]]; then
    # get length of RHO vertical dim; only rename if singleton (surface)
    local rho_zlen
    rho_zlen=$(ncdump -h "$rho_file" | awk -v d="$rho_zdim" '
      $0 ~ /dimensions:/ {inD=1}
      inD && $0 ~ ("^\\s*" d "\\s*=\\s*[0-9]+\\s*;") {
        sub(/.*=/,""); sub(/;.*/,""); gsub(/[[:space:]]/,""); print; exit }' || true)
    if [[ -n "$rho_zlen" && "$rho_zlen" -eq 1 ]]; then
      local rho_fixed="${tmpdir}/rho_renamed_${bn}"
      if ncdump -h "$rho_file" | grep -qE "variables:.*[[:space:]]${rho_zdim}([[:space:]:]|$)"; then
        ncrename -d "$rho_zdim","$target_zdim" -v "$rho_zdim","$target_zdim" "$rho_file" "$rho_fixed"
      else
        ncrename -d "$rho_zdim","$target_zdim" "$rho_file" "$rho_fixed"
      fi
      tmp_rho="$rho_fixed"
    fi
    # If RHO has multiple levels, we leave dim names as-is (no blind remap).
  fi

  # ---- Build a temporary merged working file with only needed vars
  local work_nc="${tmpdir}/work_${bn}"
  cdo -L -O merge -selname,"$scen_csv" "$f" -selname,"$rho_var" "$tmp_rho" "$work_nc"

  # ---- Compute BASE* (per file), then keep only those outputs
  local tmp_expr="${tmpdir}/expr_${bn}"
  cdo -L -O expr,"$EXPR_FORMULA" "$work_nc" "$tmp_expr"
  cdo -L -O selname,"$VARS_KEEP" "$tmp_expr" "$out_file"

  echo "Done:  $bn  -> ${out_file##*/}   (RHO=${rho_var} from $(basename "$rho_file"))"
}

export -f merge_and_compute_base
export rhopath out_dir tmpdir IS_SPECIAL SCEN_VARS RHO_VARNAME VARS_KEEP EXPR_FORMULA

# ===== Run in parallel =====
# Robust to funny filenames; for debugging, you can add:  --joblog "${out_dir}/parallel_joblog.tsv"
find "$in_dir" -maxdepth 1 -type f -name '*.nc' -print0 \
| parallel -0 -j "$((SLURM_CPUS_PER_TASK / OMP_NUM_THREADS))" merge_and_compute_base {}

