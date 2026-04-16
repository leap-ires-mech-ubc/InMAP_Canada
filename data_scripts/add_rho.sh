#!/bin/bash
#SBATCH --time=23:59:00
#SBATCH --account=def-agiang01 #def-rscholes
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --array=2-9
#SBATCH --job-name='add_rho_and_compute_base'
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -euo pipefail

# --- Allow local testing: if no SLURM_ARRAY_TASK_ID, use first argument ---
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  if [[ $# -lt 1 ]]; then
    echo "Usage: sh add_rho.sh <task_id>" >&2
    exit 1
  fi
  SLURM_ARRAY_TASK_ID="$1"
  SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-16}"
  echo "Running locally with SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"
else
  echo "Running under SLURM with SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"
fi

trap 'status=$?; echo "ERROR at line ${LINENO}: ${BASH_COMMAND}" >&2; exit $status' ERR
if [[ "${DEBUG:-0}" == "1" ]]; then set -x; fi

# ===== Modules =====
module load nco cdo  # (Add: module load parallel  if needed on your cluster)

# ===== Environment =====
export OMP_NUM_THREADS=1
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

# Paths
in_dir="/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/${prefix}"
rhopath="/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/BASEGM_surface"
out_dir="${in_dir}/base_vars"          # final outputs
mkdir -p "$out_dir"

# TMPDIR: use SLURM_TMPDIR if present, else mktemp; clean only if mktemp
tmpdir="${SLURM_TMPDIR:-$(mktemp -d)}"
trap '
  if [[ -n "${SLURM_TMPDIR:-}" && "$tmpdir" == "$SLURM_TMPDIR" ]]; then
    :  # SLURM will clean local scratch
  else
    rm -rf "$tmpdir"
  fi
' EXIT

# ===== REDO_MODE (manual setting) =====
# 0 = full run (process ALL .nc files)
# 1 = redo only failures (requires TSV summary)
REDO_MODE=1    # <<< YOU SET THIS

CHECKLOG_ROOT="/home/tfmrodge/projects/def-rscholes/tfmrodge/InMAP_Canada/checklogs"
TSV_LATEST="${CHECKLOG_ROOT}/${prefix}/coverage_${prefix}_latest.tsv"

FILES_TO_REDO=()
if (( REDO_MODE == 1 )); then
  echo "Redo mode active."
  if [[ ! -s "$TSV_LATEST" ]]; then
    echo "ERROR: REDO_MODE=1 but TSV summary not found:"
    echo "  $TSV_LATEST"
    exit 1
  fi
  echo "Using TSV summary:"
  echo "  $TSV_LATEST"

  # Read only not-ok rows; first column (basename) sanitized
  while IFS=$'\t' read -r bn status _; do
    [[ -z "$bn" ]] && continue
    if [[ "$status" != "ok" ]]; then
      # sanitize basename
      bn="${bn#"${bn%%[![:space:]]*}"}"   # ltrim
      bn="${bn%"${bn##*[![:space:]]}"}"   # rtrim
      bn="${bn##*/}"
      [[ "$bn" = /* ]] && bn="${bn#/}"
      FILES_TO_REDO+=("$bn")
    fi
  done < <(awk -F'\t' 'NR>1 {print $1 "\t" $2}' "$TSV_LATEST")

  if (( ${#FILES_TO_REDO[@]} == 0 )); then
    echo "Summary shows 0 failed files — nothing to redo. Exiting."
    exit 0
  fi

  echo "Redoing ${#FILES_TO_REDO[@]} failed files:"
  printf '   %s\n' "${FILES_TO_REDO[@]}"
else
  echo "Full run mode active — will process ALL .nc files in $in_dir"
fi

# Sanity: in_dir should exist
if [[ ! -d "$in_dir" ]]; then
  echo "ERROR: in_dir not found: $in_dir" >&2
  exit 1
fi

# ===== Per-run status directories (safe under parallel) =====
status_dir="${out_dir}/_run_${ARRAY_ID}"
fail_dir="${status_dir}/failed"
succ_dir="${status_dir}/succeeded"
#Clear any markers from previous runs
rm -rf "$status_dir"
mkdir -p "$fail_dir" "$succ_dir"
summary_file="${status_dir}/summary.txt"

# ===== Config =====
export SCEN_VARS="${SCEN_VARS:-AF,TNI1,TNO3,THN3,TAM1,TSU1,TOC1}"
export RHO_VARNAME="${RHO_VARNAME:-RHO}"
export VARS_KEEP="${VARS_KEEP:-BASEPM25,BASEPNO3,BASEPNH4,BASEPSO4,BASESOA,BASEPRIM25}"

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

# ===== Utility to record failures/successes =====
record_failure() {
  # Args: category, file_basename, reason
  local category="$1"; local bn="$2"; local reason="$3"
  local mark="${fail_dir}/${category}__${bn}.txt"
  printf "FILE=%s\nCATEGORY=%s\nREASON=%s\n" "$bn" "$category" "$reason" > "$mark"
}
record_success() {
  local bn="$1"
  : > "${succ_dir}/${bn}.ok"
}

# ===== Per-file worker (accepts FULL PATH to scenario file) =====
merge_and_compute_base() {
  local f="$1"
  local bn
  bn="$(basename "$f")"
  local out_file="${out_dir}/${bn}"

  # Optional: uncomment for debugging path resolution
  # echo "DEBUG: bn='$bn' f='$f'" >&2

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
      record_failure "rho_leap_missing" "$bn" "prev=$prev_rho next=$next_rho"
      echo "ERROR: Missing RHO inputs for leap-day average: prev=$prev_rho next=$next_rho" >&2
      return 0
    fi
    local tmp_rho_avg="${tmpdir}/RHO_2020022900${suffix}"
    if ! cdo -L -O ensmean "$prev_rho" "$next_rho" "$tmp_rho_avg"; then
      record_failure "rho_leap_ensmean_fail" "$bn" "cdo ensmean failed"
      echo "ERROR: cdo ensmean failed for leap-day average" >&2
      return 0
    fi
    rho_file="$tmp_rho_avg"
  fi

  # ---- Check RHO exists
  if [[ ! -s "$rho_file" ]]; then
    record_failure "rho_missing" "$bn" "expected $(basename "$rho_file")"
    echo "Skipped (no RHO): $bn (expected $(basename "$rho_file"))" >&2
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
    record_failure "rho_var_not_found" "$bn" "in $(basename "$rho_file")"
    echo "ERROR: Could not find RHO-like variable in $(basename "$rho_file"))" >&2
    return 0
  fi

  # ---- Ensure required scenario vars are present
  local present
  present=$( (cdo -s showname "$f" || true) | tr ' ' '\n' \
             | grep -F -x -e AF -e TNI1 -e TNO3 -e THN3 -e TAM1 -e TSU1 -e TOC1 || true )
  local missing_list=()
  for need in AF TNI1 TNO3 THN3 TAM1 TSU1 TOC1; do
    if ! grep -Fxq "$need" <<<"$present"; then
      missing_list+=("$need")
    fi
  done
  if (( ${#missing_list[@]} > 0 )); then
    record_failure "scenario_vars_missing" "$bn" "${missing_list[*]}"
    echo "Skipped (missing scenario vars ${missing_list[*]}): $bn" >&2
    return 0
  fi
  local scen_csv
  scen_csv=$(echo "$present" | paste -sd, -)

  # ---- Optional: harmonize RHO vertical dim name to target if singleton
  local target_zdim rho_zdim
  target_zdim=$(ncdump -h "$f"        | grep -E 'dimensions:' -A999 \
                  | grep -Eo '^[[:space:]]*(lev|level|level_1|height|pressure)[[:space:]]*=' \
                  | sed -E 's/^[[:space:]]*//;s/[[:space:]]*=.*$//' | head -n1 || true)
  rho_zdim=$(   ncdump -h "$rho_file" | grep -E 'dimensions:' -A999 \
                  | grep -Eo '^[[:space:]]*(lev|level|level_1|height|pressure)[[:space:]]*=' \
                  | sed -E 's/^[[:space:]]*//;s/[[:space:]]*=.*$//' | head -n1 || true)

  local tmp_rho="$rho_file"
  if [[ -n "$target_zdim" && -n "$rho_zdim" && "$target_zdim" != "$rho_zdim" ]]; then
    local rho_zlen
    rho_zlen=$(ncdump -h "$rho_file" | awk -v d="$rho_zdim" '
      $0 ~ /dimensions:/ {inD=1}
      inD && $0 ~ ("^\\s*" d "\\s*=\\s*[0-9]+\\s*;") {
        sub(/.*=/,""); sub(/;.*/,""); gsub(/[[:space:]]/,""); print; exit }' || true)
    if [[ -n "$rho_zlen" && "$rho_zlen" -eq 1 ]]; then
      local rho_fixed="${tmpdir}/rho_renamed_${bn}"
      if ncdump -h "$rho_file" | grep -qE "variables:.*[[:space:]]${rho_zdim}([[:space:]:]|$)"; then
        ncrename -O -d "$rho_zdim","$target_zdim" -v "$rho_zdim","$target_zdim" "$rho_file" "$rho_fixed"
      else
        ncrename -O -d "$rho_zdim","$target_zdim" "$rho_file" "$rho_fixed"
      fi
      tmp_rho="$rho_fixed"
    fi
  fi

  # ---- Build a temporary merged working file with only needed vars
  local work_nc="${tmpdir}/work_${bn}"
  if ! cdo -L -O merge -selname,"$scen_csv" "$f" -selname,"$rho_var" "$tmp_rho" "$work_nc"; then
    record_failure "cdo_merge_fail" "$bn" "cdo merge"
    echo "ERROR: cdo merge failed on $bn" >&2
    return 0
  fi

  # ---- Compute BASE* (per file), then keep only those outputs
  local tmp_expr="${tmpdir}/expr_${bn}"
  if ! cdo -L -O expr,"$EXPR_FORMULA" "$work_nc" "$tmp_expr"; then
    record_failure "cdo_expr_fail" "$bn" "cdo expr"
    echo "ERROR: cdo expr failed on $bn" >&2
    return 0
  fi
  if ! cdo -L -O selname,"$VARS_KEEP" "$tmp_expr" "$out_file"; then
    record_failure "cdo_selname_fail" "$bn" "cdo selname"
    echo "ERROR: cdo selname failed on $bn" >&2
    return 0
  fi

  record_success "$bn"
  echo "Done:  $bn  -> ${out_file##*/}   (RHO=${rho_var} from $(basename "$rho_file"))"
}

# Export functions and vars used by the worker
export -f merge_and_compute_base
export -f record_failure
export -f record_success
export rhopath out_dir tmpdir IS_SPECIAL RHO_VARNAME VARS_KEEP EXPR_FORMULA fail_dir succ_dir

# ===== Run in parallel (pass FULL PATHS to the worker) =====
JOBS=$(( ${SLURM_CPUS_PER_TASK:-16} / ${OMP_NUM_THREADS:-1} ))
(( JOBS <= 0 )) && JOBS=1

if (( REDO_MODE == 1 )); then
  # Prepend full path to each basename from TSV
  printf '%s\n' "${FILES_TO_REDO[@]}" \
  | sed "s|^|$in_dir/|" \
  | parallel -j "$JOBS" merge_and_compute_base {}
else
  # Full run: feed full paths directly
  find "$in_dir" -maxdepth 1 -type f -name '*.nc' -print0 \
  | parallel -0 -j "$JOBS" merge_and_compute_base {}
fi

# ===== End-of-run summary =====
total_inputs=$(find "$in_dir" -maxdepth 1 -type f -name '*.nc' | wc -l | awk '{print $1}')
succ_count=$(find "$succ_dir" -type f -name '*.ok' 2>/dev/null | wc -l | awk '{print $1}')
fail_count=$(find "$fail_dir" -type f -name '*.txt' 2>/dev/null | wc -l | awk '{print $1}')

{
  echo "===== SUMMARY for prefix ${prefix} (array ${ARRAY_ID}) ====="
  echo "Total input files:  ${total_inputs}"
  echo "Succeeded:          ${succ_count}"
  echo "Failed/Skipped:     ${fail_count}"
  echo
  if (( fail_count > 0 )); then
    echo "---- Failure breakdown by category ----"
    find "$fail_dir" -type f -name '*.txt' -printf '%f\n' \
      | awk -F'__' '{print $1}' \
      | sort | uniq -c | awk '{printf "%6s  %s\n", $1, $2}'
    echo
    echo "---- Failed files (category :: file :: reason) ----"
    while IFS= read -r f; do
      bn=$(awk -F= '/^FILE=/{print $2}' "$f")
      catg=$(awk -F= '/^CATEGORY=/{print $2}' "$f")
      rsn=$(awk -F= '/^REASON=/{print $2}' "$f")
      echo "${catg} :: ${bn} :: ${rsn}"
    done < <(find "$fail_dir" -type f -name '*.txt' | sort)
  else
    echo "No failures 🎉"
  fi
} | tee "$summary_file" >&2

echo "Detailed summary saved to: $summary_file"
