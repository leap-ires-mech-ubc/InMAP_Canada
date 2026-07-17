#!/usr/bin/env bash
#SBATCH --time=23:59:00
#SBATCH --account=def-agiang01 #def-rscholes
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --array=2,9
#SBATCH --job-name='add_rho_and_compute_base'
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -euo pipefail
# =========================
# Debug & utilities
# =========================
DEBUG="${DEBUG:-0}"
dbg(){ if (( DEBUG )); then echo "DEBUG: $*" >&2; fi; }

trap 'status=$?; echo "ERROR at line ${LINENO}: ${BASH_COMMAND}" >&2; exit $status' ERR

# =========================
# Allow local testing (no SLURM)
# =========================
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  [[ $# -lt 1 ]] && { echo "Usage: bash compute_INMAP_variables.sh <task_id>"; exit 1; }
  SLURM_ARRAY_TASK_ID="$1"
  SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-16}"
  echo "Running locally with SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"
else
  echo "Running under SLURM with SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"
fi

# =========================
# Modules & environment
# =========================
module load nco cdo || true

export OMP_NUM_THREADS=1
export CDO_PCTL_NTHREADS=$OMP_NUM_THREADS
export HDF5_DISABLE_VERSION_CHECK=1
export HDF5_USE_FILE_LOCKING=FALSE
export TMPDIR="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"

# Optional one-time debug about tools
if (( DEBUG )); then
  dbg "cdo path: $(command -v cdo || echo 'not found')"
  dbg "cdo version: $(cdo --version 2>&1 | head -n1)"
  dbg "ncdump path: $(command -v ncdump || echo 'not found')"
  dbg "ncks path: $(command -v ncks || echo 'not found')"
fi

# =========================
# Inputs / paths
# =========================
PREFIX_MAP=/home/tfmrodge/projects/def-rscholes/tfmrodge/InMAP_Canada/data_scripts/config_gemmachouts.txt
ARRAY_ID=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is not set}
prefix=$(awk -v ArrayTaskID="$ARRAY_ID" '$1==ArrayTaskID {print $2}' "$PREFIX_MAP")
[[ -z "$prefix" ]] && { echo "ERROR: No prefix for array id $ARRAY_ID in $PREFIX_MAP" >&2; exit 1; }

basegm_dir="/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/BASEGM_2015_017/"
in_dir="/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/${prefix}"

# BASE detection: treat any prefix containing "BASEGM" as BASE case
if [[ "$prefix" == *BASEGM* ]]; then
  in_dir="$basegm_dir"
  IS_BASEGM=1
else
  IS_BASEGM=0
fi

# Out dir
out_dir="${in_dir}/20260129_base_vars"
mkdir -p "$out_dir"
DEBUG_DIR="${out_dir}/_debug"; (( DEBUG )) && mkdir -p "$DEBUG_DIR"

# =========================
# Config
# =========================
SCEN_VARS="${SCEN_VARS:-AF,TNI1,TNO3,THN3,TAM1,TSU1,TOC1}"
RHO_VARNAME="${RHO_VARNAME:-RHO}"
VARS_KEEP="${VARS_KEEP:-RHO,BASEPM25,BASEPNO3,BASEPNH4,BASEPSO4,BASESOA,BASEPRIM25}"
SURFACE_LEVEL="${SURFACE_LEVEL:-0.99875}"

# Non-BASE expr formula (consistent with your current choice)
#OLD BASEPNO3=(TNI1+TNO3+THN3)*RHO;
EXPR_FORMULA=$(cat <<'EXPR'
RHO=RHO;
BASEPM25=AF;
BASEPNO3=(TNI1)*RHO;
BASEPNH4=TAM1*RHO;
BASEPSO4=TSU1*RHO;
BASESOA=TOC1*RHO;
BASEPRIM25=BASEPM25-BASESOA-BASEPSO4-BASEPNH4-BASEPNO3;
EXPR
)

# 2020 leap-day handling scenarios (borrowed RHO)
SPECIAL_SCENARIOS=("BAU_2020_105" "COVID_2020_003")
IS_SPECIAL=0; for s in "${SPECIAL_SCENARIOS[@]}"; do [[ "$prefix" == "$s" ]] && IS_SPECIAL=1; done
rhopath="/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/BASEGM_surface" #"$basegm_dir"

# =========================
# REDO switch
# =========================
REDO_MODE="${REDO_MODE:-1}"
CHECKLOG_ROOT="/home/tfmrodge/projects/def-rscholes/tfmrodge/InMAP_Canada/checklogs"
TSV_LATEST="${CHECKLOG_ROOT}/${prefix}/coverage_${prefix}_latest.tsv"

# =========================
# Per-run status / logging
# =========================
status_dir="${out_dir}/_run_${ARRAY_ID}"
fail_dir="${status_dir}/failed"
succ_dir="${status_dir}/succeeded"
rm -rf "$status_dir"; mkdir -p "$fail_dir" "$succ_dir"
summary_file="${status_dir}/summary.txt"

record_failure(){ # category, bn, reason
  printf "FILE=%s\nCATEGORY=%s\nREASON=%s\n" "$2" "$1" "$3" > "${fail_dir}/${1}_${2}.txt";
}
record_success(){ : > "${succ_dir}/$1.ok"; }

# =========================
# Year detection (for 2020 RHO borrowing)
# =========================
year_guess=$(grep -oE '[12][0-9]{3}' <<<"$prefix" | head -n1 || true)
if [[ -z "$year_guess" ]]; then
  f_any=$(ls -1 "${in_dir}"/*.nc 2>/dev/null | head -n1 || true)
  [[ -n "${f_any:-}" ]] && year_guess=$(basename "$f_any" | cut -c1-4)
fi
[[ -z "$year_guess" ]] && year_guess="2019"

# =========================
# Worker: per-file processing (FULL PATH argument)
# =========================
merge_and_compute_base() {
  local f="$1"
  local bn
  bn="$(basename "$f")"
  local out_file="${out_dir}/${bn}"

  dbg "Processing file: [$f]"
  [[ -n "${f:-}" ]] && ls -l "$f" >&2 || true

  if [[ -z "${f:-}" || ! -s "$f" ]]; then
    record_failure "missing_or_empty" "${bn:-empty}" "file missing or empty"
    echo "MISSING/EMPTY: ${bn:-empty}" >&2
    return 0
  fi

  # Find if variables exist
local present
present="$(
  cdo -s showname "$f" \
    | tr '[:space:]' '\n' \
    | awk 'NF>0'
)" || {
  record_failure "unreadable" "$bn" "cdo showname failed"
  echo "UNREADABLE: $bn (cdo showname)" >&2
  return 0
}

(( DEBUG )) && dbg "Vars($bn): $(echo "$present" | paste -sd, -)"vars


  if [[ -z "${present//[[:space:]]/}" ]]; then
    record_failure "unreadable_or_no_vars" "$bn" "No variables via ncdump -h"
    echo "UNREADABLE/EMPTY-VARS (ncdump): $bn" >&2
    return 0
  fi

  # ======================= BASEGM CASE =======================
  if (( IS_BASEGM == 1 )); then
    # Ensure key species exist (exact names)
    local missing=()
    for v in AF TNI1 TNO3 THN3 TAM1 TSU1 TOC1 RHO; do
      grep -Fxq "$v" <<<"$present" || missing+=("$v")
    done
    if (( ${#missing[@]} > 0 )); then
      record_failure "base_vars_missing" "$bn" "${missing[*]}"
      echo "MISSING VAR(s) in BASE file $bn: ${missing[*]}" >&2
      return 0
    fi

    # 1) Extract surface at sigma=0.99875 (chem surface)
    local sfc="${TMPDIR:-${SLURM_TMPDIR:-/tmp}}/sfc_${bn}"
    if ! cdo -L -O sellevel,"$SURFACE_LEVEL" "$f" "$sfc"; then
      record_failure "sellevel_fail" "$bn" "level=$SURFACE_LEVEL"
      echo "SELLEVEL FAIL: $bn" >&2
      return 0
    fi

    # 2) Compute all BASE* variables via ncap2 (2‑D → safe, fast)
    #Old BASEPNO3=(TNI1+TNO3+THN3)*RHO; # 
    local expr_nc="${TMPDIR:-${SLURM_TMPDIR:-/tmp}}/expr_${bn}"
    if ! ncap2 -O \
      -s "RHO=RHO;
          BASEPM25=AF;
          BASEPNO3=TNI1*RHO;
          BASEPNH4=TAM1*RHO;
          BASEPSO4=TSU1*RHO;
          BASESOA=TOC1*RHO;
          BASEPRIM25=BASEPM25-BASESOA-BASEPSO4-BASEPNH4-BASEPNO3;" \
      "$sfc" "$expr_nc"; then
      record_failure "ncap2_expr_fail" "$bn" "ncap2 expr"
      echo "NCAP2 FAIL: $bn" >&2
      return 0
    fi

        # Clip tiny negative BASEPRIM25 values
    if ! ncap2 -O -s 'where(BASEPRIM25<0) BASEPRIM25=0' "$expr_nc" "$expr_nc"; then
      record_failure "clip_prim25_fail" "$bn" "ncap2 clip"
      echo "CLIP FAIL: $bn" >&2
      return 0
    fi

    # 3) Keep only the desired output fields
    if ! cdo -L -O selname,"$VARS_KEEP" "$expr_nc" "$out_file"; then
      record_failure "selname_fail" "$bn" "$VARS_KEEP"
      echo "SELNAME FAIL: $bn" >&2
      return 0
    fi

    record_success "$bn"
    echo "BASE OK: $bn -> ${out_file##*/} (surface=$SURFACE_LEVEL via original method)"
    return 0
  fi
  # ===================== END BASEGM CASE =====================

  # ---------- Non-BASE scenarios ----------
  local need_miss=()
  for v in AF TNI1 TNO3 THN3 TAM1 TSU1 TOC1; do
    grep -Fxq "$v" <<<"$present" || need_miss+=("$v")
  done
  if (( ${#need_miss[@]} > 0 )); then
    record_failure "scenario_vars_missing" "$bn" "${need_miss[*]}"
    echo "MISSING VAR(s) in scenario $bn: ${need_miss[*]}" >&2
    return 0
  fi

  # Timestamp parsing from basename: YYYYMMDDHH_* -> extract YYYY, MM, DD
  local ts="${bn%%_*}"
  local suffix="${bn#${ts}}"
  local yyyy="${ts:0:4}" mm="${ts:4:2}" dd="${ts:6:2}"

  # Derive RHO file path (2019 borrowing for 2020 + leap-day handling)
  local rho_file
  if [[ "$yyyy$mm$dd" == "20200229" && $IS_SPECIAL -eq 1 ]]; then
    local prev="${rhopath}/2019022800${suffix}"
    local next="${rhopath}/2019030100${suffix}"
    if [[ ! -s "$prev" || ! -s "$next" ]]; then
      record_failure "rho_leap_missing" "$bn" "prev=$(basename "$prev") next=$(basename "$next")"
      echo "RHO LEAP MISSING for $bn" >&2
      return 0
    fi
    rho_file="${TMPDIR:-${SLURM_TMPDIR:-/tmp}}/RHO_2020022900${suffix}"
    if ! cdo -L -O ensmean "$prev" "$next" "$rho_file"; then
      record_failure "rho_leap_ensmean_fail" "$bn" "cdo ensmean"
      echo "RHO ENSMEAN FAIL: $bn" >&2
      return 0
    fi
  else
    # Borrow RHO from 2019 for 2020 and 2021 (same month/day/hour/suffix)
    local rho_year="$yyyy"
    if [[ "$yyyy" == "2020" || "$yyyy" == "2021" ]]; then
      rho_year="2019"
    fi
    rho_file="${rhopath}/${rho_year}${mm}${dd}00${suffix}"
  fi

  if [[ ! -s "$rho_file" ]]; then
    record_failure "rho_missing" "$bn" "$(basename "$rho_file")"
    echo "RHO MISSING: $bn (expected $(basename "$rho_file"))" >&2
    echo "Rhopath components rhopath $rhopath rho_year $rho_year mm $mm dd $dd suffix $suffix"
    return 0
  fi

  # Identify RHO variable name inside rho_file (using CDO)
  local rho_var="$RHO_VARNAME"
  if ! cdo -s showname "$rho_file" | tr '[:space:]' '\n' | awk 'NF>0' | grep -Fxq "$rho_var"; then
    for try in RHO rho RHO_SFC air_density; do
      if cdo -s showname "$rho_file" | tr '[:space:]' '\n' | awk 'NF>0' | grep -Fxq "$try"; then
        rho_var="$try"
        break
      fi
    done
  fi
  if ! cdo -s showname "$rho_file" | tr '[:space:]' '\n' | awk 'NF>0' | grep -Fxq "$rho_var"; then
    record_failure "rho_var_not_found" "$bn" "$(basename "$rho_file")"
    echo "RHO VAR NOT FOUND in $(basename "$rho_file") for $bn" >&2
    return 0
  fi

  # Merge scenario vars from this file + RHO var from rho_file
  local scen_csv; scen_csv=$(echo "$present" | paste -sd,)
  local work_nc="${TMPDIR:-${SLURM_TMPDIR:-/tmp}}/work_${bn}"
  local expr_nc="${TMPDIR:-${SLURM_TMPDIR:-/tmp}}/expr_${bn}"

  if ! cdo -L -O merge -selname,"$scen_csv" "$f" -selname,"$rho_var" "$rho_file" "$work_nc"; then
    record_failure "cdo_merge_fail" "$bn" "cdo merge"
    echo "MERGE FAIL: $bn" >&2
    return 0
  fi
  if ! cdo -L -O expr,"$EXPR_FORMULA" "$work_nc" "$expr_nc"; then
    record_failure "cdo_expr_fail" "$bn" "cdo expr"
    echo "EXPR FAIL: $bn" >&2
    return 0
  fi
  # Clip tiny negative BASEPRIM25 values
  if ! ncap2 -O -s 'where(BASEPRIM25<0) BASEPRIM25=0' "$expr_nc" "$expr_nc"; then
    record_failure "clip_prim25_fail" "$bn" "ncap2 clip"
    echo "CLIP FAIL: $bn" >&2
    return 0
  fi
  if ! cdo -L -O selname,"$VARS_KEEP" "$expr_nc" "$out_file"; then
    record_failure "cdo_selname_fail" "$bn" "$VARS_KEEP"
    echo "SELNAME FAIL: $bn" >&2
    return 0
  fi

  record_success "$bn"
  echo "SCEN OK: $bn -> ${out_file##*/} (RHO=$(basename "$rho_file") var=$rho_var)"
}

export -f merge_and_compute_base record_failure record_success dbg
export rhopath out_dir IS_BASEGM IS_SPECIAL RHO_VARNAME VARS_KEEP EXPR_FORMULA SURFACE_LEVEL
export status_dir fail_dir succ_dir summary_file DEBUG DEBUG_DIR TMPDIR OMP_NUM_THREADS CDO_PCTL_NTHREADS HDF5_DISABLE_VERSION_CHECK HDF5_USE_FILE_LOCKING

# =========================
# Build job list (REDO vs FULL)
# =========================
declare -a JOB_FILES=()

if [[ "$REDO_MODE" == "1" ]]; then
  echo "REDO_MODE=1: Using TSV ${TSV_LATEST}"
  if [[ ! -s "$TSV_LATEST" ]]; then
    echo "ERROR: REDO_MODE=1 but TSV not found: $TSV_LATEST" >&2
    exit 1
  fi
  while IFS=$'\t' read -r bn status _; do
    [[ -z "$bn" || "$status" == "status" ]] && continue    # skip header/empties
    [[ "$status" == "ok" ]] && continue
    bn="${bn##*/}"
    [[ -n "$bn" ]] && JOB_FILES+=( "${in_dir%/}/${bn}" )
  done < "$TSV_LATEST"

  if (( ${#JOB_FILES[@]} == 0 )); then
    echo "TSV shows 0 failed files — nothing to redo. Exiting."
    exit 0
  fi
  echo "Redoing ${#JOB_FILES[@]} files:"
  printf '  %s\n' "${JOB_FILES[@]}"
else
  echo "REDO_MODE=0: Full run over ${in_dir}"
  # Avoid sort -z; guard against empties
  while IFS= read -r -d '' f; do
    [[ -n "$f" ]] && JOB_FILES+=( "$f" )
  done < <(find "$in_dir" -maxdepth 1 -type f -name '*.nc' -print0)

  if (( ${#JOB_FILES[@]} == 0 )); then
    echo "ERROR: No .nc files in $in_dir" >&2
    exit 1
  fi
fi

# Filter any accidental empties (defensive)
if (( ${#JOB_FILES[@]} > 0 )); then
  cleaned=()
  for f in "${JOB_FILES[@]}"; do
    [[ -n "$f" ]] && cleaned+=( "$f" )
  done
  JOB_FILES=("${cleaned[@]}")
fi
dbg "Total JOB_FILES after filtering: ${#JOB_FILES[@]}"

# Persist job list (only when DEBUG=1)
if (( DEBUG )); then
  JOBLIST_FILE="${DEBUG_DIR}/job_files_${ARRAY_ID}.txt"
  {
    echo "# Job files for array ${ARRAY_ID} (in_dir=${in_dir})"
    echo "# Generated: $(date -Is)"
    i=0
    for f in "${JOB_FILES[@]}"; do
      printf '%06d\t%s\n' "$i" "$f"
      ((++i))
    done
  } > "$JOBLIST_FILE"
  dbg "Wrote job list to: $JOBLIST_FILE"
fi

# =========================
# Decide sequential vs parallel
# =========================
# Compute default workers from SLURM_CPUS_PER_TASK, allow override via MAX_JOBS
computed_workers=$(( ${SLURM_CPUS_PER_TASK:-16} / ${OMP_NUM_THREADS:-1} ))
(( computed_workers <= 0 )) && computed_workers=1
JOBS="${MAX_JOBS:-$computed_workers}"

dbg "Computed workers: $computed_workers; Using JOBS=$JOBS"

if (( JOBS <= 1 )); then
  # Sequential (explicit)
  for f in "${JOB_FILES[@]}"; do
    merge_and_compute_base "$f"
  done
else
  # Parallel
  if command -v parallel >/dev/null 2>&1; then
    # Export function definitions to a string and inject into each worker shell
    FUNC_DEF="$(
      declare -f merge_and_compute_base
      declare -f record_failure
      declare -f record_success
      declare -f dbg
    )"
    export FUNC_DEF

    # Ensure env vars are visible to workers
    export rhopath out_dir IS_BASEGM IS_SPECIAL RHO_VARNAME VARS_KEEP EXPR_FORMULA SURFACE_LEVEL
    export status_dir fail_dir succ_dir summary_file DEBUG DEBUG_DIR
    export TMPDIR OMP_NUM_THREADS CDO_PCTL_NTHREADS HDF5_DISABLE_VERSION_CHECK HDF5_USE_FILE_LOCKING

    printf '%s\0' "${JOB_FILES[@]}" \
    | parallel -0 -j "$JOBS" \
        --joblog "${status_dir}/parallel.joblog" \
        --env FUNC_DEF \
        --env rhopath --env out_dir --env IS_BASEGM --env IS_SPECIAL \
        --env RHO_VARNAME --env VARS_KEEP --env EXPR_FORMULA --env SURFACE_LEVEL \
        --env status_dir --env fail_dir --env succ_dir --env summary_file \
        --env DEBUG --env DEBUG_DIR \
        --env TMPDIR --env OMP_NUM_THREADS --env CDO_PCTL_NTHREADS \
        --env HDF5_DISABLE_VERSION_CHECK --env HDF5_USE_FILE_LOCKING \
        'bash -lc "$FUNC_DEF; merge_and_compute_base \"{}\""' 
  else
    # fallback to sequential if GNU parallel not available
    for f in "${JOB_FILES[@]}"; do merge_and_compute_base "$f"; done
  fi
fi

# =========================
# Summary
# =========================
total_inputs=${#JOB_FILES[@]}
succ_count=$(find "$succ_dir" -type f -name '*.ok' 2>/dev/null | wc -l | awk '{print $1}')
fail_count=$(find "$fail_dir" -type f -name '*.txt' 2>/dev/null | wc -l | awk '{print $1}')
{
  echo "===== SUMMARY for prefix ${prefix} (array ${ARRAY_ID}) ====="
  echo "Processed inputs: $total_inputs"
  echo "Succeeded: $succ_count"
  echo "Failed/Skipped: $fail_count"
  if (( fail_count > 0 )); then
    echo
    echo "---- Failed files (category :: file :: reason) ----"
    while IFS= read -r f; do
      bn=$(awk -F= '/^FILE=/{print $2}' "$f")
      catg=$(awk -F= '/^CATEGORY=/{print $2}' "$f")
      rsn=$(awk -F= '/^REASON=/{print $2}' "$f")
      echo "${catg} :: ${bn} :: ${rsn}"
    done < <(find "$fail_dir" -type f -name '*.txt' | sort)
  fi
} | tee "$summary_file" >&2
echo "Summary: $summary_file"
