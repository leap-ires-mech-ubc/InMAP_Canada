#!/bin/bash
#SBATCH --time=0:59:59
#SBATCH --account=def-agiang01 #def-rscholes
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=350M # --mem=250G
#SBATCH --array=1
#SBATCH --job-name='check_vars_coverage'
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -euo pipefail
trap 'status=$?; echo "ERROR at line ${LINENO}: ${BASH_COMMAND}" >&2; exit $status' ERR
if [[ "${DEBUG:-0}" == "1" ]]; then set -x; fi

# ===== Modules =====
module load nco cdo
# If parallel is a module on your cluster, the next line ensures it is present:
module load parallel 2>/dev/null || true

# ===== Environment =====
# Concurrency defaults
JOBS="${SLURM_CPUS_PER_TASK:-8}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
export CDO_PCTL_NTHREADS="$OMP_NUM_THREADS"
export HDF5_DISABLE_VERSION_CHECK=1
export HDF5_USE_FILE_LOCKING=FALSE

# ===== Inputs (array-based prefix selection) =====
PREFIX_MAP=/home/tfmrodge/projects/def-rscholes/tfmrodge/InMAP_Canada/data_scripts/config_gemmachouts.txt
ARRAY_ID="${SLURM_ARRAY_TASK_ID:-${1:-}}" #ARRAY_ID=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is not set}
prefix=$(awk -v ArrayTaskID="$ARRAY_ID" '$1==ArrayTaskID {print $2}' "$PREFIX_MAP")
if [[ -z "${prefix}" ]]; then
  echo "ERROR: No prefix found for array id ${ARRAY_ID} in ${PREFIX_MAP}" >&2
  exit 1
fi

# ===== Paths (scan out_dir; report to report root) =====
#/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/BASEGM_surface
#/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/BAU_2020_105/base_vars
in_dir="/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/${prefix}"
out_dir="${in_dir}/20260129_base_vars"          # << directory to scan for files
#mkdir -p "$out_dir"

REPORT_ROOT="/home/tfmrodge/projects/def-rscholes/tfmrodge/InMAP_Canada/checklogs"
REPORT_DIR="${REPORT_ROOT%/}/${prefix}"
mkdir -p "$REPORT_DIR"

# ===== Check configuration =====
# Choose either YEAR or START/END (inclusive). Leave the other blank.
# if [[ "$prefix" == *2020* ]]; then
#   YEAR="2020"              # e.g., "2020" (empty string "" to disable)
# else
#   YEAR="2019"
# fi
# START=""                 # e.g., "2020-02-28"
# END=""                   # e.g., "2020-03-01"
# Special date range for 2020 prefixes:
if [[ "$prefix" == *2020* ]]; then
  START="2020-02-01"
  END="2021-01-31"
  #YEAR="2020"
else
  # Default: full calendar year (or your existing logic)
  START="2019-01-01"
  END="2019-12-31"
  #YEAR="2019"
fi

# Hour indexing: set to 1 for 000..023, set to 0 for 001..024
ZERO_BASED=1

# Variables to require (comma-separated exact names; no aliases/groups)
if [[ "$prefix" == *BASE* ]]; then
  VARS_CSV="RHO,BASEPM25,BASEPNO3,BASEPNH4,BASEPSO4,BASESOA,BASEPRIM25" #"RHO"           # e.g., "2020" (empty string "" to disable)
  out_dir="/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/BASEGM_2015_017/base_vars"
else
  VARS_CSV="BASEPM25,BASEPNO3,BASEPNH4,BASEPSO4,BASESOA,BASEPRIM25"
fi

IFS=',' read -r -a REQ_VARS <<< "$VARS_CSV"

# ===== Derived report filenames (UTC timestamp) =====
TS=$(date -u +%Y%m%dT%H%M%SZ)
TSV="${REPORT_DIR}/coverage_${prefix}_${TS}.tsv"
SUMMARY="${REPORT_DIR}/coverage_${prefix}_${TS}.summary.txt"
RUNLOG="${REPORT_DIR}/coverage_${prefix}_${TS}.run.log"

# ===== Tool checks =====
command -v cdo >/dev/null || { echo "ERROR: cdo not found on PATH"; exit 1; }
# Prefer GNU parallel; if missing, fallback to xargs -P
if command -v parallel >/dev/null 2>&1; then
  PARALLEL_TOOL="parallel --will-cite --jobs ${JOBS}"
else
  echo "WARN: GNU parallel not found; falling back to xargs -P ${JOBS}" | tee -a "$RUNLOG" >&2
  PARALLEL_TOOL=""
fi

# ===== Helper: hour generator =====
hours_for_day() {
  local start_h end_h
  if (( ZERO_BASED == 1 )); then start_h=0; end_h=23; else start_h=1; end_h=24; fi
  for (( h=start_h; h<=end_h; h++ )); do printf "%03d\n" "$h"; done
}

# ===== Checker (prints exactly one TSV row per file) =====
check_one() {
  local bn="$1"
  local f="${out_dir%/}/${bn}"

  if [[ ! -e "$f" ]]; then
    printf "%s\tmissing\tnot found\n" "$bn"
    return 0
  fi
  if [[ ! -s "$f" ]]; then
    printf "%s\tempty\tsize=0\n" "$bn"
    return 0
  fi
  if ! cdo -s showname "$f" >/dev/null 2>&1; then
    printf "%s\tunreadable\tcdo showname failed\n" "$bn"
    return 0
  fi

  local present
  present=$( cdo -s showname "$f" 2>/dev/null | grep -v '^cdi' | tr ' ' '\n' )

  local missing=()
  local v
  for v in "${REQ_VARS[@]}"; do
    grep -Fxq "$v" <<<"$present" || missing+=("$v")
  done

  if (( ${#missing[@]} > 0 )); then
    printf "%s\tvars_missing\t%s\n" "$bn" "$(IFS=,; echo "${missing[*]}")"
  else
    printf "%s\tok\t-\n" "$bn"
  fi
}

export -f check_one
export out_dir
export REQ_VARS
export ZERO_BASED

# ===== Build expected basenames =====
declare -a BNS=()

# if [[ -n "$YEAR" ]]; then
#   d="${YEAR}-01-01"
#   # Validate start date and existence of GNU date
#   date -d "$d" >/dev/null 2>&1 || { echo "ERROR: invalid YEAR: $YEAR"; exit 1; }
#   while : ; do
#     day=$(date -d "$d" +%Y%m%d)
#     while read -r hh; do
#       BNS+=("${day}00_${hh}.nc")
#     done < <(hours_for_day)
#     [[ "$d" == "${YEAR}-12-31" ]] && break
#     d=$(date -d "$d + 1 day" +%F)
#   done
#else
date -d "$START" >/dev/null 2>&1 || { echo "ERROR: invalid START"; exit 1; }
date -d "$END"   >/dev/null 2>&1 || { echo "ERROR: invalid END"; exit 1; }
d="$START"
while : ; do
  day=$(date -d "$d" +%Y%m%d)
  while read -r hh; do
    BNS+=("${day}00_${hh}.nc")
  done < <(hours_for_day)
  [[ "$d" == "$END" ]] && break
  d=$(date -d "$d + 1 day" +%F)
done
#fi

# ===== Run checks in parallel and write TSV =====
# Header
echo -e "basename\tstatus\tdetails" > "$TSV"

# Body: run checker in parallel; stderr to runlog
if command -v parallel >/dev/null 2>&1; then
    export -f check_one
    printf '%s\n' "${BNS[@]}" \
        | parallel --will-cite --jobs "$JOBS" check_one \
            >>"$TSV" 2>>"$RUNLOG"
else
    printf '%s\n' "${BNS[@]}" \
        | xargs -r -P "$JOBS" -I{} bash -c 'check_one "$@"' _ {} \
            >>"$TSV" 2>>"$RUNLOG"
fi

# ===== Summarize =====
# Counts
ok=$(awk -F'\t' 'NR>1 && $2=="ok"' "$TSV" | wc -l | awk '{print $1}')
miss=$(awk -F'\t' 'NR>1 && $2=="missing"' "$TSV" | wc -l | awk '{print $1}')
emp=$(awk -F'\t' 'NR>1 && $2=="empty"' "$TSV" | wc -l | awk '{print $1}')
bad=$(awk -F'\t' 'NR>1 && $2=="unreadable"' "$TSV" | wc -l | awk '{print $1}')
vbad=$(awk -F'\t' 'NR>1 && $2=="vars_missing"' "$TSV" | wc -l | awk '{print $1}')
tot=$(awk 'NR>1' "$TSV" | wc -l | awk '{print $1}')
notok=$((miss + emp + bad + vbad))

{
  echo "Prefix: ${prefix}"
  echo "Scanned directory: ${out_dir}"
  echo "Report TSV: ${TSV}"
  echo "Run log: ${RUNLOG}"
  echo
  echo "Counts"
  echo "------"
  echo "Total checked: $tot"
  echo "OK:            $ok"
  echo "Not OK total:  $notok"
  echo "  missing:     $miss"
  echo "  empty:       $emp"
  echo "  unreadable:  $bad"
  echo "  vars_missing:$vbad"
  echo
  echo "Files NOT OK (basename<TAB>status<TAB>details)"
  echo "----------------------------------------------"
  awk -F'\t' 'NR>1 && $2!="ok" {print $0}' "$TSV"
} | tee "$SUMMARY"

# (Optional) keep a "latest" symlink for convenience
ln -sfn "$(basename "$TSV")" "${REPORT_DIR}/coverage_${prefix}_latest.tsv" || true
ln -sfn "$(basename "$SUMMARY")" "${REPORT_DIR}/coverage_${prefix}_latest.summary.txt" || true