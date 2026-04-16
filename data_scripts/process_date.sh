#!/usr/bin/env bash
# process_date.sh
# Usage: process_date.sh YEAR MM DD HH [--in-gemmach DIR] [--in-rdps DIR] [--out-gemmach DIR] [--out-rdps DIR] [--gridfile FILE] [--no-rdps-remap]
# Example:
#   IN_GEMMACH_ROOT=/project/... IN_RDPS_ROOT=/project/... OUT_GEMMACH_ROOT=/scratch/... OUT_RDPS_ROOT=/scratch/... \
#   GRIDFILE=/project/.../gridfile_lcc.txt ./process_date.sh 2019 01 01 00

set -euo pipefail

# -----------------------------
# Positional args
# -----------------------------
if [[ $# -lt 4 ]]; then
  echo "Usage: $0 YEAR MM DD HH [--in-gemmach DIR] [--in-rdps DIR] [--out-gemmach DIR] [--out-rdps DIR] [--gridfile FILE] [--no-rdps-remap]"
  exit 0
fi
year="$1"; month="$2"; day="$3"; hour="$4"

# Validate date (skip if bad -> "ignore stopped input")
if ! date -d "${year}-${month}-${day}" >/dev/null 2>&1; then
  echo "Invalid date: ${year}-${month}-${day} — skipping."
  exit 1
fi

# -----------------------------
# Same suffixes/mid as originals
# -----------------------------
suf="_00_00_gemmach.nc"
suf2="_00_00_rdpsqc.nc"
mid="00_0"

# -----------------------------
# Defaults (override by env or flags)
# -----------------------------
IN_GEMMACH_ROOT="${IN_GEMMACH_ROOT:-/home/tfmrodge/scratch/GEMMACH_data/data/BASEGM_2015_017}"
IN_RDPS_ROOT="${IN_RDPS_ROOT:-/home/tfmrodge/scratch/GEMMACH_data/data/RDPS_QC}"
OUT_GEMMACH_ROOT="${OUT_GEMMACH_ROOT:-/home/tfmrodge/scratch/GEMMACH_data/data/nc3/BASEGM_2015_017}"
OUT_RDPS_ROOT="${OUT_RDPS_ROOT:-/home/tfmrodge/scratch/GEMMACH_data/data/nc3/RDPS_QC}"
GRIDFILE="${GRIDFILE:-gridfile_lcc.txt}"
DO_RDPS_REMAP="${DO_RDPS_REMAP:-1}"

# -----------------------------
# Optional CLI overrides
# -----------------------------
shift 4 || true
while [[ $# -gt 0 ]]; do
  case "$1" in
    --in-gemmach)    IN_GEMMACH_ROOT="$2"; shift 2;;
    --in-rdps)       IN_RDPS_ROOT="$2"; shift 2;;
    --out-gemmach)   OUT_GEMMACH_ROOT="$2"; shift 2;;
    --out-rdps)      OUT_RDPS_ROOT="$2"; shift 2;;
    --gridfile)      GRIDFILE="$2"; shift 2;;
    --no-rdps-remap) DO_RDPS_REMAP=0; shift 1;;
    *) echo "Unknown option: $1" >&2; exit 2;;
  esac
done

# -----------------------------
# Node-local temp workspace
# Prefer $SLURM_TMPDIR when available (Alliance best practice)
# -----------------------------
TMP_BASE="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"
WORKDIR="${TMP_BASE}/gemmach_processing_${SLURM_JOB_ID:-$$}_${year}${month}${day}${hour}"
mkdir -p "${WORKDIR}"
cleanup() { rm -rf "${WORKDIR}" 2>/dev/null || true; }
trap cleanup EXIT INT TERM

# -----------------------------
# Compose paths, mirroring your originals
# -----------------------------
# Final outputs (used also as existence checks)
outgem="${OUT_GEMMACH_ROOT}/${year}-${month}-${day}_${hour}${suf}"
file2="${OUT_RDPS_ROOT}/${year}-${month}-${day}_${hour}${suf2}"

# Raw inputs
ingem="${IN_GEMMACH_ROOT}/${year}${month}${day}${mid}${hour}.nc"
inrdps="${IN_RDPS_ROOT}/${year}${month}${day}${mid}${hour}.nc"

# Intermediates in node-local storage
tmp_gem1="${WORKDIR}/${year}-${month}-${day}_${hour}${suf}-tmp1.nc"
tmp_gem2="${WORKDIR}/${year}-${month}-${day}_${hour}${suf}-tmp2.nc"
tmp_rdps="${WORKDIR}/${year}-${month}-${day}_${hour}${suf2}-tmp.nc"

# Ensure output dirs exist
mkdir -p "${OUT_GEMMACH_ROOT}" "${OUT_RDPS_ROOT}"

# -----------------------------
# GEMMACH pipeline (close to your originals)
# 1) ncap2 compute pressure
# 2) cdo remapbil to GRIDFILE
# 3) ncatted _FillValue, ncpdq pack+compress to final
# -----------------------------
if [[ ! -f "${outgem}" ]]; then
  if [[ -f "${ingem}" ]]; then
    if [[ -f "${GRIDFILE}" ]]; then
      echo "[GEMMACH] cdo remapbil -> ${tmp_gem1}"
      cdo -b F32 remapbil,"${GRIDFILE}" "${ingem}" "${tmp_gem1}"
      #Bring back a_1 and b_1
      ncks -A -v a_1,b_1 "${ingem}" "${tmp_gem1}"
    else
      echo "[GEMMACH] WARNING: GRIDFILE not found (${GRIDFILE}); skipping."
      exit 1
    fi
    echo "[GEMMACH] ncap2 pressure -> ${outgem}"
    ncap2 --no_tmp_fl -6 -O -s 'pressure=float(exp(a_1+b_1*ln(level1))/100)' "${tmp_gem1}" "${outgem}"
    # echo "[GEMMACH] ncks -> ${outgem}"
    # #ncatted -O -a _FillValue,,o,f,-32767 "${tmp_gem2}"
    # #ncpdq -P pack --ppc default=6 -6 -O "${tmp_gem2}" "${file1}"
    # ncks --no_tmp_fl -6 -O "${tmp_gem2}" "${file1}"
  else
    echo "[GEMMACH] Missing input: ${ingem} — skipping."
  fi
else
  echo "[GEMMACH] Final exists: ${file1} — skipping."
fi

# -----------------------------
# RDPS pipeline
# If final missing and input present:
# 1) optionally cdo remapbil to GRIDFILE
# 2) ncks -6 to final
# -----------------------------
if [[ ! -f "${file2}" ]]; then
  if [[ -f "${inrdps}" ]]; then
    if [[ "${DO_RDPS_REMAP}" -eq 1 && -f "${GRIDFILE}" ]]; then
      echo "[RDPS] cdo remapbil -> ${tmp_rdps}"
      cdo -b F32 remapbil,"${GRIDFILE}" "${inrdps}" "${tmp_rdps}"
      echo "[RDPS] ncks -6 -> ${file2}"
      ncks -6 -O "${tmp_rdps}" "${file2}"
    else
      echo "[RDPS] no remap -> ${tmp_rdps}"
      exit 1
    fi
  else
    echo "[RDPS] Missing input: ${inrdps} — skipping."
  fi
else
  echo "[RDPS] Final exists: ${file2} — skipping."
fi

echo "[DONE] ${year}-${month}-${day} ${hour} — intermediates were in ${WORKDIR}"