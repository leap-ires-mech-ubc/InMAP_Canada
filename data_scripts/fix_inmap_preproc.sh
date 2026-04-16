#!/usr/bin/env bash
# Fix NaN _FillValue in InMAP inputs and fill holes (nearest-neighbour recommended).
# Usage:
#   ./fix_inmap_preproc.sh INPUT.nc OUTPUT.nc [nn|fillmiss|zero]
#
#data_scripts/fix_inmap_preproc.sh /home/tfmrodge/scratch/GEMMACH_data/data/Inmap_outputs/Preproc/20260210/20260217_inmapData_GEMMACH_BASEGM_2015_017.nc /home/tfmrodge/scratch/GEMMACH_data/data/Inmap_outputs/Preproc/20260210/20260217_inmapData_GEMMACH_BASEGM_2015_017_cdoNN.nc
#
# Default method = nn (nearest-neighbour to same grid).
# Dependencies: cdo, nco (ncap2, ncdump)


set -euo pipefail

# -------------------------
# Parse args / defaults
# -------------------------
if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: $0 INPUT.nc OUTPUT.nc [nn|fillmiss|zero]"
  exit 1
fi

IN="$1"
OUT="$2"
METHOD="${3:-nn}"  # nn | fillmiss | zero

MISSVAL="1e20"     # numeric _FillValue to use
WORKDIR="$(mktemp -d)"
trap 'rm -rf "$WORKDIR"' EXIT

# -------------------------
# Dependency checks
# -------------------------
for cmd in cdo ncap2 ncdump; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "Error: '$cmd' not found in PATH. Please load/install it and retry." >&2
    exit 2
  fi
done

echo ">>> Input:    $IN"
echo ">>> Output:   $OUT"
echo ">>> Method:   $METHOD"
echo ">>> Workdir:  $WORKDIR"
echo ">>> MISSVAL:  $MISSVAL"
echo

TMP1="$WORKDIR/tmp_setmiss.nc"
TMP2="$WORKDIR/tmp_fixnan.nc"
GRIDTXT="$WORKDIR/grid.txt"
NCOSCRIPT="$WORKDIR/replace_nans.nco"

# -------------------------
# Step 1: Set a numeric _FillValue attribute globally
#         (metadata only; data values still contain NaN at this point)
# -------------------------
echo ">>> Setting _FillValue attribute to $MISSVAL for all variables ..."
cdo -L setmissval,"$MISSVAL" "$IN" "$TMP1"

# -------------------------
# Step 2: Replace actual NaN values with _FillValue for all float/double vars
#         Build one ncap2 script that applies to all variables in a single pass.
# -------------------------
echo ">>> Building NaN-replacement script for all float/double variables ..."
# Extract variable names from header (float/double only)
mapfile -t VARS < <(ncdump -h "$TMP1" \
  | awk '
      /variables:/ {vars=1; next}
      vars && NF==0 {exit}
      vars && ($1=="float" || $1=="double") {
        # $2 is like VarName(dim1,dim2,...)
        v=$2
        sub(/\(.*/,"",v)   # strip dimensions
        gsub(/;/,"",v)     # strip stray semicolons
        print v
      }
    ' \
  | sort -u)

if [[ ${#VARS[@]} -eq 0 ]]; then
  echo "Error: No float/double variables found to process." >&2
  exit 3
fi

# Generate an NCO script that replaces NaN with _FillValue(var) for each variable.
# Using 'where(isnan(var)) var=_FillValue(var);'
: > "$NCOSCRIPT"
for v in "${VARS[@]}"; do
  echo "where(isnan(${v})) ${v}=_FillValue(${v});" >> "$NCOSCRIPT"
done

echo ">>> Replacing NaN values with _FillValue in a single ncap2 pass ..."
ncap2 -O -S "$NCOSCRIPT" "$TMP1" "$TMP2"

# -------------------------
# (Optional) Quick sanity peek: count how many missings remain per var
# -------------------------
echo ">>> Quick check: variables and their _FillValue after NaN replacement:"
ncdump -h "$TMP2" | sed -n '/variables:/,/global attributes:/p' | grep '_FillValue' || true
echo

# -------------------------
# Step 3: Fill holes
#         Preferred: nearest-neighbour to the same grid (nn).
#         Fallback: fillmiss (if remapnn struggles with staggered fields).
#         Alternate: zero (not recommended, but provided on request).
# -------------------------
case "$METHOD" in
  nn)
    echo ">>> Filling holes with nearest-neighbour on the same grid ..."
    # Create a self-grid description and remap using nearest neighbour
    cdo griddes "$TMP2" > "$GRIDTXT" || true

    set +e
    cdo -L remapnn,"$GRIDTXT" "$TMP2" "$OUT"
    status=$?
    set -e
    if [[ $status -ne 0 ]]; then
      echo "!!! remapnn failed (possibly due to staggered dims). Falling back to fillmiss ..."
      cdo -L fillmiss "$TMP2" "$OUT"
    fi
    ;;

  fillmiss)
    echo ">>> Filling holes with CDO fillmiss (spatial interpolation) ..."
    cdo -L fillmiss "$TMP2" "$OUT"
    ;;

  zero)
    echo ">>> Replacing remaining missing values with zero (NOT recommended for met fields) ..."
    cdo -L setmisstoc,0 "$TMP2" "$OUT"
    ;;

  *)
    echo "Error: Unknown method '$METHOD' (use: nn | fillmiss | zero)" >&2
    exit 4
    ;;
esac

# -------------------------
# Step 4: Keep a consistent numeric _FillValue on the final file
# -------------------------
echo ">>> Enforcing consistent _FillValue=$MISSVAL on output ..."
cdo -L setmissval,"$MISSVAL" "$OUT" "$OUT.tmp" && mv "$OUT.tmp" "$OUT"

echo ">>> Done. Output written to: $OUT"