#!/bin/bash
#SBATCH --time=23:59:59
#SBATCH --account=def-rscholes
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=3900M
#SBATCH --array=1-9
#SBATCH --job-name='annual_avg_conc'
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -euo pipefail

# === Modules ===
module load nco cdo
# === Config / Inputs ===
PREFIX_MAP=/home/tfmrodge/projects/def-rscholes/tfmrodge/InMAP_Canada/data_scripts/config_gemmachouts.txt
ARRAY_ID=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is not set}
prefix=$(awk -v ArrayTaskID="$ARRAY_ID" '$1==ArrayTaskID {print $2}' "$PREFIX_MAP")
if [[ -z "${prefix}" ]]; then
  echo "ERROR: No prefix found for array id ${ARRAY_ID} in ${PREFIX_MAP}" >&2
  exit 1
fi
# Special-case handling
if [[ "$prefix" == "BASEGM_2015_017" ]]; then
    prefix="BASEGM_surface"
fi

echo "Array task ${ARRAY_ID} → prefix: ${prefix}"

# === Paths ===
in_dir="/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/${prefix}/"
out_dir="/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/scenario_sums3"
mkdir -p "$out_dir"

# === Testing limit (set MAXFILES="" to use all files) ===
MAXFILES="" #5   # set to "" to process all files 5
if [[ -n "${MAXFILES}" ]]; then
  echo "Limiting to first ${MAXFILES} files for testing."
  # Sorted, deterministic file list
  mapfile -t files < <(ls -1 "${in_dir}"*.nc | sort | head -n "${MAXFILES}")
  stamp=$(date +%Y%m%d)_test
else
  echo "Using ALL files in ${in_dir}"
  mapfile -t files < <(ls -1 "${in_dir}"*.nc | sort)
  stamp=$(date +%Y%m%d)
fi
# Print number of files found
echo "Number of files to process: ${#files[@]}"


if [[ ${#files[@]} -eq 0 ]]; then
  echo "ERROR: No .nc files found in ${in_dir}" >&2
  exit 1
fi

# === Output ===
outfile="${out_dir}/${stamp}_${prefix}_annual_mean.nc"

# === Scratch ===
tmpdir=$(mktemp -d)
trap 'rm -rf "${tmpdir}"' EXIT

echo "Input dir : ${in_dir}"
echo "Output file: ${outfile}"
echo "Temp dir   : ${tmpdir}"
echo "File count : ${#files[@]}"

# === Parallelism hint for CDO ===
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
CDO_THREADS=$OMP_NUM_THREADS
#CDO_THREADS=${SLURM_NTASKS:-1}

# === Step 1: Merge time (hourly across the whole period) ===
# Use netCDF4 + compression
echo "[1/4] Merging time across ${#files[@]} files..."
cdo -L -P "${CDO_THREADS}" -f nc4 -z zip_5 -mergetime "${files[@]}" "${tmpdir}/merged.nc"

# === Optional: select surface level if this is a surface-only run ===
# NOTE: Adjust the level value if required by your GEM-MACH coordinate.
if [[ "$prefix" == "BASEGM_surface" ]]; then
  echo "[2/4] Selecting surface level (level=0.99875)..."
  cdo -L -P "${CDO_THREADS}" -f nc4 -z zip_5 -sellevel,0.99875 \
      "${tmpdir}/merged.nc" "${tmpdir}/merged_sfc.nc"
  mv "${tmpdir}/merged_sfc.nc" "${tmpdir}/merged.nc"
fi

# === (Optional) If RHO is missing, derive it (example) ===
# Uncomment and adapt if needed:
# echo "Checking if RHO exists..."
# if ! ncks -m "${tmpdir}/merged.nc" | grep -qE '(^|\s)RHO(,|$)'; then
#   echo "RHO not found; deriving from pressure and temperature (adjust names/units!)."
#   # Example names/units (must match your files!):
#   # - Pressure in Pa: PS
#   # - Temperature in K: T
#   # RHO = PS / (R_d * T) with R_d = 287.05 J/(kg·K)
#   ncap2 -O -s 'RHO=PS/(287.05*T)' "${tmpdir}/merged.nc" "${tmpdir}/merged.nc"
# fi

# === Step 2: Create derived hourly variables BEFORE averaging ===
# This ensures we weight by contemporaneous density (RHO) and sum components correctly.
# Combine all expressions in ONE ncap2 call to avoid repeated IO.
echo "[3/4] Computing derived hourly variables with ncap2..."
ncap2 -O -s '
  BASEPM25 = AF;
  BASEPNO3 = (TNI1 + TNO3 + THN3) * RHO;
  BASEPNH4 = TAM1 * RHO;
  BASEPSO4 = TSU1 * RHO;
  BASESOA  = TOC1 * RHO;
  BASEPRIM25 = BASEPM25 - BASESOA - BASEPSO4 - BASEPNH4 - BASEPNO3 
' "${tmpdir}/merged.nc" "${tmpdir}/derived_hourly.nc"

# === Step 3: Annual means ===
# If multiple calendar years are present in the merged stream, yearmean creates one timestep per year.
# For a single year, you’ll get one timestep.
echo "[4/4] Computing annual means with CDO (yearmean)..."
cdo -L -P "${CDO_THREADS}" -f nc4 -z zip_5 yearmean "${tmpdir}/derived_hourly.nc" "${tmpdir}/annual_mean_full.nc"

# === (Optional) Keep only variables you care about ===
# Add existing variables you also want to retain (e.g., AF itself).
VARS="BASEPM25,BASEPNO3,BASEPNH4,BASEPSO4,BASESOA"
cdo -L -P "${CDO_THREADS}" -f nc4 -z zip_5 selvar,"${VARS}" "${tmpdir}/annual_mean_full.nc" "${outfile}"

# === Final touch: ensure metadata compression ===
# You can additionally chunk/compress via ncks if desired:
# ncks -O -4 -L 1 "${outfile}" "${outfile}"

#Add projection data

# --- Projection defaults (safe with `set -u`) ---
# If these env vars are already set, they will be respected; otherwise, defaults are assigned.
: "${proj4:='+proj=lcc +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +a=6378137 +rf=298.257222101 +units=m +no_defs'}"
: "${ulx:=-4184312.05377675}"
: "${uly:=4270177.170000}"
: "${lrx:=3165687.946000}"
: "${lry:=-2029822.82977676}"

# Target file to annotate (final output from your pipeline)
target="${outfile}"


# Ensure CF conventions tag
ncatted -O -a Conventions,global,o,c,"CF-1.8" "${target}"

# Create a grid-mapping variable 'crs' if missing
if ! ncks -m "${target}" 2>/dev/null | grep -qE "(^|[[:space:]])crs([[:space:]]|=|\()"; then
  ncap2 -O -s 'crs=0' "${target}" "${target}"
fi

# Populate CF attributes for Lambert Conformal Conic on 'crs'
# (values taken from your proj4 string)
ncatted -O \
  -a grid_mapping_name,crs,o,c,"lambert_conformal_conic" \
  -a longitude_of_central_meridian,crs,o,d,-96.0 \
  -a latitude_of_projection_origin,crs,o,d,40.0 \
  -a standard_parallel,crs,o,d,50.0,70.0 \
  -a false_easting,crs,o,d,0.0 \
  -a false_northing,crs,o,d,0.0 \
  -a semi_major_axis,crs,o,d,6378137.0 \
  -a inverse_flattening,crs,o,d,298.257222101 \
  -a units,crs,o,c,"m" \
  -a proj4_params,crs,o,c,"${proj4}" \
  "${target}"

# Global metadata: projected CRS & bounds (meters)
# Use x/y min/max because these are projected coordinates (not lon/lat)
ncatted -O \
  -a geospatial_bounds_crs,global,o,c,"${proj4}" \
  -a geospatial_x_min,global,o,d,${ulx} \
  -a geospatial_x_max,global,o,d,${lrx} \
  -a geospatial_y_min,global,o,d,${lry} \
  -a geospatial_y_max,global,o,d,${uly} \
  "${target}"


# Attach the CRS to the variables you keep in $VARS (and BASEPRIM25 if you keep it)
for v in BASEPM25 BASEPNO3 BASEPNH4 BASEPSO4 BASESOA BASEPRIM25; do
  if ncks -m "${target}" 2>/dev/null | grep -qE "(^|[[:space:]])${v}[[:space:]]*\\("; then
    ncatted -O -a grid_mapping,"${v}",o,c,"crs" "${target}"
  fi
done

echo "Done. Wrote: ${outfile}"