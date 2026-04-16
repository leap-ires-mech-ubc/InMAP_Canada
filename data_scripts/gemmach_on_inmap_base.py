#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Windowed zonal statistics for GEMMACH NetCDF variables over InMAP shapefile geometries.
- Robust to missing NetCDF georeferencing by materializing per-variable GeoTIFFs (CRS + transform).
- Parallelized with joblib over geometry chunks, using rasterio windowed reads.

Author: Timothy Rodgers + M365 Copilot
Last updated: 2026-01-28
"""

# ===================== Early environment setup (HPC + PROJ/GDAL) =====================
import os
import os.path as op

# Avoid BLAS/OMP oversubscription inside joblib workers
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

# Ensure PROJ & GDAL data dirs are configured BEFORE heavy GIS imports (geopandas/rasterio)
try:
    import pyproj
    from pyproj import datadir

    proj_env = os.environ.get("PROJ_LIB")
    if not proj_env or not (op.isdir(proj_env) and op.exists(op.join(proj_env, "proj.db"))):
        # Prefer module path if available
        mod_root = os.environ.get("EBROOTPROJ")
        if mod_root and op.exists(op.join(mod_root, "share", "proj", "proj.db")):
            proj_env = op.join(mod_root, "share", "proj")
            os.environ["PROJ_LIB"] = proj_env
            os.environ["PROJ_DATA"] = proj_env
        else:
            # Fallback: bundled pyproj data (varies by wheel)
            d = op.dirname(pyproj.__file__)
            for c in (
                op.join(d, "proj_dir", "share", "proj"),
                op.join(d, "proj_data"),
                op.join(d, "share", "proj"),
            ):
                if op.exists(op.join(c, "proj.db")):
                    proj_env = c
                    os.environ["PROJ_LIB"] = c
                    os.environ["PROJ_DATA"] = c
                    break
    # Tell pyproj explicitly
    if proj_env and op.exists(op.join(proj_env, "proj.db")):
        datadir.set_data_dir(proj_env)
except Exception:
    pass

# ===================== Standard imports =====================
import sys
import time
import math
import traceback
import multiprocessing as mp
from datetime import datetime
import socket

import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import rioxarray  # needed for .rio accessor on DataArray
from joblib import Parallel, delayed

import rasterio
import rasterio.windows
import rasterio.features
from affine import Affine  # for building transform and identity comparisons


# After rasterio import, help GDAL find its data dir (optional)
try:
    gdal_data = op.join(op.dirname(rasterio.__file__), "gdal_data")
    if op.exists(gdal_data):
        os.environ.setdefault("GDAL_DATA", gdal_data)
except Exception:
    pass


# ===================== Debugger (optional attach) =====================
def start_debugpy_if_requested(default_port: int = 5678):
    """
    Enable remote debugging if ENABLE_DEBUGPY=1 is set in env.
    Use DEBUGPY_PORT to override port. Exits cleanly with a helpful message
    if the port is already in use.
    """
    if os.environ.get("ENABLE_DEBUGPY") != "1":
        return
    port = int(os.environ.get("DEBUGPY_PORT", default_port))

    # Quick check: is port already listening?
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.settimeout(0.2)
        if s.connect_ex(("127.0.0.1", port)) == 0:
            print(
                f"⚠️  Port {port} is already in use. "
                f"Kill the listener or set DEBUGPY_PORT to a free port and re-run."
            )
            sys.exit(2)

    import debugpy

    try:
        debugpy.listen(("0.0.0.0", port))
    except OSError as ex:
        print(f"❌ debugpy.listen failed on port {port}: {ex}")
        print("   Tip: export DEBUGPY_PORT=5680 (or another) and try again.")
        sys.exit(2)

    print(f"🔧 debugpy: listening on {port} (waiting for VS Code to attach...)")
    debugpy.wait_for_client()
    print("🔧 debugpy: client attached, continuing")


# Call as early as possible (you can run Python with -Xfrozen_modules=off to avoid breakpoint issues)
start_debugpy_if_requested()


# ===================== Logging helpers =====================
def now() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def eprint(msg: str):
    sys.stderr.write(f"[{now()}] {msg}\n")
    sys.stderr.flush()


def iprint(msg: str):
    print(f"[{now()}] {msg}", flush=True)


# ===================== Scenario mapping & variable list =====================
SCENARIOS = {
    1: ["BASEGM_2015_017", "BASEGM_2015_017"],
    2: ["BAU_2020_E108", "BAU_2020_105"],
    3: ["BAU_2025_E107", "BAU_2025_106"],
    4: ["BAU_2030_E107", "BAU_2030_106"],
    5: ["BAU_2035_E107", "BAU_2035_106"],
    6: ["BPT_2015_E002", "BPT_2015_002"],
    7: ["BPT_2015_E006", "BPT_2015_006"],
    8: ["BPT_2015_E015", "BPT_2015_015"],
    9: ["COVID_2020_E003", "COVID_2020_003"],
    10: ["BPT_2015_E016", "BPT_2015_016"],
}

# NetCDF variable names to process (must match exactly)
RASTCOLS = ["BASEPM25", "BASEPNO3", "BASEPNH4", "BASEPSO4", "BASESOA", "BASEPRIM25"]


# ===================== Utility: chunking =====================
def chunker(seq, size):
    """Yield successive chunks from a sequence."""
    for pos in range(0, len(seq), size):
        yield seq[pos : pos + size]


# ===================== RasterIO helpers =====================
def list_subdatasets(nc_path: str):
    """List GDAL subdatasets (variables) inside a NetCDF."""
    with rasterio.open(nc_path) as ds:
        return list(ds.subdatasets)


def resolve_subdataset_path(nc_path: str, varname: str) -> str:
    """
    Return a GDAL subdataset path matching the requested variable (case-insensitive).
    If the dataset has no subdatasets (e.g., a GeoTIFF), return nc_path (assume single-band).
    """
    subs = list_subdatasets(nc_path)
    if not subs:
        iprint(f"No subdatasets listed in '{nc_path}'. Assuming single-band raster.")
        return nc_path

    # Try exact & case-insensitive match: NETCDF:"/file.nc":varname
    candidates = [s for s in subs if s.endswith(f":{varname}")]
    if not candidates:
        var_lc = varname.lower()
        candidates = [s for s in subs if s.split(":")[-1].lower() == var_lc]

    if not candidates:
        # Fallback: partial match
        candidates = [s for s in subs if varname.lower() in s.lower()]

    if not candidates:
        raise ValueError(
            f"Variable '{varname}' not found in NetCDF subdatasets:\n"
            + "\n".join(subs[:20] + (["..."] if len(subs) > 20 else []))
        )

    if len(candidates) > 1:
        iprint(f"Multiple matches for '{varname}', using first:\n    {candidates[0]}")

    return candidates[0]


def zonal_stats_windowed_opened(src: rasterio.io.DatasetReader, geom, stats=('sum',)):
    """
    Compute zonal stats using a pre-opened rasterio dataset and windowed reads.
    Reads only the window overlapping the geometry bounding box. Robust to edge cases.
    """
    # 1) Intersect the polygon bounds with the raster bounds
    gxmin, gymin, gxmax, gymax = geom.bounds
    rb = src.bounds  # (left, bottom, right, top)

    left   = max(gxmin, rb.left)
    bottom = max(gymin, rb.bottom)
    right  = min(gxmax, rb.right)
    top    = min(gymax, rb.top)

    # 2) If no overlap, short-circuit
    if (right <= left) or (top <= bottom):
        return {s: (0.0 if s == 'sum' else np.nan) for s in stats}

    # 3) Compute the window for the intersected bounds
    try:
        window = rasterio.windows.from_bounds(left, bottom, right, top, transform=src.transform)
    except Exception as ex:
        # Extra guard: if anything odd slips through (rare), return safe defaults
        # Helps avoid "Bounds and transform are inconsistent"
        # You can also iprint a one-liner if you want to log it.
        return {s: (0.0 if s == 'sum' else np.nan) for s in stats}

    if window.width <= 0 or window.height <= 0:
        return {s: (0.0 if s == 'sum' else np.nan) for s in stats}

    # 4) Read only that window
    data = src.read(1, window=window, masked=True)
    if data.size == 0:
        return {s: (0.0 if s == 'sum' else np.nan) for s in stats}

    # 5) Rasterize a polygon mask inside this window and compute stats
    try:
        mask = rasterio.features.geometry_mask(
            [geom],
            out_shape=(data.shape[0], data.shape[1]),
            transform=src.window_transform(window),
            invert=True
        )
    except Exception:
        return {s: (0.0 if s == 'sum' else np.nan) for s in stats}

    combined = np.ma.array(data, mask=(data.mask | ~mask))
    vals = combined.compressed()

    if vals.size == 0:
        return {s: (0.0 if s == 'sum' else np.nan) for s in stats}

    out = {}
    for s in stats:
        if s == 'sum':
            out[s] = float(np.nansum(vals))
        elif s == 'mean':
            out[s] = float(np.nanmean(vals))
        elif s == 'min':
            out[s] = float(np.nanmin(vals))
        elif s == 'max':
            out[s] = float(np.nanmax(vals))
        elif s == 'count':
            out[s] = int(vals.size)
        else:
            raise ValueError(f"Unsupported stat '{s}'")
    return out


def process_chunk(nc_path: str, subds_path: str, geoms_chunk, stats=("sum",)):
    """Open the subdataset once per chunk and process all geometries in that chunk."""
    results = []
    with rasterio.open(subds_path if subds_path else nc_path) as src:
        for g in geoms_chunk:
            res = zonal_stats_windowed_opened(src, g, stats=stats)
            results.append(res)
    return results


def zonal_stats_parallel_windowed(
    shape_gdf: gpd.GeoDataFrame,
    nc_path: str,
    varname: str,
    stats=("sum",),
    n_jobs: int | None = None,
    chunk_size: int = 300,
) -> pd.DataFrame:
    """
    Parallel zonal stats for one variable in a NetCDF (via rasterio subdataset) or GeoTIFF.
    Returns a DataFrame of columns: {varname}_{stat}
    """
    if n_jobs in (None, 0, -1):
        # Respect SLURM allocation if present
        n_jobs = int(os.environ.get("SLURM_CPUS_PER_TASK", mp.cpu_count()))

    # Resolve the subdataset path for the variable (no-op for TIFFs)
    subds_path = resolve_subdataset_path(nc_path, varname)
    iprint(f"Variable '{varname}' → subdataset:\n  {subds_path}")

    geoms = list(shape_gdf.geometry)
    N = len(geoms)
    iprint(f"Computing stats for {N} polygons | var={varname} | n_jobs={n_jobs} | chunk={chunk_size}")

    results_all = []
    t0 = time.time()

    # Submit chunked jobs
    for chunk_id, geoms_chunk in enumerate(chunker(geoms, chunk_size), start=1):
        t_chunk0 = time.time()
        chunk_results = Parallel(n_jobs=n_jobs, backend="loky")(
            delayed(process_chunk)(nc_path, subds_path, [g], stats=stats)
            for g in geoms_chunk
        )
        results_all.extend(r[0] for r in chunk_results)

        dt_chunk = time.time() - t_chunk0
        iprint(f"  chunk {chunk_id:04d}: {len(geoms_chunk)} geoms in {dt_chunk:,.2f}s")

    dt = time.time() - t0
    iprint(f"Finished var={varname} in {dt:,.2f}s")

    df = pd.DataFrame(results_all)
    df.columns = [f"{varname}_{s}" for s in stats]
    return df


# ===================== NetCDF → GeoTIFF materialization =====================
def _pick_xy_dims(da: xr.DataArray) -> tuple[str | None, str | None]:
    # Prefer explicit 'x','y'; else guess
    dims = list(da.dims)
    x_dim = next((d for d in dims if d.lower() == "x"), None)
    y_dim = next((d for d in dims if d.lower() == "y"), None)
    if not x_dim or not y_dim:
        # Heuristic: assume last 2 dims are spatial
        if len(dims) >= 2:
            y_dim, x_dim = dims[-2], dims[-1]
    return x_dim, y_dim


def _squeeze_to_2d(da: xr.DataArray) -> xr.DataArray:
    # Drop/slice singleton dims like time=0, level=0 to get (y,x)
    for d in list(da.dims):
        if d not in ("y", "x") and da.sizes.get(d, 1) == 1:
            da = da.isel({d: 0}, drop=True)
    return da


def materialize_geotiff(nc_path: str, varname: str, out_dir: str, overwrite: bool = False) -> str:
    """
    Create a GeoTIFF for `varname` in `nc_path` with a valid CRS + affine transform.
    Returns the GeoTIFF path. Reuses existing file unless overwrite=True.
    """
    os.makedirs(out_dir, exist_ok=True)
    stem = op.basename(nc_path).replace(".nc", "")
    out_tif = op.join(out_dir, f"{stem}_{varname}.tif")
    if (not overwrite) and op.exists(out_tif):
        return out_tif

    ds = xr.open_dataset(nc_path)  # netCDF engine
    if varname not in ds:
        raise ValueError(f"Variable '{varname}' not found in {nc_path}")

    da = _squeeze_to_2d(ds[varname])
    x_dim, y_dim = _pick_xy_dims(da)
    if x_dim is None or y_dim is None:
        raise ValueError(f"Could not identify spatial dims for {varname} in {nc_path}. Dims={da.dims}")

    # Ensure coordinate vectors exist
    if x_dim not in da.coords or y_dim not in da.coords:
        raise ValueError(f"Missing {x_dim}/{y_dim} coordinates for {varname} in {nc_path}")

    # Tell rioxarray the spatial dims
    da = da.rio.set_spatial_dims(x_dim=x_dim, y_dim=y_dim)

    # Build transform from x/y (assumes regular grid; north-up)
    x = da.coords[x_dim].values
    y = da.coords[y_dim].values
    if x.size < 2 or y.size < 2:
        raise ValueError(f"Not enough coordinate values to derive transform in {varname}")

    dx = float(x[1] - x[0])
    dy = float(y[1] - y[0])
    transform = Affine(dx, 0, x.min() - dx / 2, 0, -abs(dy), y.max() + abs(dy) / 2)
    da = da.rio.write_transform(transform)

    # Write CRS from grid_mapping (prefer WKT)
    gm_name = da.attrs.get("grid_mapping")
    if not gm_name:
        gm_name = "crs" if "crs" in ds.variables else ("lambert_conformal_conic" if "lambert_conformal_conic" in ds.variables else None)

    crs_text = None
    if gm_name and gm_name in ds.variables:
        gm = ds[gm_name]
        crs_text = gm.attrs.get("spatial_ref") or gm.attrs.get("crs_wkt")

    if crs_text is None:
        # Fallback: your known LCC (matches InMAP LCC used in shapefile)
        crs_text = (
            'PROJCS["Lambert_Conformal_Conic",GEOGCS["GCS_unnamed ellipse",'
            'DATUM["unknown",SPHEROID["Unknown",6378137,298.257222101]],'
            'PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],'
            'PROJECTION["Lambert_Conformal_Conic_2SP"],'
            'PARAMETER["latitude_of_origin",40],PARAMETER["central_meridian",-96],'
            'PARAMETER["standard_parallel_1",50],PARAMETER["standard_parallel_2",70],'
            'PARAMETER["false_easting",0],PARAMETER["false_northing",0],'
            'UNIT["metre",1]]'
        )
    da = da.rio.write_crs(crs_text)

    # Optional nodata
    if da.rio.nodata is None:
        da = da.rio.write_nodata(np.nan)

    # Write GeoTIFF atomically
    tmp = out_tif[:-4] + "_tmp.tif"
    da.rio.to_raster(tmp, compress="LZW", tiled=True, BIGTIFF="IF_SAFER")
    os.replace(tmp, out_tif)
    return out_tif


# ===================== Main =====================
def main(argv=None):
    if argv is None:
        argv = sys.argv

    if len(argv) < 2:
        eprint("Usage: python gemmach_on_inmap_base.py <scenario_id>")
        sys.exit(2)

    try:
        sID = int(argv[1])
    except Exception:
        eprint("First argument must be an integer scenario ID (1-10).")
        sys.exit(2)

    if sID not in SCENARIOS:
        eprint(f"Scenario ID {sID} not in valid range 1..{len(SCENARIOS)}")
        sys.exit(2)

    emissions, scenario = SCENARIOS[sID]
    iprint(f"Starting zonal on scenario_id={sID} | emissions={emissions} | scenario={scenario}")

    # ---- Paths ----
    ncpath = "/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/scenario_sums_LCC"
    nc_name = f"{ncpath}/20260120_{scenario}_annual_mean_LCC.nc"

    inmapoutpath = "/home/tfmrodge/scratch/GEMMACH_data/data/Inmap_outputs/Outputs"
    shp_path = f"{inmapoutpath}/20250930_InmapOuts_{emissions}_static.shp"
    out_path = f"{inmapoutpath}/20260128_InmapOuts_{emissions}_static_withGEMMACH.shp"

    # ---- Pre-materialize TIFFs once per variable (serial; avoids worker races) ----
    tif_cache_dir = f"/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/geotiff_cache/{scenario}"
    os.makedirs(tif_cache_dir, exist_ok=True)

    tif_map: dict[str, str | None] = {}
    for rc in RASTCOLS:
        try:
            iprint(f"Materializing GeoTIFF for {rc} (if needed)...")
            tif_map[rc] = materialize_geotiff(nc_name, rc, tif_cache_dir, overwrite=False)
        except Exception as ex:
            eprint(f"Failed to materialize {rc}: {ex}")
            tif_map[rc] = None  # we will try NetCDF subdataset fallback if needed

    # ---- Read shapefile ----
    iprint(f"Reading shapefile: {shp_path}")
    t0 = time.time()
    try:
        outshp = gpd.read_file(shp_path)
    except Exception as ex:
        eprint(f"Failed to read shapefile '{shp_path}': {ex}")
        traceback.print_exc(file=sys.stderr)
        sys.exit(3)
    iprint(f"Loaded {len(outshp)} polygons in {time.time()-t0:,.2f}s | CRS={outshp.crs}")

    # ---- Determine raster CRS using materialized first variable if available ----
    first_var = RASTCOLS[0]
    try:
        open_path_0 = tif_map.get(first_var) or resolve_subdataset_path(nc_name, first_var)
        with rasterio.open(open_path_0) as src0:
            raster_crs = src0.crs
            raster_transform = src0.transform
            raster_shape = (src0.height, src0.width)
            nodata = src0.nodata
        iprint(f"Raster opened for '{first_var}' | CRS={raster_crs} | shape={raster_shape} | nodata={nodata}")
    except Exception as ex:
        eprint(f"Failed opening raster for '{first_var}': {ex}")
        traceback.print_exc(file=sys.stderr)
        sys.exit(4)

    # Reproject shapes to raster CRS if necessary
    if raster_crs is None:
        eprint("WARNING: Raster CRS is None. Proceeding without reprojection (results may be incorrect).")
        working_gdf = outshp
    else:
        if outshp.crs != raster_crs:
            iprint(f"Reprojecting shapes from {outshp.crs} to raster CRS {raster_crs} for correct spatial overlay.")
            t1 = time.time()
            working_gdf = outshp.to_crs(raster_crs)
            iprint(f"Reprojected shapes in {time.time()-t1:,.2f}s")
        else:
            working_gdf = outshp

    # ---- Determine CPU usage ----
    n_jobs_env = os.environ.get("SLURM_CPUS_PER_TASK")
    n_jobs = int(n_jobs_env) if n_jobs_env else mp.cpu_count()
    iprint(f"Using n_jobs={n_jobs} (SLURM_CPUS_PER_TASK={n_jobs_env})")

    # ---- Compute zonal stats for each variable ----
    all_cols_df = outshp.reset_index(drop=True).copy()  # to attach results while preserving original CRS/geometry
    stats = ("sum",)
    chunk_size = 200  # tune (100-500 works well)

    t_all0 = time.time()
    for rc in RASTCOLS:
        try:
            iprint(f"==== Variable {rc} ====")
            input_path = tif_map.get(rc)
            if input_path and op.exists(input_path):
                subds_path = input_path  # use GeoTIFF directly
            else:
                subds_path = resolve_subdataset_path(nc_name, rc)  # fallback to NetCDF subdataset

            # Safety: if NetCDF subdataset has no georef, force TIFF path (overwrite)
            with rasterio.open(subds_path) as _src:
                if (_src.crs is None) or (_src.transform == Affine.identity):
                    iprint("CRS/transform missing; forcing GeoTIFF path")
                    subds_path = materialize_geotiff(nc_name, rc, tif_cache_dir, overwrite=True)

            # Run zonal stats
            df_rc = zonal_stats_parallel_windowed(
                shape_gdf=working_gdf,
                nc_path=subds_path,   # TIFF or subdataset path
                varname=rc,           # ignored if TIFF
                stats=stats,
                n_jobs=n_jobs,
                chunk_size=chunk_size,
            )
            all_cols_df = pd.concat([all_cols_df, df_rc], axis=1)

        except Exception as ex:
            eprint(f"Variable '{rc}' failed: {ex}")
            traceback.print_exc(file=sys.stderr)
            # Keep schema: fill with NaNs
            for s in stats:
                all_cols_df[f"{rc}_{s}"] = np.nan

    iprint(f"All variables processed in {time.time()-t_all0:,.2f}s")

    # ---- Write output shapefile with original CRS ----
    final_gdf = gpd.GeoDataFrame(all_cols_df, geometry=outshp.geometry, crs=outshp.crs)
    iprint(f"Writing output: {out_path}")
    tW = time.time()
    try:
        final_gdf.to_file(out_path)
    except Exception as ex:
        eprint(f"Failed writing shapefile '{out_path}': {ex}")
        traceback.print_exc(file=sys.stderr)
        sys.exit(5)
    iprint(f"Wrote {out_path} in {time.time()-tW:,.2f}s")

    iprint("Done.")


if __name__ == "__main__":
    try:
        main()
    except Exception as ex:
        eprint(f"FATAL: {ex}")
        traceback.print_exc(file=sys.stderr)
        sys.exit(99)