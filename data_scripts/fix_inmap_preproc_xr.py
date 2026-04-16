#!/usr/bin/env python3
"""
Clip a NetCDF dataset to the largest all-valid (no-NaN) bounding boxes on the horizontal grid,
and additionally write a second copy where NaNs are replaced by a fixed fill value (e.g., -9999).

- Computes masks from representative variables:
    (y, x)         -> rep_var_yx
    (y, xStagger)  -> rep_var_y_xs
    (yStagger, x)  -> rep_var_ys_x

- "Validity" for 3D variables means **all vertical levels** are non-NaN at a horizontal cell.
- Finds the largest all-valid rectangle in each 2D mask using the histogram/stack maximal rectangle algorithm.
- Reconciles shared dims so the saved NetCDF has consistent x/y dimension sizes.
- Saves:
    1) "<input_basename>_clipped.nc"      (original data; NaNs kept)
    2) "<input_basename>_clipped_filled.nc" (NaNs -> FILL_VALUE; _FillValue updated)

- Verifies (after saving) with secondary variables that the clipped extents have no NaNs
  (according to the conservative vertical policy), and then verifies the filled file has no NaNs at all.

Requirements: xarray, numpy. Executed without CLI arguments.
"""

import os
import sys
import numpy as np
import xarray as xr

# ======================
# ====== CONFIG =========
# ======================

DATASET_PATH = "/home/tfmrodge/scratch/GEMMACH_data/data/Inmap_outputs/Preproc/20260210/20260217_inmapData_GEMMACH_BASEGM_2015_017.nc"

# Representative variables per dimension signature
rep_var_yx     = "TotalPM25"   # dims like (z, y, x) or (y, x)
rep_var_y_xs   = "UAvg"        # dims like (z, y, xStagger) or (y, xStagger)
rep_var_ys_x   = "VAvg"        # dims like (z, yStagger, x) or (yStagger, x)

# Secondary check variables (post-save verification)
check_var_yx   = "Temperature"  # another (z, y, x) or (y, x)
check_var_y_xs = "UDeviation"   # (z, y, xStagger) or (y, xStagger)
check_var_ys_x = "VDeviation"   # (z, yStagger, x) or (yStagger, x)

# Vertical reduction policy: a horizontal cell is valid only if **all** vertical levels are non-NaN.
VERTICAL_ALL = True  # conservative

# Load behavior
LOAD_ALL_INTO_MEMORY = True  # You asked to load fully; set False for lazy/dask.

# Filled output behavior
WRITE_FILLED_OUTPUT = True
FILL_VALUE = -9999.0
# NetCDF write format (keep NETCDF3_64BIT if you want maximum compatibility with older tools)
NETCDF_FORMAT = "NETCDF3_64BIT"  # or None (default NetCDF4) for better compression/attributes

# Optional compression for NetCDF4 (ignored by NetCDF3)
COMPRESS = False
COMP_LEVEL = 1  # 0..9

# ======================
# ====== HELPERS =======
# ======================

def ensure_var_exists(ds: xr.Dataset, varname: str):
    if varname not in ds:
        raise KeyError(f"Variable '{varname}' not found in dataset.")
    return ds[varname]

def build_2d_valid_mask(arr: xr.DataArray, y_name: str, x_name: str) -> np.ndarray:
    """
    Create a 2D boolean mask (True=valid cell) for dims (y_name, x_name) from
    a representative variable that may have additional (vertical) dims.

    Validity policy:
      - If VERTICAL_ALL=True and the array has a vertical dim ('z' or 'zStagger'),
        a cell is valid only if **all** levels along the vertical are non-NaN.
      - If no vertical dim, validity is simply arr.notnull() at 2D.
    """
    dims = set(arr.dims)
    vertical_dims = []
    for cand in ("z", "zStagger"):
        if cand in dims:
            vertical_dims.append(cand)

    if y_name not in dims or x_name not in dims:
        raise ValueError(
            f"Array '{arr.name}' must include dims '{y_name}' and '{x_name}', "
            f"but has dims {arr.dims}"
        )

    if vertical_dims:
        valid = arr.notnull().all(dim=tuple(vertical_dims)) if VERTICAL_ALL else arr.notnull().any(dim=tuple(vertical_dims))
    else:
        valid = arr.notnull()

    valid2d = valid.transpose(y_name, x_name).values
    if valid2d.dtype != bool:
        valid2d = valid2d.astype(bool)
    return valid2d

def largest_rectangle_in_histogram(heights: np.ndarray):
    """
    1D maximal rectangle under histogram (monotonic stack).
    Returns left, right_exclusive, height, area.
    """
    stack = []
    max_area = 0
    best_left = 0
    best_right = 0
    best_height = 0

    for i in range(len(heights) + 1):
        h = heights[i] if i < len(heights) else 0
        start = i
        while stack and stack[-1][1] > h:
            idx, height = stack.pop()
            area = height * (i - idx)
            if area > max_area:
                max_area = area
                best_left = idx
                best_right = i
                best_height = height
            start = idx
        stack.append((start, h))
    return best_left, best_right, best_height, max_area

def largest_all_valid_rectangle(mask2d: np.ndarray):
    """
    Find the largest area rectangle of True in a 2D boolean mask.
    Returns (y0, y1_exclusive, x0, x1_exclusive, area).
    """
    if mask2d.ndim != 2:
        raise ValueError("mask2d must be a 2D boolean array.")
    nrows, ncols = mask2d.shape
    heights = np.zeros(ncols, dtype=np.int64)
    best = (0, 0, 0, 0, 0)  # y0, y1, x0, x1, area

    for r in range(nrows):
        row = mask2d[r]
        heights = heights + 1
        heights = np.where(row, heights, 0)

        left, right, height, area = largest_rectangle_in_histogram(heights)
        if area > best[4]:
            y1 = r + 1
            y0 = y1 - height
            x0 = left
            x1 = right
            best = (y0, y1, x0, x1, area)

    return best

def constrain_mask(mask: np.ndarray, y_slice=None, x_slice=None):
    """
    Return a view of mask constrained in y/x if slices are provided.
    slices are tuples (start, end_exclusive).
    """
    m = mask
    if y_slice is not None:
        ys, ye = y_slice
        m = m[ys:ye, :]
    if x_slice is not None:
        xs, xe = x_slice
        m = m[:, xs:xe]
    return m

def nonempty_intersection(a0, a1, b0, b1, label="interval"):
    """Return intersection [max(a0,b0), min(a1,b1)) and assert it's non-empty."""
    s = max(a0, b0)
    e = min(a1, b1)
    if e <= s:
        raise RuntimeError(f"No overlap for {label}: [{a0},{a1}) ∩ [{b0},{b1}) is empty.")
    return s, e

def ds_apply_clipping(ds: xr.Dataset, clip_slices: dict) -> xr.Dataset:
    """
    clip_slices maps dim_name -> (start, end_exclusive).
    Applies isel on all dims present; leaves other dims untouched.
    """
    indexers = {}
    for dim_name, (s, e) in clip_slices.items():
        indexers[dim_name] = slice(s, e)
    return ds.isel(indexers)

# ============================
# ====== MAIN WORKFLOW =======
# ============================

def main():
    print("Opening dataset...")
    ds = xr.open_dataset(DATASET_PATH, decode_cf=True)

    if LOAD_ALL_INTO_MEMORY:
        print("Loading dataset fully into memory (as requested)...")
        ds.load()

    # Original globals (for updated x0,y0 calc)
    x0 = float(ds.attrs.get("x0", 0.0))
    y0 = float(ds.attrs.get("y0", 0.0))
    dx = float(ds.attrs.get("dx", 1.0))
    dy = float(ds.attrs.get("dy", 1.0))

    # Build masks from representatives
    print("Building 2D validity masks from representative variables...")
    var_yx   = ensure_var_exists(ds, rep_var_yx)
    var_y_xs = ensure_var_exists(ds, rep_var_y_xs)
    var_ys_x = ensure_var_exists(ds, rep_var_ys_x)

    mask_yx   = build_2d_valid_mask(var_yx,   "y",        "x")
    mask_y_xs = build_2d_valid_mask(var_y_xs, "y",        "xStagger")
    mask_ys_x = build_2d_valid_mask(var_ys_x, "yStagger", "x")

    # First pass: Largest rectangles independently
    print("Finding largest all-valid rectangles (independent pass)...")
    y0a, y1a, x0a, x1a, area_a = largest_all_valid_rectangle(mask_yx)
    y0b, y1b, xs0, xs1, area_b = largest_all_valid_rectangle(mask_y_xs)
    ys0, ys1, x0c, x1c, area_c = largest_all_valid_rectangle(mask_ys_x)

    # Guard: all groups must have a valid rectangle
    if area_a == 0:
        raise RuntimeError("No all-valid rectangle found for (y,x) representative.")
    if area_b == 0:
        raise RuntimeError("No all-valid rectangle found for (y,xStagger) representative.")
    if area_c == 0:
        raise RuntimeError("No all-valid rectangle found for (yStagger,x) representative.")

    print("\nIndependent rectangles:")
    print(f"(y,x):        y=[{y0a}:{y1a})  x=[{x0a}:{x1a})  size=({y1a-y0a},{x1a-x0a})")
    print(f"(y,xStagger): y=[{y0b}:{y1b})  xS=[{xs0}:{xs1}) size=({y1b-y0b},{xs1-xs0})")
    print(f"(yStagger,x): yS=[{ys0}:{ys1}) x=[{x0c}:{x1c})  size=({ys1-ys0},{x1c-x0c})")

    # Reconcile 'y' across (y,x) and (y,xStagger)
    y_final0, y_final1 = nonempty_intersection(y0a, y1a, y0b, y1b, label="y dimension")
    print(f"\nAfter reconciling y: y=[{y_final0}:{y_final1})  (height={y_final1-y_final0})")

    # Recompute within constrained y for those two groups
    yx_mask_yc   = constrain_mask(mask_yx,   y_slice=(y_final0, y_final1))
    yxs_mask_yc  = constrain_mask(mask_y_xs, y_slice=(y_final0, y_final1))

    y0a2r, y1a2r, x0a2, x1a2, area_a2 = largest_all_valid_rectangle(yx_mask_yc)
    y0b2r, y1b2r, xs0_2, xs1_2, area_b2 = largest_all_valid_rectangle(yxs_mask_yc)

    # Guard: must remain valid after y reconciliation
    if (y1a2r - y0a2r) == 0 or (x1a2 - x0a2) == 0:
        raise RuntimeError("After y reconciliation, (y,x) rectangle is empty.")
    if (y1b2r - y0b2r) == 0 or (xs1_2 - xs0_2) == 0:
        raise RuntimeError("After y reconciliation, (y,xStagger) rectangle is empty.")

    # Convert relative y back to absolute indices
    y0a2 = y_final0 + y0a2r
    y1a2 = y_final0 + y1a2r
    # (We intentionally do not carry y0b2/y1b2; we enforce a single y-slice.)

    # Now reconcile 'x' across (y,x) and (yStagger,x)
    x_final0, x_final1 = nonempty_intersection(x0a2, x1a2, x0c, x1c, label="x dimension")
    print(f"After reconciling x: x=[{x_final0}:{x_final1})  (width={x_final1-x_final0})")

    # Recompute within constrained (y, x) to finalize for (y,x)
    yx_mask_yx_c = constrain_mask(mask_yx, y_slice=(y_final0, y_final1), x_slice=(x_final0, x_final1))

    if not yx_mask_yx_c.any():
        raise RuntimeError("Final (y,x) constrained region is empty.")

    if not np.all(yx_mask_yx_c):
        # Shrink to largest rectangle within constrained area
        y0a3r, y1a3r, x0a3r, x1a3r, _ = largest_all_valid_rectangle(yx_mask_yx_c)
        if (y1a3r - y0a3r) == 0 or (x1a3r - x0a3r) == 0:
            raise RuntimeError("Unable to find a valid rectangle after reconciling 'y' and 'x'.")
        x0_base = x_final0
        y0a2 = y_final0 + y0a3r
        y1a2 = y_final0 + y1a3r
        x_final0 = x0_base + x0a3r
        x_final1 = x0_base + x1a3r

    # For (yStagger,x): finalize with x fixed; y indices are absolute since we didn't slice y
    ysx_mask_xc = constrain_mask(mask_ys_x, x_slice=(x_final0, x_final1))
    ys0_2r, ys1_2r, x0c2_r, x1c2_r, _ = largest_all_valid_rectangle(ysx_mask_xc)
    ys0_2 = ys0_2r
    ys1_2 = ys1_2r

    # Final slices (exclusive end)
    y_slice        = (y0a2, y1a2)
    x_slice        = (x_final0, x_final1)
    xStagger_slice = (xs0_2, xs1_2)
    yStagger_slice = (ys0_2, ys1_2)

    print("\n=== Final Slices (exclusive end) ===")
    print(f"(y,x):        y=[{y_slice[0]}:{y_slice[1]}]  x=[{x_slice[0]}:{x_slice[1]}]  size=({y_slice[1]-y_slice[0]}, {x_slice[1]-x_slice[0]})")
    print(f"(y,xStagger): y=[{y_slice[0]}:{y_slice[1]}]  xS=[{xStagger_slice[0]}:{xStagger_slice[1]}]  size=({y_slice[1]-y_slice[0]}, {xStagger_slice[1]-xStagger_slice[0]})")
    print(f"(yStagger,x): yS=[{yStagger_slice[0]}:{yStagger_slice[1]}]  x=[{x_slice[0]}:{x_slice[1]}]  size=({yStagger_slice[1]-yStagger_slice[0]}, {x_slice[1]-x_slice[0]})")

    # Apply clipping to the whole dataset
    clip_slices = {}
    for dim_name, sl in [("y", y_slice), ("x", x_slice), ("xStagger", xStagger_slice), ("yStagger", yStagger_slice)]:
        if dim_name in ds.dims:
            clip_slices[dim_name] = sl

    print("\nClipping dataset...")
    ds_clipped = ds_apply_clipping(ds, clip_slices)

    # Sanity: slice widths must match resulting dims
    calc_ny = y_slice[1] - y_slice[0]
    calc_nx = x_slice[1] - x_slice[0]
    actual_ny = int(ds_clipped.dims.get("y", 0))
    actual_nx = int(ds_clipped.dims.get("x", 0))
    assert calc_ny == actual_ny, f"Computed ny={calc_ny} but ds_clipped has y={actual_ny}"
    assert calc_nx == actual_nx, f"Computed nx={calc_nx} but ds_clipped has x={actual_nx}"

    # Update globals for base (y,x) grid only
    new_x0 = float(x0 + x_slice[0] * dx)
    new_y0 = float(y0 + y_slice[0] * dy)
    ny_new = actual_ny
    nx_new = actual_nx

    ds_clipped.attrs["x0"] = new_x0
    ds_clipped.attrs["y0"] = new_y0
    ds_clipped.attrs["nx"] = nx_new
    ds_clipped.attrs["ny"] = ny_new

    print("\n=== Updated global attributes ===")
    print(f"x0: {new_x0}")
    print(f"y0: {new_y0}")
    print(f"nx: {nx_new}")
    print(f"ny: {ny_new}")

    # Save (original clipped with NaNs preserved)
    base, ext = os.path.splitext(DATASET_PATH)
    if ext == "":
        ext = ".nc"
    out_path = f"{base}_clipped{ext}"
    print(f"\nSaving clipped dataset to: {out_path}")

    if COMPRESS and NETCDF_FORMAT is None:
        encoding = {name: {"zlib": True, "complevel": COMP_LEVEL} for name in ds_clipped.data_vars}
    else:
        encoding = {}

    ds_clipped.to_netcdf(
        out_path,
        encoding=encoding,
        format=NETCDF_FORMAT if NETCDF_FORMAT else None
    )

    # ============================
    # Post-save validation checks
    # ============================
    print("\nRe-opening clipped dataset for validation...")
    ds2 = xr.open_dataset(out_path, decode_cf=True)
    if LOAD_ALL_INTO_MEMORY:
        ds2.load()

    def check_no_nans_in_region(dscheck: xr.Dataset, varname, dims):
        v = ensure_var_exists(dscheck, varname)
        y_like, x_like = dims
        v2 = v
        dims_set = set(v.dims)
        for vd in ("z", "zStagger"):
            if vd in dims_set:
                v2 = v2.notnull().all(dim=vd)
        if v2 is v:
            v2 = v2.notnull()
        if y_like in v2.dims and x_like in v2.dims:
            v2 = v2.transpose(y_like, x_like)
            return bool(v2.values.all())
        return True

    ok_yx   = check_no_nans_in_region(ds2, check_var_yx,   ("y", "x"))         if check_var_yx   in ds2 else True
    ok_yxs  = check_no_nans_in_region(ds2, check_var_y_xs, ("y", "xStagger"))  if check_var_y_xs in ds2 else True
    ok_ysx  = check_no_nans_in_region(ds2, check_var_ys_x, ("yStagger", "x"))  if check_var_ys_x in ds2 else True

    print("\n=== Validation (secondary variables) ===")
    print(f"{check_var_yx:>15} (y,x)        no-NaN in clipped region: {ok_yx}")
    print(f"{check_var_y_xs:>15} (y,xStagger) no-NaN in clipped region: {ok_yxs}")
    print(f"{check_var_ys_x:>15} (yStagger,x) no-NaN in clipped region: {ok_ysx}")

    if not (ok_yx and ok_yxs and ok_ysx):
        print("\nWARNING: One or more secondary checks found NaNs within the clipped region.")
        print("Consider different representatives or further constraints.")

    # ============================
    # Save a filled (no-NaN) copy
    # ============================
    if WRITE_FILLED_OUTPUT:
        out_path_filled = f"{base}_clipped_filled{ext}"
        print(f"\nCreating filled dataset (NaNs -> {FILL_VALUE}) and saving to: {out_path_filled}")

        # Create a shallow dataset with variables replaced one-by-one to avoid doubling memory too much
        ds_filled = xr.Dataset(coords=ds_clipped.coords, attrs=dict(ds_clipped.attrs))
        fill_encoding = {}

        for vname, da in ds_clipped.data_vars.items():
            # Only replace for float variables; for others, pass through
            if np.issubdtype(da.dtype, np.floating):
                filled_da = da.where(da.notnull(), other=FILL_VALUE)
                # Ensure dtype stays floating
                if not np.issubdtype(filled_da.dtype, np.floating):
                    filled_da = filled_da.astype("float32")
                # Update attributes
                new_attrs = dict(da.attrs)
                new_attrs["_FillValue"] = float(FILL_VALUE)
                new_attrs["missing_value"] = float(FILL_VALUE)
                filled_da.attrs = new_attrs
                ds_filled[vname] = filled_da
                # Ensure writer uses the same fill value on-disk
                fill_encoding[vname] = {"_FillValue": float(FILL_VALUE)}
                if COMPRESS and NETCDF_FORMAT is None:
                    fill_encoding[vname].update({"zlib": True, "complevel": COMP_LEVEL})
            else:
                # Keep as-is (coords/no-float data-vars)
                ds_filled[vname] = da
                if COMPRESS and NETCDF_FORMAT is None:
                    fill_encoding[vname] = {"zlib": True, "complevel": COMP_LEVEL}

        # Save
        ds_filled.to_netcdf(
            out_path_filled,
            encoding=fill_encoding,
            format=NETCDF_FORMAT if NETCDF_FORMAT else None
        )

        # Verify no NaNs remain in filled file (data variables only)
        ds3 = xr.open_dataset(out_path_filled, decode_cf=True)
        if LOAD_ALL_INTO_MEMORY:
            ds3.load()
        has_any_nan = False
        for vname, da in ds3.data_vars.items():
            if np.issubdtype(da.dtype, np.floating):
                if np.isnan(da.values).any():
                    print(f"WARNING: Filled file still has NaNs in variable '{vname}'")
                    has_any_nan = True
        print(f"Filled file contains NaNs? {has_any_nan}")

    print("\nDone.")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print("ERROR:", e)
        sys.exit(1)