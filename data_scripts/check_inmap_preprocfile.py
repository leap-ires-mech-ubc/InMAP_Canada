#check_inmap_preprocfile.py

# compare_inmapdata.py
import xarray as xr, numpy as np
import pandas as pd
from collections import defaultdict
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

OLD = '/home/tfmrodge/scratch/GEMMACH_data/data/Inmap_outputs/Preproc/20260210/inmapData_GEMMACH_Aug1-15.nc' #"/home/tfmrodge/scratch/GEMMACH_data/inmapData_GEMMACH_BASEGM_2015_017_complete.nc"
NEW = '/home/tfmrodge/scratch/GEMMACH_data/data/Inmap_outputs/Preproc/20260210/inmapData_GEMMACH_June16-30.nc' #"/home/tfmrodge/scratch/GEMMACH_data/data/Inmap_outputs/Preproc/20260210/20260217_inmapData_GEMMACH_BASEGM_2015_017_complete.nc"
print(OLD)
print(NEW)
# # If your files are large, you can remove .load() to keep things lazy.

# # Variables to inspect (present in your headers)
# CRIT_YX   = ["TotalPM25","Kxxyy","Kzz","Temperature","ParticleDryDep",
#              "NH3DryDep","NOxDryDep","SO2DryDep","VOCDryDep",
#              "SPartitioning","NHPartitioning","NOPartitioning",
#              "aSOA","bSOA","aVOC","bVOC","gNH","gNO","gS","pNH","pNO","pS","WindSpeed",
#              "WindSpeedInverse","WindSpeedMinusOnePointFour","WindSpeedMinusThird","alt","S1","Sclass","SO2oxidation"]
# CRIT_Y_XS = ["UAvg","UDeviation"]
# CRIT_YS_X = ["VAvg","VDeviation"]
# CRIT_W    = ["WAvg","LayerHeights"]
# EXISTENCE_ONLY = ["Dz"]  # will check Dz through consistency with LayerHeights

# CHECK_GROUPS = [("yx", CRIT_YX, "y", "x"),
#                 ("y_xs", CRIT_Y_XS, "y", "xStagger"),
#                 ("ys_x", CRIT_YS_X, "yStagger", "x"),
#                 ("w", CRIT_W, "y", "x")]

# def load_ds(p):
#     ds = xr.open_dataset(p, decode_cf=True)
#     ds.load()
#     return ds

# def dims_summary(ds):
#     return {k:int(v) for k,v in ds.dims.items()}

# def attrs_summary(ds):
#     keys = ["x0","y0","dx","dy","nx","ny","data_version"]
#     return {k: ds.attrs.get(k, None) for k in keys}

# def has_var(ds, v):
#     return v in ds.data_vars

# def nan_mask_summary(da, ydim, xdim, pad=1):
#     # Interior finite mask (requires all finite across vertical dims).
#     v_dims = [d for d in da.dims if d in ("z","zStagger")]
#     a = da
#     if v_dims:
#         a = xr.apply_ufunc(np.isfinite, a).all(dim=tuple(v_dims))
#     else:
#         a = xr.apply_ufunc(np.isfinite, a)
#     a2 = a
#     # Overall finite%
#     finite_pct = float(a2.mean().item())*100.0
#     # Interior finite%
#     if ydim in a2.dims and xdim in a2.dims and a2.sizes[ydim] > 2 and a2.sizes[xdim] > 2:
#         a_int = a2.isel({ydim: slice(pad, a2.sizes[ydim]-pad),
#                          xdim: slice(pad, a2.sizes[xdim]-pad)})
#         finite_pct_int = float(a_int.mean().item())*100.0
#     else:
#         finite_pct_int = np.nan
#     # Border ring NaN fraction
#     if ydim in a2.dims and xdim in a2.dims and a2.sizes[ydim] > 2 and a2.sizes[xdim] > 2:
#         top    = a2.isel({ydim: slice(0,1)}).mean()
#         bottom = a2.isel({ydim: slice(-1,None)}).mean()
#         left   = a2.isel({xdim: slice(0,1)}).mean()
#         right  = a2.isel({xdim: slice(-1,None)}).mean()
#         # Convert to NaN fraction (1 - finite)
#         border_nan_frac = float(1.0 - np.mean([float(top), float(bottom), float(left), float(right)]))
#     else:
#         border_nan_frac = np.nan
#     return finite_pct, finite_pct_int, border_nan_frac

# def stats_summary(da):
#     a = da
#     # reduce vertical dims when needed for min/max (nan-safe)
#     while any(d in a.dims for d in ("z","zStagger")):
#         d = "z" if "z" in a.dims else "zStagger"
#         a = a.min(dim=d, skipna=True)
#     mn = float(a.min(skipna=True).item()) if a.size else np.nan
#     mx = float(a.max(skipna=True).item()) if a.size else np.nan
#     all_zero = bool((np.nan_to_num(a.values)==0).all()) if a.size else False
#     return mn, mx, all_zero

# def compare_two(ds_old, ds_new):
#     report_rows = []
#     issues = []

#     # 0) Heads-up summaries
#     print("OLD dims:", dims_summary(ds_old))
#     print("NEW dims:", dims_summary(ds_new))
#     print("OLD attrs:", attrs_summary(ds_old))
#     print("NEW attrs:", attrs_summary(ds_new))

#     # Stagger sanity
#     for tag, ds in [("OLD", ds_old), ("NEW", ds_new)]:
#         xd = ds.dims.get("x",0); xsd = ds.dims.get("xStagger",0)
#         yd = ds.dims.get("y",0); ysd = ds.dims.get("yStagger",0)
#         zd = ds.dims.get("z",0); zsd = ds.dims.get("zStagger",0)
#         if xsd != xd+1: issues.append(f"{tag}: xStagger ({xsd}) != x+1 ({xd+1})")
#         if ysd != yd+1: issues.append(f"{tag}: yStagger ({ysd}) != y+1 ({yd+1})")
#         if zsd != zd+1: issues.append(f"{tag}: zStagger ({zsd}) != z+1 ({zd+1})")

#     # 1) Variable presence / dtype / units mismatch
#     all_names = set()
#     for _, names, _, _ in CHECK_GROUPS:
#         all_names |= set(names)
#     all_names |= set(EXISTENCE_ONLY)

#     for v in sorted(all_names):
#         exists_old = has_var(ds_old, v)
#         exists_new = has_var(ds_new, v)
#         if not exists_old or not exists_new:
#             issues.append(f"Var presence mismatch: {v}: old={exists_old}, new={exists_new}")
#             continue
#         t_old = str(ds_old[v].dtype)
#         t_new = str(ds_new[v].dtype)
#         if t_old != t_new:
#             issues.append(f"dtype differs {v}: old={t_old}, new={t_new}")
#         u_old = ds_old[v].attrs.get("units", None)
#         u_new = ds_new[v].attrs.get("units", None)
#         if u_old != u_new:
#             issues.append(f"units differ {v}: old={u_old}, new={u_new}")

#     # 2) Finite/NaN patterns & ranges
#     def add_row(tag, v, scope, fin, fin_int, border_nan, mn, mx, all_zero, shape):
#         report_rows.append(dict(
#             file=tag, var=v, scope=scope, shape=str(shape),
#             finite_pct=fin, finite_pct_interior=fin_int,
#             border_nan_frac=border_nan, vmin=mn, vmax=mx, all_zero=all_zero
#         ))

#     for tag, ds in [("old", ds_old), ("new", ds_new)]:
#         for scope, names, ydim, xdim in CHECK_GROUPS:
#             for v in names:
#                 if not has_var(ds, v) or ydim not in ds[v].dims or xdim not in ds[v].dims:
#                     continue
#                 fin, fin_int, bnan = nan_mask_summary(ds[v], ydim, xdim, pad=1)
#                 mn, mx, all_zero = stats_summary(ds[v])
#                 add_row(tag, v, scope, fin, fin_int, bnan, mn, mx, all_zero, ds[v].shape)

#     # 3) Dz vs LayerHeights consistency (optional but very telling)
#     for tag, ds in [("old", ds_old), ("new", ds_new)]:
#         if "LayerHeights" in ds and "Dz" in ds:
#             dz_from_lh = ds["LayerHeights"].diff("zStagger")
#             dz_from_lh = dz_from_lh.rename({"zStagger":"z"})
#             # Max abs diff (nan-safe)
#             diff = (dz_from_lh - ds["Dz"]).astype("float64")
#             mad = float(np.nanmax(np.abs(diff.values)))
#             if mad > 1e-6:  # tolerance can be adjusted
#                 issues.append(f"Dz mismatch in {tag}: max|Dz - diff(LayerHeights)| = {mad}")

#     # Final report
#     df = pd.DataFrame(report_rows)
#     return df, issues

# if __name__ == "__main__":
#     ds_old = load_ds(OLD)
#     ds_new = load_ds(NEW)
#     df, issues = compare_two(ds_old, ds_new)

#     # Sort by variables with interior finite% < 100 in the new file (most suspicious)
#     suspicious = df[(df["file"]=="new") & (df["finite_pct_interior"] < 100)]
#     print("\n=== Suspicious (new file interior finite% < 100%) ===")
#     print(suspicious.sort_values(["finite_pct_interior","var"]).to_string(index=False))

#     print("\n=== Differences & Issues ===")
#     for s in issues:
#         print("-", s)

#     # Optional: write CSV for a full comparison table
#     out_csv = "inmap_preproc_compare_report.csv"
#     df.to_csv(out_csv, index=False)
#     print("\nWrote:", out_csv)
# import xarray as xr
# import numpy as np

old = xr.open_dataset(OLD)
new = xr.open_dataset(NEW)

def nz(da):  # nan-safe min/max
    return float(da.min(skipna=True).item()), float(da.max(skipna=True).item())

def count_where(cond):
    # cond is a boolean DataArray
    return int(cond.sum().item())

def scan(ds, tag):
    print(f"\n=== {tag} ===")
    # 1) Dz: non-positive, tiny cells
    if "Dz" in ds:
        mn, mx = nz(ds["Dz"])
        bad_le0 = count_where(ds["Dz"] <= 0)
        tiny = count_where(ds["Dz"] < 0.05)  # threshold can be tuned
        print(f"Dz min/max = {mn:.6g}/{mx:.6g}  |  Dz<=0 cells: {bad_le0}  |  Dz<0.05m cells: {tiny}")

    # 2) LayerHeights monotonic & positive increments
    if "LayerHeights" in ds:
        dz_from_lh = ds["LayerHeights"].diff("zStagger").rename({"zStagger":"z"})
        lh_bad_le0 = count_where(dz_from_lh <= 0)
        print(f"diff(LayerHeights)<=0 cells across z: {lh_bad_le0}")

    # 3) Kzz / Kxxyy non-negative
    for v in ["Kzz","Kxxyy"]:
        if v in ds:
            mn, mx = nz(ds[v])
            neg = count_where(ds[v] < 0)
            print(f"{v} min/max = {mn:.6g}/{mx:.6g}  |  negatives: {neg}")

    # 4) PBL height strictly positive (if present)
    if "Pblh" in ds:
        mn, mx = nz(ds["Pblh"])
        bad = count_where(ds["Pblh"] <= 0)
        print(f"Pblh min/max = {mn:.6g}/{mx:.6g}  |  Pblh<=0 cells: {bad}")

    # 5) Inverse density positive
    if "alt" in ds:
        mn, mx = nz(ds["alt"])
        bad = count_where(ds["alt"] <= 0)
        print(f"alt (1/rho) min/max = {mn:.6g}/{mx:.6g}  |  alt<=0 cells: {bad}")

    # 6) Wind fields finite and reasonable
    for v in ["WindSpeed","WindSpeedInverse","WindSpeedMinusOnePointFour","WindSpeedMinusThird","UAvg","VAvg","WAvg"]:
        if v in ds:
            mn, mx = nz(ds[v])
            any_nan = bool(np.isnan(ds[v].values).any())
            print(f"{v} min/max = {mn:.6g}/{mx:.6g}  |  has NaN in file? {any_nan}")

scan(old, "OLD")
scan(new, "NEW")
# import xarray as xr
# import numpy as np

# ds = xr.open_dataset(NEW)

# # Boolean mask of problematic cells across all z
# bad_any = (ds["Dz"] <= 0).any(dim="z")  # (y,x): True means at least one layer has Dz<=0

# # Optionally: count # of bad layers per column
# bad_count = (ds["Dz"] <= 0).sum(dim="z")  # (y,x) integer

# # Write to GeoTIFF for QGIS
# bad_any = bad_any.astype("uint8")  # convert booleans to 0/1
# #bad_any.rio.set_crs("EPSG:32610", inplace=True) # or your LCC CRS
# bad_any.rio.to_raster("Dz_le0_anylayer.tif")

# #bad_count.rio.set_crs("EPSG:32610", inplace=True)
# bad_count.rio.to_raster("/home/tfmrodge/scratch/GEMMACH_data/data/Inmap_outputs/Preproc/20260210/Dz_le0_count_old.tif")