#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
if os.environ.get("ENABLE_DEBUGPY") == "1":
    import debugpy
    debugpy.listen(("0.0.0.0", 5678))
    print("🔧 debugpy: listening on 5678 (waiting for VS Code to attach...)")
    debugpy.wait_for_client()
    print("🔧 debugpy: client attached, continuing")
import sys
import numpy as np
import xarray as xr
import rioxarray
from rasterio.transform import from_origin
from datetime import datetime

def now():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def ip(msg):
    print(f"[{now()}] {msg}", flush=True)

# ---- GEMMACH scenario table (same as your main pipeline) ----
SCENARIOS = {
    1: ['BASEGM_2015_017','BASEGM_2015_017'],
    2: ['BAU_2020_E108','BAU_2020_105'],
    3: ['BAU_2025_E107','BAU_2025_106'],
    4: ['BAU_2030_E107','BAU_2030_106'],
    5: ['BAU_2035_E107','BAU_2035_106'],
    6: ['BPT_2015_E002','BPT_2015_002'],
    7: ['BPT_2015_E006','BPT_2015_006'],
    8: ['BPT_2015_E015','BPT_2015_015'],
    9: ['COVID_2020_E003','COVID_2020_003'],
    10:['BPT_2015_E016','BPT_2015_E016']
}

# ---- Your LCC grid spec ----
xfirst = -4184312.05377675
xinc   =  10002.6608465054
yfirst = -2029822.82977676
yinc   =  10002.6608464873
xsize  =  744
ysize  =  669

# CRS for GEMMACH (Lambert Conformal Conic)
proj_str = (
    "+proj=lcc +lat_1=50 +lat_2=70 +lat_0=40 "
    "+lon_0=-96 +x_0=0 +y_0=0 "
    "+a=6378137 +rf=298.257222101 +units=m +no_defs"
)

def write_lcc_georef(nc_in, nc_out):
    ip(f"Opening {nc_in}")
    ds = xr.open_dataset(nc_in, engine="h5netcdf")

    # Build x/y coordinates
    x = xfirst + xinc * np.arange(xsize)
    y = yfirst + yinc * np.arange(ysize)

    ds = ds.assign_coords(x=("x", x))
    ds = ds.assign_coords(y=("y", y))

    # Affine transform — y_up requires negative pixel height
    y_top = yfirst + yinc * (ysize - 1)
    transform = from_origin(xfirst, y_top, xinc, -abs(yinc))

    # Write CRS + transform to each spatial variable
    for var in ds.data_vars:
        da = ds[var]
        if ("y" in da.dims) and ("x" in da.dims):
            da = da.rio.write_crs(proj_str, inplace=False, grid_mapping_name="crs")
            da = da.rio.write_transform(transform, inplace=False)
            ds[var] = da

    ip(f"Writing output → {nc_out}")
    encoding = {v: {"zlib": True, "complevel": 2} for v in ds.data_vars}
    ds.to_netcdf(nc_out, engine="h5netcdf", encoding=encoding)
    ip(f"Done writing {nc_out}")


def main():
    if len(sys.argv) < 2:
        print("Usage: python gemmach_reproject_lcc.py <scenario_id>")
        sys.exit(1)

    sID = int(sys.argv[1])
    emissions, scenario = SCENARIOS[sID]

    ip(f"=== Running LCC projection fix for scenario {scenario} ===")

    # paths #/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/scenario_sums2/20260120_BAU_2020_105_annual_mean_LCC.nc
    src_dir = "/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/scenario_sums2"
    dst_dir = "/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/scenario_sums_LCC"

    os.makedirs(dst_dir, exist_ok=True)

    nc_in  = f"{src_dir}/20260120_{scenario}_annual_mean_LCC.nc"
    nc_out = f"{dst_dir}/20260120_{scenario}_annual_mean_LCC.nc"
    write_lcc_georef(nc_in, nc_out)

    ip("All done.")


if __name__ == "__main__":
    main()
