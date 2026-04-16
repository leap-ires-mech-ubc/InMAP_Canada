import xarray as xr
import geopandas as gpd
from rasterstats import zonal_stats
from rasterio.enums import Resampling
import rioxarray as rio
import pandas as pd
import numpy as np
import sys
from joblib import Parallel, delayed
import multiprocessing

ncpus = multiprocessing.cpu_count()

def compute_stats_for_geom(geom, image, affine, stats, crs, nodata, all_touched=False):
    """Compute stats for a single geometry"""
    single_shape = gpd.GeoDataFrame(geometry=[geom], crs=crs)
    result = zonal_stats(
        single_shape,
        image,
        affine=affine,
        stats=stats,
        nodata=nodata,
        all_touched=all_touched,  # <-- now configurable
        raster_out=False
    )[0]
    return result


def ZonalStatsParallel(
    shape,
    raster,
    rastcols,
    stats=['sum'],
    n_jobs=4,
    upsample_factor=1,                 # <-- settable outside the loop
    resampling=Resampling.bilinear,    # <-- settable; bilinear recommended for intensive vars
    all_touched=False                  # <-- settable; closer to area-weighted when True
):
    geometries = shape.geometry
    crs = shape.crs
    df_concat = shape.copy()

    for rc in rastcols:
        # --- 2D extraction ---
        data_array = raster[rc]
        if all(dim in data_array.dims for dim in ['time', 'level1']):
            data_array = data_array.isel(time=0, level1=0)
        elif 'time' in data_array.dims:
            data_array = data_array.isel(time=0)
        elif 'level1' in data_array.dims:
            data_array = data_array.isel(level1=0)
        squeezable = [d for d in data_array.dims if data_array.sizes[d] == 1]
        if len(squeezable) > 0:
            data_array = data_array.squeeze(drop=True)


        # Spatial dims
        data_array = data_array.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=False)

        # Ensure x ascending (L->R) and y descending (top->bottom)
        if 'x' in data_array.coords and float(data_array.x[0]) > float(data_array.x[-1]):
            data_array = data_array.sortby('x', ascending=True)
        if 'y' in data_array.coords and float(data_array.y[0]) < float(data_array.y[-1]):
            data_array = data_array.sortby('y', ascending=False)

        # Optional upsampling (for closer to area-weighted results with all_touched)
        if upsample_factor and upsample_factor > 1:
            xres, yres = data_array.rio.resolution()  # yres typically negative for north-up
            new_res = (xres / upsample_factor, yres / upsample_factor)
            data_array = data_array.rio.reproject(
                dst_crs=data_array.rio.crs,
                resolution=new_res,
                resampling=resampling,
            )
            # Keep coord order sane after reproject
            if float(data_array.x[0]) > float(data_array.x[-1]):
                data_array = data_array.sortby('x', ascending=True)
            if float(data_array.y[0]) < float(data_array.y[-1]):
                data_array = data_array.sortby('y', ascending=False)

        # Transform after sorting / upsampling
        affine = data_array.rio.transform(recalc=True)

        # Sanity checks
        if affine.a <= 0:
            raise ValueError(f"[{rc}] Unexpected pixel width in transform (a={affine.a}).")
        if affine.e >= 0:
            raise ValueError(f"[{rc}] Unexpected pixel height in transform (e={affine.e}). Expected negative.")

        # Extract 2D image (rows, cols) = (y, x)
        image = data_array.values
        if image.ndim != 2:
            raise ValueError(f"[{rc}] Unexpected number of dims ({image.ndim}). Dims: {data_array.dims}")
        assert image.shape == (data_array.sizes['y'], data_array.sizes['x'])

        # Nodata handling: pass through to zonal_stats (NaN mask still useful)
        nodata_value = data_array.rio.nodata
        if nodata_value is not None:
            image = np.where(image == nodata_value, np.nan, image)

        # --- Parallel per geometry ---
        results = Parallel(n_jobs=n_jobs)(
            delayed(compute_stats_for_geom)(geom, image, affine, stats, crs, nodata_value, all_touched)
            for geom in geometries
        )

        # Collect
        zonal_df = pd.DataFrame(results)
        zonal_df.columns = [f"{rc}_{s}" for s in stats]
        df_concat = pd.concat([df_concat.reset_index(drop=True), zonal_df], axis=1)
        print(f"Done {rc}")

    return gpd.GeoDataFrame(df_concat, geometry='geometry', crs=crs)

# --- Main ---
scenarios = {
    1: ['BASEGM_2015_017','BASEGM_2015_017'],
    2: ['BAU_2020_E108','BAU_2020_105'],
    3: ['BAU_2025_E107','BAU_2025_106'],
    4: ['BAU_2030_E107','BAU_2030_106'],
    5: ['BAU_2035_E107','BAU_2035_106'],
    6: ['BPT_2015_E002','BPT_2015_002'],
    7: ['BPT_2015_E006','BPT_2015_006'],
    8: ['BPT_2015_E015','BPT_2015_015'],
    9: ['COVID_2020_E003','COVID_2020_003'],
    10: ['BPT_2015_E016','BPT_2015_016']
}

ncpath = "/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/20260129_scenario_sums/"
sID = int(sys.argv[1])
print(scenarios[sID][1])
scenario = scenarios[sID][1]
emissions = scenarios[sID][0]

inmapoutpath = "/home/tfmrodge/scratch/GEMMACH_data/data/Inmap_outputs/Outputs"
shp = f"{inmapoutpath}/20260310_InmapOuts_{emissions}_static.shp" #/home/tfmrodge/scratch/GEMMACH_data/data/Inmap_outputs/Outputs/20260310_InmapOuts_BASEGM_2015_017_static.shp
nc_name = f"{ncpath}/20260317_{scenario}_annual_mean_LCC.nc"

outshp = gpd.read_file(shp)
nc_file = xr.open_dataset(nc_name, engine='h5netcdf')

# Assign raster CRS and reproject vectors to raster
proj_str = (
    "+proj=lcc +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 "
    "+x_0=0 +y_0=0 +a=6378137 +rf=298.257222101 +units=m +no_defs"
)
nc_file = nc_file.rio.write_crs(proj_str)
if outshp.crs != nc_file.rio.crs:
    outshp = outshp.to_crs(nc_file.rio.crs)

# Clean shapes
outshp = outshp[outshp.is_valid].copy()
outshp = outshp.explode(index_parts=False, ignore_index=True)

RASTCOLS = ["BASEPM25", "BASEPNO3", "BASEPNH4", "BASEPSO4", "BASESOA", "BASEPRIM25"]

# Run
upsample_factor=10
outarea = ZonalStatsParallel(outshp, nc_file, RASTCOLS, stats=['mean'],
    n_jobs=-1,upsample_factor=upsample_factor,
    resampling=Resampling.bilinear,
    all_touched=True
)

# Strip any known stat suffixes: _mean, _sum, _min, _max, etc.
stats = ["mean", "sum", "min", "max", "median"]
suffixes = [f"_{s}" for s in stats]

rename_map = {
    col: col.split("_")[0]
    for col in outarea.columns
    if any(col.endswith(suf) for suf in suffixes)
}

# --- NEW: drop any existing columns that would collide after rename (case-insensitive) ---
if rename_map:
    target_names = set(rename_map.values())
    target_names_lower = {t.lower() for t in target_names}

    # columns that are going to be renamed (sources), keep them for now
    sources_being_renamed = set(rename_map.keys())

    cols_to_drop = [
        col for col in outarea.columns
        # don't drop geometry column or the sources that will be renamed
        if col not in sources_being_renamed
        and col != outarea.geometry.name
        # would collide with one of the target names (case-insensitive)
        and col.lower() in target_names_lower
    ]

    if cols_to_drop:
        outarea = outarea.drop(columns=cols_to_drop)
        print(f"Dropped pre-existing columns to avoid shapefile name collisions: {cols_to_drop}")

# Perform the rename (e.g., BASEPM25_mean -> BASEPM25)
outarea = outarea.rename(columns=rename_map)
#Checks - XX=da['BASEPM25']-da['BASESOA']-da['BASEPSO4']-da['BASEPNH4']-da['BASEPNO3']
# (XX.isnull()) | (XX.values == da['BASEPRIM25'].values) | (da['BASEPRIM25'].isnull()).sum()
#da4['BASEPNO3']=(da4['TNI1']+da4['TNO3']+da4['THN3'])*da4['RHO'];
#da5['BASEPRIM25'] = da5['BASEPM25']-da5['BaseSOA']-da5['BasePSO4']-da5['BASEPNH4']-da5['BASEPNO3']

# Already in raster CRS; 
out_path = f"{inmapoutpath}/20260310_InmapOuts_{emissions}_static_withGEMMACH.shp"
outarea.to_file(out_path)
print(f"output file as: {out_path}")