import xarray as xr
import geopandas as gpd
from rasterstats import zonal_stats
import rioxarray as rio
import pandas as pd
import numpy as np
import sys
from joblib import Parallel, delayed
import multiprocessing
#module load python/3.10.2 scipy-stack mpi4py
#source inmapenv/bin/activate

#salloc --time=0:59:59 --mem-per-cpu=3G --ntasks=15 --account=def-rscholes
#pip install --no-index --upgrade pip 
#pip install --no-index xarray geopandas

ncpus = multiprocessing.cpu_count()
def compute_stats_for_geom(geom, image, affine, stats, crs):
    # Wrap geometry in a GeoDataFrame with correct CRS
    single_shape = gpd.GeoDataFrame(geometry=[geom], crs=crs)
    result = zonal_stats(
                single_shape,
                image,
                affine=affine,
                stats=stats,
                nodata=nodata,       #
                all_touched=False,    # or True if needed for your use case
                raster_out=False
            )[0]

    #result = zonal_stats(single_shape, image, stats=stats)[0]#No affine since same pcr
    return result

def ZonalStatsParallel(shape, raster, rastcols, stats=['sum'], n_jobs=4):
    # shape - geopandas shapefile object
    # image- xarray dataset
    # affine - raster affine (xrdataset.rio.transform() to get)
    # stats - stats as list, i.e. ['min', 'mean', 'max'] ; ['sum']
    # the result is final_gdf as GeoDataFrame
    #affine = raster.rio.transform()
    geometries = shape.geometry
    crs = shape.crs  # Automatically pull CRS from shapefile
    df_concat = shape.copy()
    for rc in rastcols: 
        data_array = raster[rc]
    # Reduce to 2D
    if all(dim in data_array.dims for dim in ['time', 'level1']):
        data_array = data_array.isel(time=0, level1=0)
    elif 'time' in data_array.dims:
        data_array = data_array.isel(time=0)
    elif 'level1' in data_array.dims:
        data_array = data_array.isel(level1=0)

    # Tell rioxarray which are spatial dims
    data_array = data_array.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=False)

    # Ensure x is ascending (left→right) and y is descending (top→bottom):
    if 'x' in data_array.coords and data_array.x[0] > data_array.x[-1]:
        data_array = data_array.sortby('x', ascending=True)
    if 'y' in data_array.coords and data_array.y[0] < data_array.y[-1]:
        # If y increases northward, flip to top→bottom
        data_array = data_array.sortby('y', ascending=False)
    # Recalculate transform after sorting
    affine = data_array.rio.transform(recalc=True)
    # Sanity checks on transform signs (north-up raster)
    if affine.a <= 0:
        raise ValueError(f"Unexpected pixel width in transform (a={affine.a}).")
    if affine.e >= 0:
        # If this triggers, y is still "upside down"
        raise ValueError(f"Unexpected pixel height in transform (e={affine.e}). Expected negative.")

    # Extract the 2D numpy array in (rows, cols) = (y, x) order
    image = data_array.values
    if image.ndim != 2:
        raise ValueError(f"Unexpected number of dims ({image.ndim}) for {rc}. Dims: {data_array.dims}")

    # Optional: check alignment between array shape and coords
    assert image.shape == (data_array.sizes['y'], data_array.sizes['x'])

    # Handle nodata (better: tell zonal_stats via nodata=…)
    nodata_value = data_array.rio.nodata
    if nodata_value is not None:
        image = np.where(image == nodata_value, np.nan, image)


        # Parallelize over geometries
        
        results = Parallel(n_jobs=ncpus)(
            delayed(compute_stats_for_geom)(geom, image, affine, stats, outshp.crs, nodata_value)
            for geom in geometries
        )


        # Convert results to DataFrame
        zonal_df = pd.DataFrame(results)
        zonal_df.columns = [f"{rc}_{s}" for s in stats]

        # Concatenate with original shape
        df_concat = pd.concat([df_concat.reset_index(drop=True), zonal_df], axis=1)
        print(f"Done {rc}")

    return gpd.GeoDataFrame(df_concat, geometry='geometry', crs=crs)

scenarios={1:['BASEGM_2015_017','BASEGM_2015_017'],2:['BAU_2020_E108','BAU_2020_105'],3:['BAU_2025_E107','BAU_2025_106'],4:['BAU_2030_E107','BAU_2030_106'],
5:['BAU_2035_E107','BAU_2035_106'],6:['BPT_2015_E002','BPT_2015_002'],7:['BPT_2015_E006','BPT_2015_006'],
8:['BPT_2015_E015','BPT_2015_015'],9:['COVID_2020_E003','COVID_2020_003'],10:['BPT_2015_E016','BPT_2015_016']}

#shp = gpd.read_file('/home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/BASEGM_2015_E010_Sumlcc.shp')
ncpath = '/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/scenario_sums_LCC/'
#ncpath = '/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/scenario_sums'
#fpath = '/home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/'
sID = int(sys.argv[1])
print(scenarios[sID][1])
scenario = scenarios[sID][1]
emissions= scenarios[sID][0]
#shp = gpd.read_file('/home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/BASEGM_2015_E010_Sumlcc.shp')
inmapoutpath = '/home/tfmrodge/scratch/GEMMACH_data/data/Inmap_outputs/Outputs'
shp = f"{inmapoutpath}/20250930_InmapOuts_{emissions}_static.shp"
#nc_area = '/home/tfmrodge/scratch/GEMMACH_data/data/Inmap_outputs/Emissions/BASEGM_2015_E010_area_lcc.nc'
#major = '/home/tfmrodge/scratch/GEMMACH_data/data/BASEGM_2015_E010/BASEGM_2015_E010_major_lcc.nc'
ncpath = "/home/tfmrodge/scratch/GEMMACH_data/data/GEMMACH_outputs/scenario_sums_LCC"
#nc_name = f"{ncpath}/{scenario}_Cav_georeferenced.nc"
nc_name = f"{ncpath}/20260120_{scenario}_annual_mean_LCC.nc"
outshp = gpd.read_file(shp)
nc_file = xr.open_dataset(nc_name, engine = 'h5netcdf')
# Assign the CRS manually
# proj_str = (
#     "+proj=lcc +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 "
#     "+x_0=0 +y_0=0 +a=6378137 +rf=298.257222101 +units=m +no_defs"
# )
# nc_file = nc_file.rio.write_crs(proj_str)
# if nc_file.rio.crs != outshp.crs:
#     print("Reprojecting raster to match shapefile CRS...")
#     nc_file = nc_file.rio.reproject(outshp.crs)
#rastcols=['BASEPM25','BASEPNO3','BASEPNH4','BasePSO4','BaseSOA','BasePrimPM25']

proj_str = (
    "+proj=lcc +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 "
    "+x_0=0 +y_0=0 +a=6378137 +rf=298.257222101 +units=m +no_defs"
)
nc_file = nc_file.rio.write_crs(proj_str)  # define raster CRS (no resampling)
if outshp.crs != nc_file.rio.crs:
    outshp = outshp.to_crs(nc_file.rio.crs)  # reproject vectors to raster CRS
RASTCOLS = ["BASEPM25", "BASEPNO3", "BASEPNH4", "BASEPSO4", "BASESOA", "BASEPRIM25"]
outshp = outshp[outshp.is_valid].copy()
outshp = outshp.explode(index_parts=False, ignore_index=True)
outarea = ZonalStatsParallel(outshp, nc_file,RASTCOLS, ['mean'],n_jobs=-1)
#outarea = outarea.to_crs('EPSG:4326')
outarea.to_crs(crs=outshp.crs,inplace=True) #Set CRS to make the same
out_path = f"{inmapoutpath}/20260128_InmapOuts_{emissions}_static_withGEMMACH.shp"
outarea.to_file(out_path)
# outmajor = ZonalStats(shp, major, 'sum')
# outmajor.to_file('/home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/BASEGM_2015_E010_major_lcc.shp')
#gpd.read_file('/home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/BAU_2020_E108_majorpts.shp')