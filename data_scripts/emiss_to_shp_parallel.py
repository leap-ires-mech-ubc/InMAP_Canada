import xarray as xr
import geopandas as gpd
from rasterstats import zonal_stats
import rioxarray as rio
import pandas as pd
import numpy as np
import sys
from joblib import Parallel, delayed


#salloc --time=3:0:0 --mem-per-cpu=10G --ntasks=10 --account=def-rscholes
#pip install --no-index --upgrade pip 
#pip install --no-index xarray geopandas

from joblib import Parallel, delayed
import geopandas as gpd
import pandas as pd
from rasterstats import zonal_stats

def compute_stats_for_geom(geom, image, affine, stats, crs):
    # Wrap geometry in a GeoDataFrame with correct CRS
    single_shape = gpd.GeoDataFrame(geometry=[geom], crs=crs)
    result = zonal_stats(single_shape, image, affine=affine, stats=stats)[0]
    #result = zonal_stats(single_shape, image, stats=stats)[0]#No affine since same pcr
    return result

def ZonalStatsParallel(shape, raster, rastcols, stats=['sum'], n_jobs=4):
    # shape - geopandas shapefile object
    # image- xarray dataset
    # affine - raster affine (xrdataset.rio.transform() to get)
    # stats - stats as list, i.e. ['min', 'mean', 'max'] ; ['sum']
    # the result is final_gdf as GeoDataFrame
    affine = raster.rio.transform()
    geometries = shape.geometry
    crs = shape.crs  # Automatically pull CRS from shapefile

    df_concat = shape.copy()

    for rc in rastcols:
        image = raster[rc][0, :, :].values
        #Handle nodata
        nodata_value = raster[rc].rio.nodata
        if nodata_value is not None:
            image = np.where(image == nodata_value, np.nan, image)


        # Parallelize over geometries
        results = Parallel(n_jobs=n_jobs)(
            delayed(compute_stats_for_geom)(geom, image, affine, stats, crs)
            for geom in geometries
        )

        # Convert results to DataFrame
        zonal_df = pd.DataFrame(results)
        zonal_df.columns = [f"{rc}_{s}" for s in stats]

        # Concatenate with original shape
        df_concat = pd.concat([df_concat.reset_index(drop=True), zonal_df], axis=1)
        print(f"Done {rc}")

    return gpd.GeoDataFrame(df_concat, geometry='geometry', crs=crs)



domajor=True
doarea=True
scenarios={1:['BASEGM_2015_E010','BASEGM_2015_017'],2:['BAU_2020_E108','BAU_2020_105'],3:['BAU_2025_E107','BAU_2025_106'],4:['BAU_2030_E107','BAU_2030_106'],
5:['BAU_2035_E107','BAU_2035_106'],6:['BPT_2015_E002','BPT_2015_002'],7:['BPT_2015_E006','BPT_2015_006'],
8:['BPT_2015_E015','BPT_2015_015'],9:['COVID_2020_E003','COVID_2020_003'],10:['BPT_2015_E016','BPT_2015_016']}

#shp = gpd.read_file('/home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/BASEGM_2015_E010_Sumlcc.shp')
ncpath = '/home/tfmrodge/scratch/GEMMACH_data/data/emissions_nc/20260223/'
fpath = '/home/tfmrodge/scratch/GEMMACH_data/data/emissions_nc/20260223/'
sID = int(sys.argv[1])
print(scenarios[sID][0])
scenario = scenarios[sID][0]

#shp = gpd.read_file('/home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/BASEGM_2015_E010_Sumlcc.shp')
shp = '/home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/shp_template/InMAPGEMMACH_area_template_lcc.shp'
#nc_area = '/home/tfmrodge/scratch/GEMMACH_data/data/Inmap_outputs/Emissions/BASEGM_2015_E010_area_lcc.nc'
#major = '/home/tfmrodge/scratch/GEMMACH_data/data/BASEGM_2015_E010/BASEGM_2015_E010_major_lcc.nc'
#nc_major = f"{ncpath}{scenario}_major_lcc.nc"

outshp = gpd.read_file(shp)
if doarea:
    nc_area = f"{ncpath}{scenario}_area_lcc.nc"
    in_area = xr.open_dataset(nc_area,engine='rasterio')
    rastcols=['NH3','NOx','PM25','SOx','VOC']
    outarea = ZonalStatsParallel(outshp, in_area,rastcols, ['sum'],n_jobs=64)
    outarea = outarea.to_crs('EPSG:4326')
    outarea.to_file(f"{fpath}202600223_{scenario}_area_EPSG4326.shp")
if domajor:
    nc_major = f"{ncpath}{scenario}_major_lcc.nc"
    in_major = xr.open_dataset(nc_major,engine='rasterio')
    rastcols=['NH3','NOx','PM25','SOx','VOC']
    outmajor = ZonalStatsParallel(outshp, in_major,rastcols, ['sum'],n_jobs=64)
    outmajor = outmajor.to_crs('EPSG:4326')
    outmajor.to_file(f"{fpath}202600223_{scenario}_majorgrid_EPSG4326.shp")
#outarea.set_crs(crs=outshp.crs,inplace=True) #Set CRS to make the same

# outmajor = ZonalStats(shp, major, 'sum')
# outmajor.to_file('/home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/BASEGM_2015_E010_major_lcc.shp')
#gpd.read_file('/home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/BAU_2020_E108_majorpts.shp')