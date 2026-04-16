import xarray as xr
import geopandas as gpd
from rasterstats import zonal_stats
import rioxarray as rio
import pandas as pd
import numpy as np
import sys
from joblib import Parallel, delayed


#salloc --time=1:0:0 --mem-per-cpu=3G --ntasks=4 --account=def-rscholes
#pip install --no-index --upgrade pip 
#pip install --no-index xarray geopandas

def ZonalStats(shape, raster,rastcols, stats=['sum']):
    # shape - geopandas shpaefile object
    # raster - xarray dataset
    # affine - raster affine (xrdataset.rio.transform() to get)
    # stats - stats as list, i.e. ['min', 'mean', 'max'] ; ['sum']
    # the result is final_gdf as GeoDataFrame
    affine = raster.rio.transform() 
    df_concat = shape
    for rc in rastcols:
        image = raster[rc][0,:,:]
        zonalst = pd.DataFrame(pd.DataFrame(zonal_stats(shape, image.values,affine=affine, stats = stats)).values,columns=[rc])
        df_concat = pd.concat([zonalst, df_concat], axis=1)
        print('done '+rc)
    final_gdf = gpd.GeoDataFrame(df_concat, geometry=df_concat.geometry)
    return final_gdf

def compute_zonal_stat(shape, image, affine, rc, stats):
    # shape - geopandas shpaefile object
    # image- xarray dataset
    # affine - raster affine (xrdataset.rio.transform() to get)
    # stats - stats as list, i.e. ['min', 'mean', 'max'] ; ['sum']
    # the result is final_gdf as GeoDataFrame
    zonalst = pd.DataFrame(zonal_stats(shape, image.values, affine=affine, stats=stats))
    zonalst.columns = [rc]
    return zonalst

def ZonalStats(shape, raster, rastcols, stats=['sum'], n_jobs=4):
    affine = raster.rio.transform()
    df_concat = shape.copy()

    # Prepare inputs for parallel processing
    results = Parallel(n_jobs=n_jobs)(
        delayed(compute_zonal_stat)(
            shape,
            raster[rc][0, :, :],
            affine,
            rc,
            stats
        ) for rc in rastcols
    )

    # Concatenate results
    for zonalst in results:
        df_concat = pd.concat([zonalst, df_concat], axis=1)

    final_gdf = gpd.GeoDataFrame(df_concat, geometry=df_concat.geometry)




scenarios={1:['BASEGM_2015_E010','BASEGM_2015_017'],2:['BAU_2020_E108','BAU_2020_105'],3:['BAU_2025_E107','BAU_2025_106'],4:['BAU_2030_E107','BAU_2030_106'],
5:['BAU_2035_E107','BAU_2035_106'],6:['BPT_2015_E002','BPT_2015_002'],7:['BPT_2015_E006','BPT_2015_006'],
8:['BPT_2015_E015','BPT_2015_015'],9:['COVID_2020_E003','COVID_2020_003'],10:['BPT_2015_E016','BPT_2015_016']}

#shp = gpd.read_file('/home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/BASEGM_2015_E010_Sumlcc.shp')
ncpath = '/home/tfmrodge/scratch/GEMMACH_data/data/emissions_nc/'
fpath = '/home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/'
sID = int(sys.argv[1])
print(scenarios[sID][0])
scenario = scenarios[sID][0]

#shp = gpd.read_file('/home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/BASEGM_2015_E010_Sumlcc.shp')
shp = '/home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/shp_template/InMAPGEMMACH_area_template_lcc.shp'
#nc_area = '/home/tfmrodge/scratch/GEMMACH_data/data/Inmap_outputs/Emissions/BASEGM_2015_E010_area_lcc.nc'
#major = '/home/tfmrodge/scratch/GEMMACH_data/data/BASEGM_2015_E010/BASEGM_2015_E010_major_lcc.nc'
nc_area = f"{ncpath}{scenario}_area_lcc.nc"
outshp = gpd.read_file(shp)
in_area = xr.open_dataset(nc_area,engine='rasterio')
rastcols=['NH3','NOx','PM25','SOx','VOC']
outarea = ZonalStats(outshp, in_area,rastcols, ['sum'])
outarea = outarea.to_crs('EPSG:4326')
outarea.set_crs(crs=outshp.crs,inplace=True) #Set CRS to make the same
outarea.to_file(f"{fpath}20260223_{scenario}_area_EPSG4326.shp")
# outmajor = ZonalStats(shp, major, 'sum')
# outmajor.to_file('/home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/BASEGM_2015_E010_major_lcc.shp')