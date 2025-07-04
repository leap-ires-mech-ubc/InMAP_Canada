import xarray as xr
import geopandas as gpd
from rasterstats import zonal_stats
import rioxarray as rio
import pandas as pd
import numpy as np

#module load python/3.10.2 scipy-stack
#source home/inmapenv/bin/activate
#salloc --time=1:0:0 --mem-per-cpu=3G --ntasks=4 --account=def-rscholes
#pip install --no-index --upgrade pip 
#pip install --no-index xarray geopandas

def ZonalStats(shape, raster,rastcols, stats=['sum']):
    # shape - geopandas shpaefile object
    # raster - xarray dataset
    # affine - raster affine (xrdataset.rio.transform() to get)
    # stats - stats as list, i.e. ['min', 'mean', 'max'] ; ['sum']
    # the result is final_gdf as GeoDataFrame
    
    #df = pd.DataFrame(zonalSt)
    affine = raster.rio.transform() 
    df_concat = shape
    for rc in rastcols:
        image = raster[rc][0,:,:]
        zonalst = pd.DataFrame(pd.DataFrame(zonal_stats(shape, image.values,affine=affine, stats = stats)).values,columns=[rc])
        df_concat = pd.concat([zonalst, df_concat], axis=1)
        print('done '+rc)
    final_gdf = gpd.GeoDataFrame(df_concat, geometry=df_concat.geometry)
    return final_gdf
#shp = gpd.read_file('/home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/BASEGM_2015_E010_Sumlcc.shp')
shp = '/home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/shp_template/BASEGM_polygon_template_wgs84.shp'
nc_area = '/home/tfmrodge/scratch/GEMMACH_data/data/Inmap_outputs/Emissions/BASEGM_2015_E010_area_lcc.nc'
#major = '/home/tfmrodge/scratch/GEMMACH_data/data/BASEGM_2015_E010/BASEGM_2015_E010_major_lcc.nc'

outshp = gpd.read_file(shp)
in_area = xr.open_dataset(nc_area,engine='rasterio')
#raster = in_area.sel(x=[500,501,502,503],y=[500,501,502,503],method='nearest')
#in_area = in_area[rastcols]
#affine = in_area.rio.transform()
# for rc in rastcols:
#     img = in_area[rc][0,:,:]
#     zonal_stats(outshp, img.values,affine=affine, stats ='sum')
#outmajor = xr.open_dataset(area)
 #zonal_stats(outshp, in_area, stats ='sum')
rastcols=['NH3','NOx','PM25','SOx','VOC']
outarea = ZonalStats(outshp, in_area,rastcols, ['sum'])
outarea.set_crs(crs=outshp.crs,inplace=True) #Set CRS to make the same
outarea.to_file('/home/tfmrodge/scratch/GEMMACH_data/data/Inmap_outputs/Emissions/Emissions_shp/20250624_BASEGM_2015_E010_area.shp')
# outmajor = ZonalStats(shp, major, 'sum')
# outmajor.to_file('/home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/BASEGM_2015_E010_major_lcc.shp')