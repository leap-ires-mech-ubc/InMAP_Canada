import xarray as xr
import geopandas as gpd
from rasterstats import zonal_stats
import rioxarray as rio

#module load python/3.10.2 scipy-stack
#source inmapenv/bin/activate
#salloc --time=1:0:0 --mem-per-cpu=3G --ntasks=4 --account=def-rscholes
#pip install --no-index --upgrade pip 
#pip install --no-index xarray geopandas

# def ZonalStats(shape, raster, stats):
#     # shape - shapefile path
#     # raster - raster path
#     # stats - stats as list, f.e. 'min mean max' ; 'min'
#     # the result is final_gdf as GeoDataFrame

#     shape_gdf = gpd.read_file(shape)
#     zonalSt = zonal_stats(shape, raster, stats = stats)
#     df = pd.DataFrame(zonalSt)
#     df_concat = pd.concat([df, shape_gdf], axis=1)
#     final_gdf = gpd.GeoDataFrame(df_concat, geometry=df_concat.geometry)
#     return final_gdf

directory = 
for file in os.listdir(directory):
    if sumdf == None:
        sumdf = pd.read_csv('file')
    newdf = pd.read_csv('file')
    if len(newdf)!=len(sumdf)
        print(file+' Not the same length as sum'
    sumdf.loc[:,'EA3':'EPC1'] = sumdf.loc[:,'EA3':'EPC1']+ newdf.loc[:,'EA3':'EPC1']
sumdf = gpd.GeoDataFrame(sumdf, geometry=gpd.GeoSeries.from_xy(maj12['longitude'], maj12['latitude']))
#maj11 = pd.read_csv('/home/tfmrodge/scratch/GEMMACH_data/BAU_2020_E108_major_01_1_pts.csv',chunksize=1000).get_chunk()
#maj11 = gpd.GeoDataFrame(maj11, geometry=gpd.GeoSeries.from_xy(maj11['longitude'], maj11['latitude']))
#maj12 = pd.read_csv('/home/tfmrodge/scratch/GEMMACH_data/BAU_2020_E108_major_01_2_pts.csv',chunksize=1000).get_chunk()
#maj12 = gpd.GeoDataFrame(maj12, geometry=gpd.GeoSeries.from_xy(maj12['longitude'], maj12['latitude']))
# maj11.loc[:,'EA3':'EPC1'].columns
#shp = gpd.read_file('/home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/BASEGM_2015_E010_Sumlcc.shp')
shp = '/home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/BASEGM_2015_E010_Sumlcc.shp'
area = '/home/tfmrodge/scratch/GEMMACH_data/data/BASEGM_2015_E010/BASEGM_2015_E010_area_lcc.nc'
major = '/home/tfmrodge/scratch/GEMMACH_data/data/BASEGM_2015_E010/BASEGM_2015_E010_major_lcc.nc'

outshp = gpd.read_file(shp)
outarea = xr.open_dataset(area)
outmajor = xr.open_dataset(area)


outarea = ZonalStats(shp, area, 'sum')
outarea.to_file('/home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/BASEGM_2015_E010_area_lcc.shp')
outmajor = ZonalStats(shp, major, 'sum')
outmajor.to_file('/home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/BASEGM_2015_E010_major_lcc.shp')