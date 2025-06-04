#import xarray as xr
import geopandas as gpd
import pandas as pd
#from rasterstats import zonal_stats
#import rioxarray as rio
import glob
import os
#import zipfile
import calendar
from datetime import datetime
import sys
#Problem files
# fpath='/home/tfmrodge/scratch/GEMMACH_data/data/BAU_2025_E107/BAU_2025_E107_major_07_2_pts.csv'
#fpath='/home/tfmrodge/scratch/GEMMACH_data/data/BPT_2015_E002/BPT_2015_E002_major_03_4_pts.csv'
#reader = pd.read_csv(fpath,sep=',',chunksize=1e6)
#specific_rows = [978853]
#df = pd.read_csv(fpath, skiprows = lambda x: x not in specific_rows)
#module load python/3.10.2 scipy-stack mpi4py
#source inmapenv/bin/activate
#
#
# scenarios=['BAU_2025_E107','BAU_2030_E107',
#  'BAU_2035_E107','BPT_2015_E002','BPT_2015_E006','BPT_2015_E015','COVID_2020_E003']
#scenarios = ['BAU_2025_E107']
scenarios = [sys.argv[1]]
print(scenarios)
fpath='/home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/'
cols = ['latitude', 'longitude', 'Height', 'Diam', 'Temp','Velocity']
for scenario in scenarios:
    sumdf=gpd.read_file(fpath+scenario+'_majorpts.shp')
    print(sumdf.columns)
    try:#Not sure why but doing in one line left NaNs
        sumdf.loc[:,'Height']=sumdf.loc[:,'Height_1']
        sumdf.loc[:,'Velocity']=sumdf.loc[:,'Velocity_1']
    except KeyError:
        pass
    sumdf.loc[:,'ename'] = sumdf[cols].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
    grpdf=sumdf.groupby('ename').first()
    emisscols = ['EA3', 'EA2', 'EETH', 'EC38', 'ETOL', 'EARO', 'EHCH', 'EALD',
       'ECRE', 'EISO', 'ECO', 'ENO', 'ENO2', 'ENH3', 'ESO2', 'EHON', 'EMEK',
       'ESO4', 'ECM1', 'EEC1', 'ENT1', 'EAM1', 'ESU1', 'EPC1', 'VOC', 'NOx',
       'SOx', 'PM25', 'NH3']#,'ename']
    grpdf.loc[:,emisscols] = sumdf.loc[:,emisscols+['ename']].groupby('ename').sum()
    sumdf = grpdf
    sumdf=sumdf.set_crs('EPSG:4326') #Still need to set crs
    sumdf = sumdf.reset_index() #Reset to numerical index
    print('saving')
    sumdf.to_file(fpath+scenario+'_majorpts_20240225.shp')
