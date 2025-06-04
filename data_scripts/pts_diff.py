import xarray as xr
import geopandas as gpd
import rasterstats  as rs
import rioxarray as rio
import numpy as np
import pandas as pd
import sys

#module load python/3.10.2 scipy-stack mpi4py
#source inmapenv/bin/activate
#salloc --time=3:00:00 --mem-per-cpu=8G --nodes=1 --ntasks=10 --account=def-rscholes #def-agiang01 #def-rscholes
#pip install --no-index --upgrade pip 
#pip install --no-index xarray geopandas
scenarios={1:['BASEGM_2015_017','BASEGM_2015_017'],2:['BAU_2020_E108','BAU_2020_105'],3:['BAU_2025_E107','BAU_2025_106'],4:['BAU_2030_E107','BAU_2030_106'],
5:['BAU_2035_E107','BAU_2035_106'],6:['BPT_2015_E002','BPT_2015_002'],7:['BPT_2015_E006','BPT_2015_006'],
8:['BPT_2015_E015','BPT_2015_015'],9:['COVID_2020_E003','COVID_2020_003'],10:['BPT_2015_E016','BPT_2015_016']}

#shp = gpd.read_file('/home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/BASEGM_2015_E010_Sumlcc.shp')
fpath = '/home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/'
sID = int(sys.argv[1])
print(scenarios[sID])
scenario = scenarios[sID][0]
majorptsbase=gpd.read_file('/home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/BASEGM_2015_E010_majorpts_20240225.shp')
#Set index to the unique emitter identifier
majorptsbase=majorptsbase.set_index('ename')
#print(majorptsbase.columns)
majorptsscenario = gpd.read_file(fpath+scenario+'_majorpts_20240225.shp')
majorptsscenario=majorptsscenario.set_index('ename')
majorptsdiff = majorptsscenario.copy(deep=True)

#All diffs are scenario - base, then split.
#We want all points in scenario not in base. So, need to subtract for interect of indices.
#If an emitter in in scenario not base, we keep it as-is (no operation), assuming it is a new emissions
mask = majorptsdiff.index.intersection(majorptsbase.index)
majorptsdiff.loc[mask,['NH3', 'NOx', 'PM25', 'SOx', 'VOC']] = (majorptsdiff.loc[mask,['NH3', 'NOx', 'PM25', 'SOx', 'VOC']]
                                                          - majorptsbase.loc[mask,['NH3', 'NOx', 'PM25', 'SOx', 'VOC']])
#If an emitter is in "base" not "scenario", add as a negative to diff, assuming it is a lost emissions source
mask = majorptsbase.index.difference(majorptsdiff.index)
majorptsdiff = pd.concat([majorptsdiff,majorptsbase.loc[mask,:]])
majorptsdiff.loc[mask,['NH3', 'NOx', 'PM25', 'SOx', 'VOC']] = -majorptsdiff.loc[mask,['NH3', 'NOx', 'PM25', 'SOx', 'VOC']]
majorptsdiff.drop('index',axis=1).reset_index(inplace=True)
#Save pts diff after tidying
try:#Not sure why but doing in one line left NaNs
    majorptsdiff.loc[:,'Height']=majorptsdiff.loc[:,'Height_1']
    majorptsdiff.loc[:,'Velocity']=majorptsdiff.loc[:,'Velocity_1']
except KeyError:
    pass
keepcols = ['latitude', 'longitude','Height', 'Diam', 'Temp','Velocity','geometry']
evals = ['NH3', 'NOx', 'PM25', 'SOx', 'VOC']
majorptsdiff = majorptsdiff.loc[:,keepcols+evals]
majorptsdiff.to_file(fpath+scenario+'_majorptsdiff.shp')
#Split by mask
print(majorptsdiff.loc[:,['NH3', 'NOx', 'PM25', 'SOx', 'VOC']].sum(axis=0))

for val in evals:
        #Now do major
        negmajor = majorptsdiff.loc[:,[val,'geometry']].copy(deep=True)
        posmajor = majorptsdiff.loc[:,[val,'geometry']].copy(deep=True)
        negmajor.loc[:,val]=0.0
        posmajor.loc[:,val]=0.0
        negmask = majorptsdiff.loc[:,val]<0
        #Negative emissions make positive
        negmajor.loc[negmask,val] = -1*majorptsdiff.loc[negmask,val]
        posmajor.loc[~negmask,val] = majorptsdiff.loc[~negmask,val]
        #Save to file
        negmajor.to_file(fpath+scenario+'_'+val+'_negmajor.shp')
        posmajor.to_file(fpath+scenario+'_'+val+'_posmajor.shp')
