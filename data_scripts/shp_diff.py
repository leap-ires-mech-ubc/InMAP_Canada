#This script will calculate the difference as the scenario minus the base case, 
#if you have the complete emissions for a scenario and a base case. It will also
# split them into negative and positive files for each emission category.

import geopandas as gpd
import sys

#module load python/3.10.2 scipy-stack mpi4py
#source inmapenv/bin/activate
#salloc --time=1:00:00 --mem-per-cpu=4G --nodes=1 --ntasks=10 --account=def-agiang01
#pip install --no-index --upgrade pip 
#pip install --no-index xarray geopandas
scenarios={1:['BASEGM_2015_017','BASEGM_2015_017'],2:['BAU_2020_E108','BAU_2020_105'],3:['BAU_2025_E107','BAU_2025_106'],4:['BAU_2030_E107','BAU_2030_106'],
5:['BAU_2035_E107','BAU_2035_106'],6:['BPT_2015_E002','BPT_2015_002'],7:['BPT_2015_E006','BPT_2015_006'],
8:['BPT_2015_E015','BPT_2015_015'],9:['COVID_2020_E003','COVID_2020_003'],10:['BPT_2015_E016','BPT_2015_016']}

#shp = gpd.read_file('/home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/BASEGM_2015_E010_Sumlcc.shp')
fpath = '/home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/'
sID = int(sys.argv[1])
print(scenarios[sID])
areabase=gpd.read_file('/home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/BASEGM_2015_E010_area.shp')
#majorptsbase=gpd.read_file('/home/tfmrodge/scratch/GEMMACH_data/data/Inmap_outputs/Emissions/Emissions_shp/BASEGM_2015_E010_majorpts.shp')
majorbase=gpd.read_file('/home/tfmrodge/scratch/GEMMACH_data/Emissions_shp/BASEGM_2015_E010_major.shp')
scenario = scenarios[sID][0]
#Area emissions
areascenario = gpd.read_file(fpath+scenario+'_area.shp')
areadiff = areascenario.copy(deep=True)
#All emissions are calculated as scenario - base then split by mask
areadiff.loc[:,['NH3', 'NOx', 'PM25', 'SOx', 'VOC']] = areascenario.loc[:,['NH3', 'NOx', 'PM25', 'SOx', 'VOC']] - areabase.loc[:,['NH3', 'NOx', 'PM25', 'SOx', 'VOC']]
#
# if sID in [2,3,4,5]:
#     areadiff.loc[:,['NH3', 'NOx', 'PM25', 'SOx', 'VOC']] = areascenario.loc[:,['NH3', 'NOx', 'PM25', 'SOx', 'VOC']] - areabase.loc[:,['NH3', 'NOx', 'PM25', 'SOx', 'VOC']]
# else:#Others are mostly reductions - base - scenario
#     areadiff.loc[:,['NH3', 'NOx', 'PM25', 'SOx', 'VOC']] = areabase.loc[:,['NH3', 'NOx', 'PM25', 'SOx', 'VOC']] - areascenario.loc[:,['NH3', 'NOx', 'PM25', 'SOx', 'VOC']]
areadiff.loc[:,'geometry'] = areascenario.loc[:,'geometry']
areadiff.to_file(fpath+scenario+'_areadiff.shp')
print('Saved area')
#Major Areas
majorscenario = gpd.read_file(fpath+scenario+'_major.shp')
majordiff = majorscenario.copy(deep=True)
if len(majorscenario) != len(majorbase):
    print('Lengths different base is'+len(majorbase)+' scenario is '+len(majorscenario))
#if sID in [2,3,4,5]:#All emissions are calculated as scenario - base then split by mask
majordiff.loc[:,['NH3', 'NOx', 'PM25', 'SOx', 'VOC']] = majorscenario.loc[:,['NH3', 'NOx', 'PM25', 'SOx', 'VOC']] - majorbase.loc[:,['NH3', 'NOx', 'PM25', 'SOx', 'VOC']]
# else:
#     majordiff.loc[:,['NH3', 'NOx', 'PM25', 'SOx', 'VOC']] = majorbase.loc[:,['NH3', 'NOx', 'PM25', 'SOx', 'VOC']] - majorscenario.loc[:,['NH3', 'NOx', 'PM25', 'SOx', 'VOC']]
majordiff.loc[:,'geometry'] = majorscenario.loc[:,'geometry']
majordiff.to_file(fpath+scenario+'_majordiff.shp')
print(majordiff.loc[:,['NH3', 'NOx', 'PM25', 'SOx', 'VOC']].sum(axis=0))
for val in ['NH3','NOx','PM25','SOx','VOC']:
        negarea = areadiff.loc[:,[val,'geometry']].copy(deep=True)
        posarea = areadiff.loc[:,[val,'geometry']].copy(deep=True)
        negarea.loc[:,val]=0.0
        posarea.loc[:,val]=0.0
        negmask = areadiff.loc[:,val]<0
        #Negative emissions make positive
        negarea.loc[negmask,val] = -1*areadiff.loc[negmask,val]
        posarea.loc[~negmask,val] = areadiff.loc[~negmask,val]
        #Save to file
        negarea.to_file(fpath+scenario+'_'+val+'_negarea.shp')
        posarea.to_file(fpath+scenario+'_'+val+'_posarea.shp')
        #Now do major
        negmajor = majordiff.loc[:,[val,'geometry']].copy(deep=True)
        posmajor = majordiff.loc[:,[val,'geometry']].copy(deep=True)
        negmajor.loc[:,val]=0.0
        posmajor.loc[:,val]=0.0
        negmask = majordiff.loc[:,val]<0
        #Negative emissions make positive
        negmajor.loc[negmask,val] = -1*majordiff.loc[negmask,val]
        posmajor.loc[~negmask,val] = majordiff.loc[~negmask,val]
        #Save to file
        negmajor.to_file(fpath+scenario+'_'+val+'_negmajor.shp')
        posmajor.to_file(fpath+scenario+'_'+val+'_posmajor.shp')

