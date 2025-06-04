#Plot emissions
import matplotlib
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import geopandas as gpd
from geopandas.tools import sjoin
#from geopandas.tools import sjoin
import inmap_funcs as ihf
import sys
#Plot emissions
#module load python/3.10.2 scipy-stack/2023b mpi4py/3.1.3
#source $HOME/inmapenv/bin/activate
#salloc --time=2:00:00 --mem-per-cpu=4G --nodes=1 --ntasks=10 --account=def-rscholes
#fpath = 'D:/Globus Connect Personal/Data/InMAP/Outputs/' #'D:/GlobusConnect/InMap/test/Inmap_outs/'
#Evaluate baseline against preprocessor /home/tfmrodge/scratch/GEMMACH_data/data/Inmap_outputs/Emissions/Emissions_shp/BPT_2015_E015_majorptsdiff_20240302.cpg
fpath = '/home/tfmrodge/scratch/GEMMACH_data/data/Inmap_outputs/Emissions/Emissions_shp/Diff/'#Beluga
ptspath = '/home/tfmrodge/scratch/GEMMACH_data/data/Inmap_outputs/Emissions/Emissions_shp/'#Beluga
figpath='/home/tfmrodge/scratch/GEMMACH_data/Figs/20240306/'
#List of scenarios with emissions and output naming schemes
scenarios={1:['BASEGM_2015_017','BASEGM_2015_017'],2:['BAU_2020_E108','BAU_2020_105'],3:['BAU_2025_E107','BAU_2025_106'],4:['BAU_2030_E107','BAU_2030_106'],
5:['BAU_2035_E107','BAU_2035_106'],6:['BPT_2015_E002','BPT_2015_002'],7:['BPT_2015_E006','BPT_2015_006'],
8:['BPT_2015_E015','BPT_2015_015'],9:['COVID_2020_E003','COVID_2020_003'],10:['BPT_2015_E016','BPT_2015_016']}
listvals=['NH3','NOx','PM25','SOx','VOC']
#prefix='20240122_InmapOuts_'
suffix='diff.shp' #BAU_2025_E107_majorpts_20240225.shx
sIDs = [int(sys.argv[1])] #[6] #
plotouts=True
loadfile=True #False
dopts=True #False
sumemissions=True
crs='ESRI:102002'
cmap=matplotlib.cm.RdBu_r
provinces = gpd.read_file('/home/tfmrodge/projects/def-agiang01/tfmrodge/InMAP_Canada/data_scripts/provinces_lcc2.gpkg')#.to_crs(inmap_outs.crs)
for sID in sIDs:
    scenario = scenarios[sID][0]#Set the scenario name  /home/tfmrodge/scratch/GEMMACH_data/data/Inmap_outputs/Emissions/Emissions_shp/Diff/BPT_2015_E006_areadiff.shp
    if loadfile:
        #inmap_outs = ihf.load_inmap(fpath+prefix+scenario+suffix,basefile=None,crs='ESRI:102002',clipped=False)#'GEMMACH_BASEGM_2015_017_50x50.shp'
        print('loading '+scenario)
        areaemissions = gpd.read_file(fpath+scenario+'_area'+suffix).to_crs(crs)
        #areaemissions = sjoin(areaemissions, provinces.loc[:,['PRENAME','geometry']], how='left',predicate='intersects')
        if dopts: #
            majoremissions=gpd.read_file(ptspath+scenario+'_majorptsdiff_20240302.shp').to_crs(crs)
        else:
            majoremissions=gpd.read_file(fpath+scenario+'_major'+suffix).to_crs(crs)
            #majoremissions = sjoin(majoremissions, provinces.loc[:,['PRENAME','geometry']], how='left',predicate='intersects')
        allemissions = areaemissions.copy(deep=True)
        if dopts:
            pass
        else:
            allemissions.loc[:,['NH3','NOx','PM25','SOx','VOC']]=allemissions.loc[:,['NH3','NOx','PM25','SOx','VOC']]+majoremissions.loc[:,['NH3','NOx','PM25','SOx','VOC']]
        print('loaded and summed for '+scenario)
        #Add provinces
        #allemissions = sjoin(allemissions, provinces.loc[:,['PRENAME','geometry']], how='left',predicate='intersects')
        #inmap_outs.loc[:,['BasePM25','TotalPM25']].mean()
        #Reference (GEMMACH) - INMAP
        #'Base' data is the GEMMACH output, INMAp is the other one 
    if plotouts:
        emisslist=[areaemissions, majoremissions,allemissions] #area, major, combined
        fig = ihf.plot_emissions(emisslist,provinces,legend=True,lgdshk = 0.3,lnwdth = 0.05,
                             alpha = 1.0,cmap=cmap,listvals=listvals,figpath=figpath,scenario=scenario+'_diff',diff=True,dopts=dopts)
        fig.show()
        print('plotted '+scenario)
        #for fig in figs:#Returns a dict of figs
        #    fig.show()
        
    #inmap_outs.loc[:,'delta_TotPM25'] = inmap_outs.loc[:,'TotalPM25'] - inmap_outs.loc[:,'BasePM25'] 
#Correct baseline values - some were incorrectly specified in preproc file
#inmap_outs.loc[:,].mean()

