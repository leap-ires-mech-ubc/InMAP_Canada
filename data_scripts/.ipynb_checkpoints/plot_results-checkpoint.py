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
#Plot Concentrations vs GEMMACH and generate summary statistics
#module load python/3.10.2 scipy-stack/2023b mpi4py/3.1.3
#source $HOME/inmapenv/bin/activate
fpath = '/home/tfmrodge/scratch/GEMMACH_data/data/Inmap_outputs/Outputs/'#Beluga
figpath='/home/tfmrodge/scratch/GEMMACH_data/Figs/20240306/'
#List of scenarios with emissions and output naming schemes
scenarios={1:['BASEGM_2015_017','BASEGM_2015_017'],2:['BAU_2020_E108','BAU_2020_105'],3:['BAU_2025_E107','BAU_2025_106'],4:['BAU_2030_E107','BAU_2030_106'],
5:['BAU_2035_E107','BAU_2035_106'],6:['BPT_2015_E002','BPT_2015_002'],7:['BPT_2015_E006','BPT_2015_006'],
8:['BPT_2015_E015','BPT_2015_015'],9:['COVID_2020_E003','COVID_2020_003'],10:['BPT_2015_E016','BPT_2015_016']}
pairs=[['BasePM25','TotalPM25'],['BasePNO3','PNO3'],['BasePNH4',
        'PNH4'],['BasePSO4','PSO4'],['BaseSOA','SOA'],['BasePrimPM25','PrimPM25']]
stats = ['RMSE', 'MeanBias', 'MeanError', 'MeanFractionalBias', 
        'MeanFractionalError', 'ModelRatio', 'Regression']
remap = {'Canada':'Canada','Ontario':'Ontario','Nova Scotia':'Atlantic', 'New Brunswick':'Atlantic', 'Quebec':'Quebec',
       'Prince Edward Island':'Atlantic', 'Newfoundland and Labrador':'Atlantic',
       'British Columbia':'British Columbia', 'Saskatchewan':'Prairies', 'Alberta':'Prairies', 'Manitoba':'Prairies',
       'Northwest Territories':'North', 'Yukon':'North', 'Nunavut':'North'}
geoareas = ['Canada','British Columbia', 'Quebec', 'Nunavut', 'Prince Edward Island',
 'Saskatchewan', 'Yukon', 'Manitoba', 'Ontario', 'New Brunswick', 'Northwest Territories',
 'Alberta', 'Newfoundland and Labrador', 'Nova Scotia']
georegions = ['Canada','North','British Columbia','Prairies','Ontario','Quebec','Atlantic']
prefix='20240305_InmapOuts_'
suffix='_diff_combined.shp'
sIDs = [int(sys.argv[1])]  #[8] # [6,7,8]
calcdelta=True
plotouts=False #True #False #
loadfile=True #False #
calcstats=True
#Run stats on regions not provinces
regions=False #True
if regions:
    geoareas = georegions
cmap=matplotlib.cm.RdBu_r
if loadfile:
    provinces = gpd.read_file('/home/tfmrodge/projects/def-agiang01/tfmrodge/InMAP_Canada/data_scripts/provinces_lcc2.gpkg')#.to_crs(inmap_outs.crs)
for sID in sIDs:
    scenario = scenarios[sID][0]#Set the scenario name 
    if loadfile:
        inmap_outs = ihf.load_inmap(fpath+prefix+scenario+suffix,basefile=None,crs='ESRI:102002',clipped=False)#'GEMMACH_BASEGM_2015_017_50x50.shp'
        #Add provinces
        inmap_outs = sjoin(inmap_outs, provinces.loc[:,['PRENAME','geometry']], how='left',predicate='intersects')
        #inmap_outs.loc[:,['BasePM25','TotalPM25']].mean()
        #Reference (GEMMACH) - INMAP
        #'Base' data is the GEMMACH output, INMAP is the other one
        if calcdelta:
            for vals in [['BasePM25','TotalPM25'],['BasePNO3','PNO3'],['BasePNH4',
                              'PNH4'],['BasePSO4','PSO4'],['BaseSOA','SOA'],['BasePrimPM25','PrimPM25']]:
                delname = 'delta_'+vals[0][4:]
                if vals[0] == 'BasePrimPM25':
                    #pdb.set_trace()
                    inmap_outs.loc[:,vals[0]] = (inmap_outs.BasePM25) - (inmap_outs.loc[:,['BasePNO3','BasePNH4','BasePSO4','BaseSOA']].sum(axis=1))
                #if sID in [6,7,8,9]:#Put in for the ones that are base - scenario (default is scenario-base)
                #    inmap_outs.loc[:,vals[0]]=-1*inmap_outs.loc[:,vals[0]]
                inmap_outs.loc[:,delname] = inmap_outs.loc[:,vals[1]] - inmap_outs.loc[:,vals[0]] 
    if plotouts:
        if sID in [6]:
            #BC
            xylims=[(-2200000,-1300000),(1190000,1800000)]
            fig = ihf.plot_pollutants(inmap_outs,provinces,legend=True,lgdshk = 0.3,lnwdth = 0.05,alpha = 1.0,
                              cmap=cmap,listvals=pairs,figpath=figpath,scenario=scenario+'_BC_'+prefix[:8],diff=True,xylims=xylims)
        if sID in [7,8]:
            #ON-QC Corridor
            xylims=[(1088000,1900000),(332000,1250000)]
            fig = ihf.plot_pollutants(inmap_outs,provinces,legend=True,lgdshk = 0.3,lnwdth = 0.05,alpha = 1.0,
                              cmap=cmap,listvals=pairs,figpath=figpath,scenario=scenario+'_ONQC_'+prefix[:8],diff=True,xylims=xylims)
        #All of Canada - plot for all.
        xylims=[(-2579201.070414297,3165870.),(76856.48815160134,4270028.)]
        fig = ihf.plot_pollutants(inmap_outs,provinces,legend=True,lgdshk = 0.3,lnwdth = 0.05,alpha = 1.0,
                              cmap=cmap,listvals=pairs,figpath=figpath,scenario=scenario+'_'+prefix[:8],diff=True,xylims=xylims)
        fig.show()
    if calcstats:
        #popwt = None
        #popwt=True #None#inmap_outs.loc[:,'TotalPop']
        #statdf function - calculates summary stats for different pairs of values in different regions
        statdf = ihf.summstats(inmap_outs,pairs,stats,geoareas,popwt=None,geoname='PRENAME',popcol='TotalPop')
        statdf.to_csv(figpath+scenario+'_summstats.csv')
        statdf = ihf.summstats(inmap_outs,pairs,stats,geoareas,popwt=True,geoname='PRENAME',popcol='TotalPop')
        statdf.to_csv(figpath+scenario+'_popsummstats.csv')