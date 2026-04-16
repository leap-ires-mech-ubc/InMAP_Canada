# ============================================
# Clean + Optimized Driver Script for InMAP-GEMMACH Evaluation
# ============================================

import matplotlib
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import geopandas as gpd
from geopandas.tools import sjoin_nearest
import inmap_funcs as ihf
import sys

# --------------------------------------------------------
# File paths and configuration
# --------------------------------------------------------
fpath = '/home/tfmrodge/scratch/GEMMACH_data/data/Inmap_outputs/Outputs/20260310/withGEMMACH/'
figpath = '/home/tfmrodge/scratch/GEMMACH_data/output_summaries/20260319/'

scenarios = {
    1:['BASEGM_2015_017','BASEGM_2015_017'],
    2:['BAU_2020_E108','BAU_2020_105'],
    3:['BAU_2025_E107','BAU_2025_106'],
    4:['BAU_2030_E107','BAU_2030_106'],
    5:['BAU_2035_E107','BAU_2035_106'],
    6:['BPT_2015_E002','BPT_2015_002'],
    7:['BPT_2015_E006','BPT_2015_006'],
    8:['BPT_2015_E015','BPT_2015_015'],
    9:['COVID_2020_E003','COVID_2020_003'],
    10:['BPT_2015_E016','BPT_2015_016']
}

pairs = [
    ['BASEPM25','TotalPM25'],
    ['BASEPNO3','PNO3'],
    ['BASEPNH4','PNH4'],
    ['BASEPSO4','PSO4'],
    ['BASESOA','SOA'],
    ['BASEPRIM25','PrimPM25']
]

stats = [
    'RMSE','MeanBias','MeanError',
    'MeanFractionalBias','MeanFractionalError',
    'ModelRatio','Regression'
]

remap = {
    'Canada':'Canada','Ontario':'Ontario',
    'Nova Scotia':'Atlantic','New Brunswick':'Atlantic',
    'Quebec':'Quebec','Prince Edward Island':'Atlantic',
    'Newfoundland and Labrador':'Atlantic','British Columbia':'British Columbia',
    'Saskatchewan':'Prairies','Alberta':'Prairies','Manitoba':'Prairies',
    'Northwest Territories':'North','Yukon':'North','Nunavut':'North'
}

geoareas = [
    'Canada','British Columbia','Quebec','Nunavut','Prince Edward Island',
    'Saskatchewan','Yukon','Manitoba','Ontario','New Brunswick',
    'Northwest Territories','Alberta','Newfoundland and Labrador','Nova Scotia'
]

georegions = ['Canada','North','British Columbia','Prairies','Ontario','Quebec','Atlantic']

prefix='20260310_InmapOuts_'
suffix='_static_withGEMMACH.shp'

# Allow CLI: "3" or "3,7,8"
arg = sys.argv[1]
if "," in arg:
    sIDs = [int(x) for x in arg.split(",")]
else:
    sIDs = [int(arg)]

calcdelta = True
plotouts = True
loadfile = True
calcstats = True
regions = False   # Toggle for province vs region aggregation

cmap = matplotlib.cm.RdBu_r

# --------------------------------------------------------
# Load province boundaries
# --------------------------------------------------------
provinces = gpd.read_file(
    '/home/tfmrodge/projects/def-rscholes/tfmrodge/InMAP_Canada/data_scripts/provinces_lcc2.gpkg'
)

# --------------------------------------------------------
# Loop over scenarios
# --------------------------------------------------------
for sID in sIDs:
    scenario = scenarios[sID][0]
    # --------------------
    # Load and spatially assign province names
    # --------------------
    inmap_outs = ihf.load_inmap(
        fpath + prefix + scenario + suffix,
        basefile=None, #Use if you want to calculate everything as difference from GEMMACH_BASE
        crs='ESRI:102002',
        clipped=False
    )

    # Assign nearest province (avoids duplicates)
    inmap_outs = gpd.sjoin_nearest(
        inmap_outs,
        provinces[['PRENAME','geometry']],
        how='left'
    )

    # --------------------
    # Optional region remapping
    # --------------------
    if regions:
        inmap_outs['Region'] = inmap_outs['PRENAME'].map(remap)
        geoname = 'Region'
        geoareas = georegions
    else:
        geoname = 'PRENAME'

    # --------------------
    # Delta calculation (INMAP - GEM-MACH) for given scenario
    # --------------------
    if calcdelta:
        for ref_name, test_name in pairs:
            ind = ref_name.replace("BASE", "")     # safer than slicing
            inmap_outs['delta_' + ind] = inmap_outs[test_name] - inmap_outs[ref_name]

    # --------------------
    # Plotting block
    # --------------------
    if plotouts:
        xylims_list = [((-2579201.0704,3165870.), (76856.48,4270028.))]
        figsuffix = ["_plot_CAN.tif"]
        if sID in [6]: #BC
            xylims_list.append(((-2200000, -1300000), (1190000, 1800000)))
            figsuffix.append("_plot_BC.tif")
        if sID in [7,8]: # ON-QC corridor
            xylims_list.append(((1088000, 1900000), (332000, 1250000)))
            figsuffix.append("_plot_ONQC.tif")
        for i, xylims in enumerate(xylims_list):
            fig = ihf.plot_pollutants(
                inmap_outs, provinces,
                legend=True, lgdshk=0.3, lnwdth=0.05, alpha=1.0,
                cmap=cmap, listvals=pairs, figpath=figpath,
                scenario=scenario, diff=True, xylims=xylims
            )
            figname = figpath + scenario + figsuffix[i]
            fig.savefig(figname, dpi=300)
            plt.close(fig)

    # --------------------
    # Stats
    # --------------------
    if calcstats:
        # Unweighted
        statdf = ihf.summstats(
            inmap_outs, pairs, stats, geoareas,
            popwt=None, geoname=geoname, popcol='TotalPop'
        )
        statdf.to_csv(figpath + scenario + '_summstats.csv')

        # Population-weighted
        statdf = ihf.summstats(
            inmap_outs, pairs, stats, geoareas,
            popwt=True, geoname=geoname, popcol='TotalPop'
        )
        statdf.to_csv(figpath + scenario + '_popsummstats.csv')