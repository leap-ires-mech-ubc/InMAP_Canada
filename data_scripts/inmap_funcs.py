# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 14:54:03 2023

@author: trodge01
Functions to help analyze inmap data
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import statsmodels.api as sm
import geopandas as gpd
import pandas as pd
from geopandas.tools import sjoin

def load_inmap(outfile,basefile=None,crs='ESRI:102002',clipped=False):
    ''' Load inmap data or read in as a change from the baseline (units are still concentration) 
    outfile should be the filepath to the output shapefile you want to load
    basefile is the baseline data you want to subtract to get the change in emissions for a scenario. Loading this way means all values are deltas
    Can also have basefile already loaded in as a gdf to reduce time
    crs is the crs you want to convert to. If None, will not change CRS
    clipped will clip to the given gdf geometry, note this is SLOW
    '''
    if crs != None:
        outgdf = gpd.read_file(outfile).to_crs(crs)
    else:
        outgdf = gpd.read_file(outfile)
    #pdb.set_trace()
    if type(basefile) != type(None):
        if type(basefile) == 'str': #For a string, load the file and convert to same CRS
            if crs != None:
                basegdf = gpd.read_file(basefile).to_crs(crs)
            else:
                basegdf = gpd.read_file(basefile)
        else:
            basegdf = basefile
        for i in ['TotalPM25','PNO3','PNH4','PSO4','PrimPM25','SOA']:
            #Difference between scenarios = scenario - baseline, +'ve means increased concentration in scenario
            outgdf.loc[:,i] =  outgdf.loc[:,i] - basegdf.loc[:,i] 
    if type(clipped) != bool:
        outgdf = gpd.clip(outgdf,clipped.geometry)
    return outgdf

def summstats(df,pairs,stats,geoareas,popwt=None,geoname='PRENAME',popcol='TotalPop'):
    '''
    Calculate summary stats for InMAP outputs vs a reference. This was built for
    InMAP-Canada, and has specific inputs for that - not the most generalized
    df - data frame with at least 2 columns containing the reference vs test "pairs" to be evaluated, 
    and a column that matches geoname to do regional analysis. 
    pairs - List of 2-item lists containing the names of the reference and the test values, in that order
    stats - List of the stats you want to calculate. Input to the calcstat() function, so must comply
    geoareas - List of the regions you want to calculate, must be entries in the geoname column. Using "all"
    or "Canada" in this will take the stats across all values, with no filtering
    popwt - switch to do population weighting or not, with the column name from popcol
    geoname - string, column for the geoareas filtering
    '''
    statdf = pd.DataFrame(columns=['Location']+stats[:-1])
    statdfs = {}
    indname = ''
    for geoarea in geoareas:
        if (geoarea == 'Canada') | (geoarea == 'All'):
            pltdata = df
        else:
            pltdata = df.loc[df.PRENAME==geoarea]
        if popwt is not None: #Weight if you want to by popcol
            popwt = pltdata.loc[:,popcol]
        for pair in pairs: #Run through the pairs
            #Check for duplicate index names
            # if indname in statdf.index:
            #     indname = pair[0][4:]+'_'+pair[1]
            # else:
            #     indname = pair[0][4:]
            indname = pair[0][4:]
            statdf.loc[indname,'Location'] = geoarea
            for stat in stats:
                try:
                    if stat != 'Regression':
                        statdf.loc[pair[0][4:],stat] = calcstat(stat,pltdata.loc[:,pair[0]],pltdata.loc[:,pair[1]],popwt)
                    else:
                        try:
                            m,r2 = calcstat(stat,pltdata.loc[:,pair[0]],pltdata.loc[:,pair[1]],popwt)
                            statdf.loc[pair[0][4:],'Slope']=m
                            statdf.loc[pair[0][4:],'r²']=r2
                        except np.linalg.LinAlgError:
                            statdf.loc[pair[0][4:],'Slope']=np.nan
                            statdf.loc[pair[0][4:],'r²']=np.nan
                except ValueError:
                    print('Could not compute '+stat)
                    statdf.loc[pair[0][4:],stat] = np.nan
        #pdb.set_trace()
        statdfs[geoarea] = statdf.copy(deep=True)  
    #Put it all into one big dataframe
    statdf = pd.concat(statdfs.values()).reset_index()
    return statdf

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    ''' https://stackoverflow.com/questions/7404116/defining-the-midpoint-of-a-colormap-in-matplotlib
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower offset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax / (vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper offset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    #plt.register_cmap(cmap=newcmap)

    return newcmap

#plt.scatter(inmap_outs.loc[:,'BasePM25'],inmap_outs.loc[:,'TotalPM25']/1000)
def calcstat(stat,ref,test,popwt=None):
    '''Wrapper function for stats. Currently allows:
        RMSE, MeanBias (MB), MeanError (ME), MeanFractionalBias (MFB), 
        MeanFractionalError (MFE), ModelRatio (MR), regression (reg)
        Note that regression gives an output with the slope and r²
        '''
    if stat == 'RMSE':
        stat =  RMSE(ref,test,popwt=popwt)
    elif (stat == 'MeanBias') | (stat == 'MB'):
        stat =  MeanBias(ref,test,popwt=popwt)
    elif (stat == 'MeanError') | (stat == 'ME'):
        stat =  MeanError(ref,test,popwt=popwt)
    elif (stat == 'MeanFractionalBias') | (stat == 'MFB'):
        stat =  MeanFractionalBias(ref,test,popwt=popwt)
    elif (stat == 'MeanFractionalError') | (stat == 'MFE'):
        stat =  MeanFractionalError(ref,test,popwt=popwt)
    elif (stat == 'ModelRatio') | (stat == 'MR'):
        stat =  ModelRatio(ref,test,popwt=popwt)
    elif (stat == 'ModelRatio') | (stat == 'MR'):
        stat =  ModelRatio(ref,test,popwt=popwt)
    elif (stat == 'Regression') | (stat == 'reg'):
        stat =  Regression(ref,test,popwt=popwt)
    elif (stat == 'numobs') | (stat == 'n'):
        stat =  len(test)
    return stat
    
#Define functions to calculate model-model statistics
def RMSE(ref,test,popwt=None):
    #Calculate the RMSE for two series. These must be the same length. 
    #ref - reference value
    #test - test value, non-reference
    #popwt = None if area-weighted, series with population of each cell if pop weighted
    delta = test - ref
    if popwt is None:
        RMSE = np.sqrt((delta**2/len(delta)).sum())
    else:
        RMSE = np.sqrt(((popwt*delta**2)/(popwt.sum())).sum())
        # (inmap_outs.loc[:,'TotalPop']*inmap_outs.loc[:,'delta_TotPM25']**2)/(inmap_outs.loc[:,'TotalPop'].sum())).sum())
    return RMSE

def MeanBias(ref,test,popwt=None):
    #Calculate the mean bias (MB) for two series. These must be the same length. 
    #ref - reference value
    #test - test value, non-reference
    #popwt = None if area-weighted, series with population of each cell if pop weighted
    delta = test - ref
    if popwt is None:
        MB = (delta/len(delta)).sum()
    else:
        MB = ((popwt*delta)/(popwt.sum())).sum()
    return MB

def MeanError(ref,test,popwt=None):
    #Calculate the mean error (ME) for two series. These must be the same length. 
    #ref - reference value
    #test - test value, non-reference
    #popwt = None if area-weighted, series with population of each cell if pop weighted
    delta = np.abs(test - ref)
    if popwt is None:
        ME = (delta/len(delta)).sum()
    else:
        ME = ((popwt*delta)/(popwt.sum())).sum()
    return ME

def MeanFractionalBias(ref,test,popwt=None):
    #Calculate the mean Fractional Bias (MFB) for two series. These must be the same length. 
    #ref - reference value
    #test - test value, non-reference
    #popwt = None if area-weighted, series with population of each cell if pop weighted
    if popwt is None:
        MFB = ((2*(test - ref)/(test+ref))/len(test)).sum()
    else:#Not sure this makes sense
        MFB = ((popwt*(2*(test - ref)/(test+ref)))/(popwt.sum())).sum()
    return MFB

def MeanFractionalError(ref,test,popwt=None):
    #Calculate the mean Fractional Error (MFE) for two series. These must be the same length. 
    #ref - reference value
    #test - test value, non-reference
    #popwt = None if area-weighted, series with population of each cell if pop weighted
    if popwt is None:
        MFE = ((2*np.abs(test - ref)/(test+ref))/len(test)).sum()
    else:#Not sure this makes sense
        MFE = ((popwt*(2*np.abs(test - ref)/(test+ref)))/(popwt.sum())).sum()
    return MFE

def ModelRatio(ref,test,popwt=None):
    #Calculate the model ratio (MR) for two series. These must be the same length. 
    #ref - reference value
    #test - test value, non-reference
    #popwt = None if area-weighted, series with population of each cell if pop weighted
    if popwt is None:
        MR = ((test/ref)/len(test)).sum()
    else:#Not sure this makes sense
        MR = (((popwt*(test/ref)))/(popwt.sum())).sum()
    return MR

def Regression(ref,test,popwt=None):
    #Calculate the model ratio (MR) for two series. These must be the same length. 
    #ref - reference value
    #test - test value, non-reference
    #popwt = None if area-weighted, series with population of each cell if pop weighted
    if popwt is None:
        res = sm.OLS(ref,test,hasconst=False).fit()
        m,r2 = res.params[0],res.rsquared
    else:#For this, a weighted least squares
        res = sm.WLS(ref,test,weights=popwt).fit()
        m,r2 = res.params[0],res.rsquared
    Reg = [m,r2]
    return Reg

def MeanVal(ref,test,popwt=None):
    #Calculate the mean values (MB) for two series. These must be the same length. 
    #ref - reference value
    #test - test value, non-reference
    #popwt = None if area-weighted, series with population of each cell if pop weighted
    if popwt is None:
        MV_ref = (ref).sum()/len(ref)
        MV_test = (ref).sum()/len(ref)
    else:
        MV_ref = ((popwt*ref)/(popwt.sum())).sum()
    MVs=[MV_ref,MV_test]
    return MVs

def plot_emissions(emissions,provinces,legend=True,lgdshk = 0.3,lnwdth = 0.05,alpha = 1.0,cmap='YlOrRd',listvals=None,
                    figpath='/home/tfmrodge/scratch/GEMMACH_data/Figs/',scenario='test',diff=False,xylims=None,dopts=False):
    #pdb.set_trace()
    if type(listvals)==type(None):
            listvals=['NH3','NOx','PM25','SOx','VOC']
    #figs ={}
    for ind,val in enumerate(listvals):
        if len(emissions) ==3:#If emissions is a list of dataframes plot as area, major, combined
            fig,axs = plt.subplots(1,3,figsize = (12,12),dpi=300,sharex=True,sharey=True)
            axs=np.reshape(axs,3)
            area=emissions[0]
            major=emissions[1]
            esum=emissions[2]
        else:
            fig,axs = plt.subplots(2,3,figsize = (12,12),dpi=300,sharex=True,sharey=True)
            axs=np.reshape(axs,6)
            esum=emissions
        for ax in axs:
            provinces.geometry.boundary.plot(ax=ax, color=None, edgecolor='black',linewidth=0.1)
        #Add provinces to esum as it always exists - everything is the same grid so can just do one
        esum = sjoin(esum, provinces.loc[:,['PRENAME','geometry']], how='left',predicate='intersects')
        if dopts:
            major = sjoin(major, provinces.loc[:,['PRENAME','geometry']], how='left',predicate='intersects')
        #Use the ranges to set the vlims and the cmap if it needs to be shifted
        # if diff:
        #     vlim = [min(esum.loc[:,val]),max(esum.loc[:,val])]
        #     cmap = shiftedColorMap(cmap, start=0, midpoint=1-vlim[1]/(vlim[1]+np.abs(vlim[0])), stop=1.0, name='shiftedcmap')
        # else:
        #     vlim = [0,max(max(esum.loc[:,val]),max(esum.loc[:,val]))]
        #Plot as area, major, combined for each pollutant
        if len(emissions)==3:
            #Use provinces to set where canada is, use that for vlim and cmaps
            vlim1 = [min(area.loc[~esum.PRENAME.isna(),val]),max(area.loc[~esum.PRENAME.isna(),val])]
            if dopts:
                vlim2 = [min(major.loc[~major.PRENAME.isna(),val]),max(major.loc[~major.PRENAME.isna(),val])]
            else:
                vlim2 = [min(major.loc[~esum.PRENAME.isna(),val]),max(major.loc[~esum.PRENAME.isna(),val])]
            vlim3 = [min(esum.loc[~esum.PRENAME.isna(),val]),max(esum.loc[~esum.PRENAME.isna(),val])]
            cmap1 = shiftedColorMap(cmap, start=0, midpoint=1-vlim1[1]/(vlim1[1]+np.abs(vlim1[0])), stop=1.0, name='shiftedcmap')
            cmap2 = shiftedColorMap(cmap, start=0, midpoint=1-vlim2[1]/(vlim2[1]+np.abs(vlim2[0])), stop=1.0, name='shiftedcmap')
            cmap3 = shiftedColorMap(cmap, start=0, midpoint=1-vlim3[1]/(vlim3[1]+np.abs(vlim3[0])), stop=1.0, name='shiftedcmap')
            try:
                area.plot(val,legend=legend,ax=axs[0],legend_kwds={'shrink':lgdshk},linewidth=lnwdth,alpha=alpha,cmap=cmap1,vmin=vlim1[0],vmax=vlim1[1])
            except ValueError:
                print('No emissions to plot')
            try:
                major.plot(val,legend=legend,ax=axs[1],legend_kwds={'shrink':lgdshk},linewidth=lnwdth,alpha=alpha,cmap=cmap2,vmin=vlim2[0],vmax=vlim2[1])
            except ValueError:
                print('No emissions to plot')
            try:
                esum.plot(val,legend=legend,ax=axs[2],legend_kwds={'shrink':lgdshk},linewidth=lnwdth,alpha=alpha,cmap=cmap3,vmin=vlim3[0],vmax=vlim3[1])
                if dopts:
                    major.plot(val,legend=legend,ax=axs[1],legend_kwds={'shrink':lgdshk},linewidth=lnwdth,alpha=alpha,cmap=cmap2,vmin=vlim2[0],vmax=vlim2[1])
            except ValueError:
                print('No emissions to plot')
        else:
            vlim = [0,max(max(esum.loc[:,val]),max(esum.loc[:,val]))]
            cmap = shiftedColorMap(cmap, start=0, midpoint=1-vlim[1]/(vlim[1]+np.abs(vlim[0])), stop=1.0, name='shiftedcmap')
            try:
                esum.plot(val,legend=legend,ax=axs[ind],legend_kwds={'shrink':lgdshk},linewidth=lnwdth,alpha=alpha,cmap=cmap,vmin=vlim[0],vmax=vlim[1])
            except ValueError:
                print('No emissions to plot')
        axs[0].set_title(val[0])
        axs[0].set_xticks([]);
        axs[0].set_yticks([]);
        #Set limits
        if xylims is None:
            axs[0].set_xlim(-2579201.070414297, 3165870.);
            axs[0].set_ylim(76856.48815160134, 4270028.);
        else:
            axs[0].set_xlim(xylims[0])
            axs[0].set_ylim(xylims[1])
        fig.savefig(figpath+scenario+'_EmissPlot_'+val+'.tif',format='tif')
        #figs[ind]=fig
    return fig

def plot_pollutants(inmap_outs,provinces,legend=True,lgdshk = 0.3,lnwdth = 0.05,alpha = 1.0,cmap='YlOrRd',listvals=None,
                    figpath='/home/tfmrodge/scratch/GEMMACH_data/Figs/',scenario='test',diff=False,xylims=None):
    #pdb.set_trace()
    if type(listvals)==type(None):
            listvals=[['BasePM25','TotalPM25'],['BasePNO3','PNO3'],['BasePNH4',
                      'PNH4'],['BasePSO4','PSO4'],['BaseSOA','SOA'],['BasePrimPM25','PrimPM25']]
    #figs ={}
    for ind,vals in enumerate(listvals):
        fig,axs = plt.subplots(1,3,figsize = (12,12),dpi=300,sharex=True,sharey=True)
        for ax in axs:
            provinces.geometry.boundary.plot(ax=ax, color=None, edgecolor='black',linewidth=0.1)
        delname = 'delta_'+vals[0][4:]
        #Set the row - 6 sets of three
        try:
            axs.shape[1]
            ax = axs[ind]
        except IndexError:
            ax = axs
        #Use the ranges to set the vlims
        if diff:
            vlim1 = [min(min(inmap_outs.loc[:,vals[0]]),min(inmap_outs.loc[:,vals[1]]))
                     ,max(max(inmap_outs.loc[:,vals[0]]),max(inmap_outs.loc[:,vals[1]]))]
            plt_cmap = shiftedColorMap(cmap, start=0, midpoint=1-vlim1[1]/(vlim1[1]+np.abs(vlim1[0])), stop=1.0, name='shiftedcmap')
        else:
            vlim1 = [0,max(max(inmap_outs.loc[:,vals[0]]),max(inmap_outs.loc[:,vals[1]]))]
            plt_cmap =cmap
        vlim2 = [min(min(inmap_outs.loc[:,delname]),min(inmap_outs.loc[:,delname])),
                 max(max(inmap_outs.loc[:,delname]),max(inmap_outs.loc[:,delname]))]
        vlim = [(vlim1),(vlim2)]
        #Define shifted colormap
        cmap2 = shiftedColorMap(matplotlib.cm.RdBu_r, start=0, midpoint=1-vlim[1][1]/(vlim[1][1]+np.abs(vlim[1][0])), stop=1.0, name='shiftedcmap')
        #Plot as reference value, predicted value, and difference
        try:
            inmap_outs.plot(vals[0],legend=legend,ax=ax[0],legend_kwds={'shrink':lgdshk},linewidth=lnwdth,alpha=alpha,cmap=plt_cmap,vmin=vlim[0][0],vmax=vlim[0][1])
            inmap_outs.plot(vals[1],legend=legend,ax=ax[1],legend_kwds={'shrink':lgdshk},linewidth=lnwdth,alpha=alpha,cmap=plt_cmap,vmin=vlim[0][0],vmax=vlim[0][1])
            inmap_outs.plot(delname,legend=legend,ax=ax[2],legend_kwds={'shrink':lgdshk},linewidth=lnwdth,alpha=alpha,cmap=cmap2,vmin=vlim[1][0],vmax=vlim[1][1])
        except ValueError:
                print('No emissions to plot')
        ax[0].set_title(vals[0][4:])
        axs[0].set_xticks([])
        axs[0].set_yticks([])
        #Set limits
        if xylims is None:
            axs[0].set_xlim(-2579201.070414297, 3165870.)
            axs[0].set_ylim(76856.48815160134, 4270028.)
        else:
            axs[0].set_xlim(xylims[0])
            axs[0].set_ylim(xylims[1])
        #inmap_outs
        #state_outline = gpd.read_file('/Users/rivkahgf/Downloads/evaldata_v1.6.1/states.shp')
        #fig.savefig(figpath+'BaseCase_EvalPlot.pdf',format='pdf')
        fig.savefig(figpath+scenario+'_EvalPlot_'+vals[0][4:]+'.tif',format='tif')
        #figs[ind]=fig
    return fig