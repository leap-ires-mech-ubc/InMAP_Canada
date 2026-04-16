# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 14:54:03 2023

@author: trodge01.
#2025-10-15 Modified code with help of CoPilot. Fixed bug where the regression was backwards and MV was calculated incorrectly.
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
import math

def load_inmap(outfile,basefile=None,crs='ESRI:102002',clipped=False):
    ''' Load inmap data or read in as a change from the baseline (units are still concentration) 
    outfile (str): Filepath to the output shapefile you want to load
    basefile (str or gdf): Baseline data you want to subtract to get the change in emissions for a scenario. Loading this way means all values are deltas
        Can also have basefile already loaded in as a gdf to reduce time
    crs (None or str): Crs you want to convert to. If None, will not change CRS
    clipped (FALSE or gdf): Clip to the given gdf geometry (e.g. a shapefile of Canada, ), note this is SLOW
    '''
    #First, convert to the given CRS
    if crs != None:
        outgdf = gpd.read_file(outfile).to_crs(crs)
    else:
        outgdf = gpd.read_file(outfile)
    #Next, load the basefile if it is present. If basefile is loaded, all values are differences
    if type(basefile) != type(None):
        if isinstance(basefile, str): #For a string, load the file and convert to same CRS
            if crs != None:
                basegdf = gpd.read_file(basefile).to_crs(crs)
            else:
                basegdf = gpd.read_file(basefile)
        else:
            basegdf = basefile
        for i in ['TotalPM25','PNO3','PNH4','PSO4','PrimPM25','SOA']:
            #Difference between scenarios = scenario - baseline, +'ve means increased concentration in scenario
            outgdf.loc[:,i] =  outgdf.loc[:,i] - basegdf.loc[:,i] 
    #Clip to clipped geometry, if requested
    if type(clipped) != bool:
        outgdf = gpd.clip(outgdf,clipped.geometry)
    return outgdf

def summstats(df,pairs,stats,geoareas,popwt=None,geoname='PRENAME',popcol='TotalPop'):
    '''
    Calculate summary stats for InMAP outputs vs a reference. This was built for
    InMAP-Canada, and has specific inputs for that.
    df (pandas dataframe): Data frame with at least 2 columns containing the reference vs test "pairs" to be evaluated, 
    and a column that matches geoname to do regional analysis. 
    pairs (list): List of 2-item lists containing the names of the reference and the test values, in that order
    stats (list): List of the stats you want to calculate. Input to the calcstat() function, so must comply
    geoareas (list): List of the regions you want to calculate, must be entries in the geoname column. Using "all"
    or "Canada" in this will take the stats across all values, with no filtering
    popwt (bool): switch to do population weighting or not, with the column name from popcol
    geoname (string): column for the geoareas filtering
    '''
    statdfs = {}
    indname = ''
    for geoarea in geoareas:
        statdf = pd.DataFrame(columns=['Location']+stats[:-1])
        if geoarea in ['Canada', 'All']:
            pltdata = df
        else:
            pltdata = df.loc[df.PRENAME==geoarea]
        if popwt is not None: #Weight if you want to by popcol
            popwt = pltdata.loc[:,popcol]
        for pair in pairs: #Run through the pairs
            #
            indname = pair[0][4:]
            statdf.loc[indname,'Location'] = geoarea
            for stat in stats:
                try: #Regression is set up slightly differently
                    if stat != 'Regression':
                        statdf.loc[pair[0][4:],stat] = calcstat(stat,pltdata.loc[:,pair[0]],pltdata.loc[:,pair[1]],popwt)
                    else:
                        try:#regression can be a problem if it doesn't work - in that case, return nan
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

def calcstat(stat,ref,test,popwt=None):
    '''Wrapper function for stats. Currently allows:
        RMSE, MeanBias (MB), MeanError (ME), MeanFractionalBias (MFB), 
        MeanFractionalError (MFE), ModelRatio (MR), Ratio of Means (RoM), regression (reg), 
        number of observations (numobs)
        Note that regression gives an output with the slope and r²
        stat (str): Stats that will be calculated. Must match
        the strings in the if statements below
        ref (pd series): reference value (x)
        test (pd series): test value (y)
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
    elif (stat == 'RatioofMeans') | (stat == 'RoM'):
        stat =  RatioOfMeans(ref,test,popwt=popwt)
    elif (stat == 'Regression') | (stat == 'reg'):
        stat =  Regression(ref,test,popwt=popwt)
    elif (stat == 'numobs') | (stat == 'n'):
        stat =  len(test)
    return stat
    
#Define functions to calculate model-model statistics
# ---- Internal helpers (no interface changes) ----
def _mask_valid(ref, test, popwt=None, require_positive_weights=True):
    """
    Mask to finite ref/test (and weights if provided). Returns masked arrays.
    Prints and returns (None, None, None) if nothing valid remains.
    """
    try:
        ref = np.asarray(ref, dtype=float)
        test = np.asarray(test, dtype=float)
    except Exception as e:
        print(f"[mask] Could not convert inputs to float arrays: {e}")
        return None, None, None

    if ref.shape != test.shape:
        print(f"[mask] Shape mismatch: ref{ref.shape} vs test{test.shape}")
        return None, None, None
    valid = np.isfinite(ref) & np.isfinite(test)

    w = None
    if popwt is not None:
        try:
            w = np.asarray(popwt, dtype=float)
        except Exception as e:
            print(f"[mask] Could not convert weights to float array: {e}")
            return None, None, None
        if w.shape != ref.shape:
            print(f"[mask] Weight shape mismatch: weights{w.shape} vs data{ref.shape}")
            return None, None, None
        w_valid = np.isfinite(w)
        if require_positive_weights:
            w_valid &= (w > 0)
        valid &= w_valid

    if not np.any(valid):
        print("[mask] No valid data points after masking (check NaNs/infs/weights).")
        return None, None, None

    ref_m = ref[valid]
    test_m = test[valid]
    w_m = w[valid] if w is not None else None
    return ref_m, test_m, w_m

def _weighted_mean(x, w=None):
    if w is None:
        return np.nanmean(x)  # nan-safe mean
    sw = np.sum(w)
    if not np.isfinite(sw) or sw <= 0:
        return np.nan
    return np.sum(w * x) / sw

#Define functions to calculate model-model statistics
def RMSE(ref,test,popwt=None):
    #Calculate the RMSE for two series. These must be the same length. 
    #ref - reference value
    #test - test value, non-reference
    #popwt = None if area-weighted, series with population of each cell if pop weighted
    ref, test, w = _mask_valid(ref, test, popwt)
    if ref is None:
        print("[RMSE] Invalid inputs; returning np.nan.")
        return np.nan
    delta2 = (test - ref)**2
    try:
        return np.sqrt(_weighted_mean(delta2, w))
    except Exception as e:
        print(f"[RMSE] Error computing RMSE: {e}")
        return np.nan

def MeanBias(ref,test,popwt=None):
    #Calculate the mean bias (MB) for two series. These must be the same length. 
    #ref - reference value
    #test - test value, non-reference
    #popwt = None if area-weighted, series with population of each cell if pop weighted
    ref, test, w = _mask_valid(ref, test, popwt)
    if ref is None:
        print("[MeanBias] Invalid inputs; returning np.nan.")
        return np.nan
    delta = (test - ref)
    try:
        return _weighted_mean(delta, w)
    except Exception as e:
        print(f"[MeanBias] Error computing MB: {e}")
        return np.nan
def MeanError(ref,test,popwt=None):
    #Calculate the mean error (ME) for two series. These must be the same length. 
    #ref - reference value
    #test - test value, non-reference
    #popwt = None if area-weighted, series with population of each cell if pop weighted
    ref, test, w = _mask_valid(ref, test, popwt)
    if ref is None:
        print("[MeanError] Invalid inputs; returning np.nan.")
        return np.nan
    delta = np.abs(test - ref)
    try:
        return _weighted_mean(delta, w)
    except Exception as e:
        print(f"[MeanError] Error computing ME: {e}")
        return np.nan

def MeanFractionalBias(ref,test,popwt=None,threshold=1e-9):
    #Calculate the mean Fractional Bias (MFB) for two series. These must be the same length. 
    #ref - reference value
    #test - test value, non-reference
    #popwt = None if area-weighted, series with population of each cell if pop weighted
    #threshold for exclusion - very small values inflate this stat
    ref, test, w = _mask_valid(ref, test, popwt)
    if ref is None:
        print("[MeanFractionalBias] Invalid inputs; returning np.nan.")
        return np.nan
    denom = test + ref
    threshold = threshold * np.nanmedian(np.abs(denom))
    mask = np.isfinite(denom) &  (np.abs(denom) > threshold)
    excl = np.size(denom) - np.count_nonzero(mask)
    if excl > 0:
        print(f"[MeanFractionalBias] Excluded {excl} pairs with |test+ref| <= {threshold}.")
    if not np.any(mask):
        print("[MeanFractionalBias] All denominators are zero/invalid; returning np.nan.")
        return np.nan
    try:
        mfb_terms = 2.0 * (test[mask] - ref[mask]) / denom[mask]
        w_mask = w[mask] if w is not None else None
        return _weighted_mean(mfb_terms, w_mask)
    except Exception as e:
        print(f"[MeanFractionalBias] Error computing MFB: {e}")
        return np.nan

def MeanFractionalError(ref,test,popwt=None,threshold=1e-9):
    #Calculate the mean Fractional Error (MFE) for two series. These must be the same length. 
    #ref - reference value
    #test - test value, non-reference
    #popwt = None if area-weighted, series with population of each cell if pop weighted
    #threshold for exclusion - very small values inflate this stat
    ref, test, w = _mask_valid(ref, test, popwt)
    if ref is None:
        print("[MeanFractionalError] Invalid inputs; returning np.nan.")
        return np.nan
    denom = test + ref
    threshold = threshold * np.nanmedian(np.abs(denom))
    mask = np.isfinite(denom) &  (np.abs(denom) > threshold)
    excl = np.size(denom) - np.count_nonzero(mask)
    if excl > 0:
        print(f"[MeanFractionalError] Excluded {excl} pairs with |test+ref| <= {threshold}.")
    if not np.any(mask):
        print("[MeanFractionalError] All denominators are zero/invalid; returning np.nan.")
        return np.nan
    try:
        mfe_terms = 2.0 * np.abs(test[mask] - ref[mask]) / denom[mask]
        w_mask = w[mask] if w is not None else None
        return _weighted_mean(mfe_terms, w_mask)
    except Exception as e:
        print(f"[MeanFractionalError] Error computing MFE: {e}")
        return np.nan

def ModelRatio(ref, test, popwt=None):
    """
    Calculate the model ratio (MR):
        MR = mean(test/ref)

    If popwt is provided, computes weighted mean of ratios.

    ref  : reference values (x)
    test : test values (y)
    popwt: population or other weights (None for unweighted)
    """
    # Mask invalid values first
    ref, test, w = _mask_valid(ref, test, popwt)
    if ref is None:
        print("[ModelRatio] Invalid inputs; returning np.nan.")
        return np.nan

    # Avoid divide-by-zero problems
    valid = (ref != 0) & np.isfinite(ref) & np.isfinite(test)
    if not np.any(valid):
        print("[ModelRatio] No valid ref>0 values; returning np.nan.")
        return np.nan

    ref_v = ref[valid]
    test_v = test[valid]
    ratios = test_v / ref_v

    # Unweighted mean of ratios
    if w is None:
        return np.nanmean(ratios)

    # Weighted mean of ratios
    w_v = w[valid]
    w_sum = np.sum(w_v)
    if not np.isfinite(w_sum) or w_sum <= 0:
        print("[ModelRatio] Sum of weights is invalid; returning np.nan.")
        return np.nan
    return np.sum(w_v * ratios) / w_sum

def Regression(ref, test, popwt=None, origin=True):
    """
    Regression between test (y) and ref (x).

    Parameters
    ----------
    ref : array-like
        Reference values (x-axis)
    test : array-like
        Test values (y-axis)
    popwt : array-like or None
        Optional population weights
    origin : bool
        If True: regression is forced through origin (y = m*x)
        If False: standard regression with intercept (y = a + b*x)

    Returns
    -------
    slope, r2 : float, float
        Regression slope and coefficient of determination.
        If origin=True, R² is computed from the constrained model.
    """
    # Mask invalid data (NaN, inf, nonpositive weights)
    ref, test, w = _mask_valid(ref, test, popwt)
    if ref is None:
        print("[Regression] Invalid inputs; returning NaN.")
        return [np.nan, np.nan]
    # Reshape exogenous variable for statsmodels
    X = ref.reshape(-1, 1)
    y = test
    try:
        if origin:
            # Forced-through-origin regression: y = m*x
            if w is None:
                model = sm.OLS(y, X)
            else:
                model = sm.WLS(y, X, weights=w)

            res = model.fit()
            slope = float(res.params[0])
            # R² for origin regression is computed differently:
            # 1 - SSE / SST, but SST is centered on zero, not mean(y)
            y_pred = slope * ref
            ss_res = np.sum((y - y_pred) ** 2)
            ss_tot = np.sum(y ** 2)  # centered at zero
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
            return [slope, r2]
        else:
            # Standard regression: y = a + b*x
            X2 = sm.add_constant(X, has_constant="add")
            if w is None:
                model = sm.OLS(y, X2)
            else:
                model = sm.WLS(y, X2, weights=w)
            res = model.fit()
            intercept = float(res.params[0])
            slope = float(res.params[1])
            r2 = float(res.rsquared)
            return [slope, r2]
    except Exception as e:
        print(f"[Regression] Error during regression: {e}")
        return [np.nan, np.nan]

def MeanVal(ref,test,popwt=None):
    #Calculate the mean values (MB) for two series. These must be the same length. 
    #ref - reference value
    #test - test value, non-reference
    #popwt = None if area-weighted, series with population of each cell if pop weighted
    ref, test, w = _mask_valid(ref, test, popwt)
    if ref is None:
        print("[MeanVal] Invalid inputs; returning [np.nan, np.nan].")
        return [np.nan, np.nan]
    try:
        if w is None:
            MV_ref = np.nanmean(ref)
            MV_test = np.nanmean(test)
        else:
            MV_ref = _weighted_mean(ref, w)
            MV_test = _weighted_mean(test, w)
        return [MV_ref, MV_test]
    except Exception as e:
        print(f"[MeanVal] Error computing means: {e}")
        return [np.nan, np.nan]

def RatioOfMeans(ref, test, popwt=None):
    #Calculate the ratio of means (RoM) for two series. These must be the same length.
    #ref - reference value
    #test - test value, non-reference
    #popwt = None if area-weighted, series with population of each cell if pop weighted
    ref, test, w = _mask_valid(ref, test, popwt)
    if ref is None:
        print("[RatioOfMeans] Invalid inputs; returning np.nan.")
        return np.nan
    try:
        mean_ref = _weighted_mean(ref, w)
        mean_test = _weighted_mean(test, w)
        if mean_ref == 0 or not np.isfinite(mean_ref):
            print("[RatioOfMeans] Mean of reference is zero or invalid; returning np.nan.")
            return np.nan
        return mean_test / mean_ref
    except Exception as e:
        print(f"[RatioOfMeans] Error computing ratio of means: {e}")
        return np.nan

def plot_emissions(emissions,provinces,legend=True,lgdshk = 0.3,lnwdth = 0.05,alpha = 1.0,cmap=matplotlib.cm.YlOrRd,listvals=None,sjoinem=True,
                    figpath='/home/tfmrodge/scratch/GEMMACH_data/Figs/',scenario='test',diff=False,xylims=None,dopts=False):
    '''
    Function to plot emissions on a map. Can take either area or point emissions.
    emissions (gdf or list of 3 gdfs): Emissions to be plotted. If a list of 3 emissions is given, this will assume that
    emissions are in the order of area, major, combined (sum of area and major, or both on same plot if points)
    provinces (gdf): Outline to plot behind the emissions. 
    legend (bool): Plot with or without a legend
    lgdshk,lnwdth,alpha: Values that set the legendsize, linewidth, and alpha for plotting in matplotlib
    cmap (matplotlib cmap object): Base color map for matplotlib
    listvals (None or list of strings): Values you want to plot. None will plot NH3, NOx, PM25, SOX, VOC
    figpath (str): Path to save output figures
    scenario (str): Tag for naming files
    xylims (list): Map limits in same coordinates as emissions, provinces
    dopts (bool): Flag if the major emissions are pts (True) or area (False)
    '''
    if type(listvals)==type(None):
            listvals=['NH3','NOx','PM25','SOx','VOC']
    #figs ={}
    
    if len(emissions) ==3:#If emissions is a list of dataframes plot as area, major, combined
        triplot=True
        area=emissions[0]
        major=emissions[1]
        esum=emissions[2]
    else:
        triplot=False
        fig,axs = plt.subplots(math.ceil(len(listvals)/3),3,figsize = (12,6),dpi=300,sharex=True,sharey=True)
        axs=np.reshape(axs,3*math.ceil(len(listvals)/3))
        esum=emissions
    #Add provinces to esum as it always exists - everything is the same grid so can just do one
    if sjoinem:
        esum = sjoin(esum, provinces.loc[:,['PRENAME','geometry']], how='left',predicate='intersects')
        if (dopts & triplot):
            major = sjoin(major, provinces.loc[:,['PRENAME','geometry']], how='left',predicate='intersects')
        
    for ind,val in enumerate(listvals):
        if triplot:
            fig,axs = plt.subplots(1,3,figsize = (12,12),dpi=300,sharex=True,sharey=True)
            axs=np.reshape(axs,3)
            for ax in axs:
                provinces.geometry.boundary.plot(ax=ax, color=None, edgecolor='black',linewidth=0.1)
        else:
            provinces.geometry.boundary.plot(ax=axs[ind], color=None, edgecolor='black',linewidth=0.1)
        
        #Plot as area, major, combined for each pollutant
        if triplot:
            #Use provinces to set where canada is, use that for vlim and cmaps
            #Use the ranges to set the vlims and the cmap if it needs to be shifted
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
                area.plot(val,legend=legend,ax=axs[0],legend_kwds={'shrink':lgdshk,'label':"g/s"},linewidth=lnwdth,alpha=alpha,
                          cmap=cmap1,vmin=vlim1[0],vmax=vlim1[1])
                axs[0].set_title('Area')
            except ValueError:
                print('No emissions to plot')
            try:
                major.plot(val,legend=legend,ax=axs[1],legend_kwds={'shrink':lgdshk,'label':"g/s"},linewidth=lnwdth,alpha=alpha,
                           cmap=cmap2,vmin=vlim2[0],vmax=vlim2[1])
                axs[1].set_title('Major')
            except ValueError:
                print('No emissions to plot')
            try:
                esum.plot(val,legend=legend,ax=axs[2],legend_kwds={'shrink':lgdshk,'label':"g/s"},linewidth=lnwdth,alpha=alpha,
                          cmap=cmap3,vmin=vlim3[0],vmax=vlim3[1])
                axs[2].set_title('Sum')
                if dopts:
                    major.plot(val,legend=legend,ax=axs[1],legend_kwds={'shrink':lgdshk,'label':"g/s"},linewidth=lnwdth,alpha=alpha,
                               cmap=cmap2,vmin=vlim2[0],vmax=vlim2[1])
            except ValueError:
                print('No emissions to plot')
        else:
            #vlim = [0,max(max(esum.loc[~esum.PRENAME.isna(),val]),max(esum.loc[~esum.PRENAME.isna(),val]))]
            #vlim = [min(-1e-9,esum.loc[~esum.PRENAME.isna(),val]),max(-1e-9,esum.loc[~esum.PRENAME.isna(),val])]
            #To make sure 0 is centred properly, set min and max limits as 10% of max absolute value
            if sjoinem:
                minlim=max(1e-8,max(abs(esum.loc[~esum.PRENAME.isna(),val])))*0.1
                vlim = [min(-1*minlim,min(esum.loc[~esum.PRENAME.isna(),val])),max(minlim,max(esum.loc[~esum.PRENAME.isna(),val]))]
            else:
                minlim=max(1e-8,max(abs(esum.loc[:,val])))*0.1
                vlim = [min(-1*minlim,min(esum.loc[:,val])),max(minlim,max(esum.loc[:,val]))]
            cmap2 = shiftedColorMap(cmap, start=0, midpoint=1-vlim[1]/(vlim[1]+np.abs(vlim[0])), stop=1.0, name='shiftedcmap')
            try:
                if dopts:
                    esum.plot(val,legend=legend,ax=axs[ind],legend_kwds={'shrink':lgdshk,'label':"g/s"},linewidth=lnwdth,alpha=alpha,
                            cmap=cmap2,vmin=vlim[0],vmax=vlim[1],markersize=5.0)
                else:
                    esum.plot(val,legend=legend,ax=axs[ind],legend_kwds={'shrink':lgdshk,'label':"g/s"},linewidth=lnwdth,alpha=alpha,
                            cmap=cmap2,vmin=vlim[0],vmax=vlim[1])
                axs[ind].set_title(val)
            except ValueError:
                print('No emissions to plot')
        # if legend:
        #     for ax in axs:
        #         ax.get_legend().set_title("μg/m³")
        axs[0].set_xticks([]);
        axs[0].set_yticks([]);
        #Turn off the last subplot
        
        #Set limits
        if xylims is None:
            axs[0].set_xlim(-2579201.070414297, 3165870.);
            axs[0].set_ylim(76856.48815160134, 4270028.);
        else:
            axs[0].set_xlim(xylims[0])
            axs[0].set_ylim(xylims[1])
        if triplot:
            fig.savefig(figpath+scenario+'_EmissPlot_'+val+'.tif',format='tif')
        #figs[ind]=fig
    if not triplot:
        print(ind)
        if ind < len(axs):
            axs[len(axs)-1].set_axis_off()
        fig.savefig(figpath+scenario+'_EmissPlot_.tif',format='tif')
    return fig

def plot_pollutants(inmap_outs,provinces,legend=True,lgdshk = 0.3,lnwdth = 0.05,alpha = 1.0,cmap=matplotlib.cm.YlOrRd,listvals=None,
                    figpath='/home/tfmrodge/scratch/GEMMACH_data/Figs/',scenario='test',diff=False,xylims=None):
    '''
    Function to plot InMAP results on a map. 
    inmap_outs (gdf): File to be plotted. Should contain the values in listvals
    provinces (gdf): Outline to plot behind the concentrations. 
    legend (bool): Plot with or without a legend
    lgdshk,lnwdth,alpha: Values that set the legendsize, linewidth, and alpha for plotting in matplotlib
    cmap (matplotlib cmap object): Base color map for matplotlib. Note that currently the difference one is hardcoded as RdBu_r
    listvals (None or list of strings): Concentrations you want to plot. None will plot all
    figpath (str): Path to save output figures
    scenario (str): Tag for naming files
    diff (bool): Flag to say if you want to plot emissions as absolute values or as a difference
    xylims (list): Map limits in same coordinates as emissions, provinces
    '''
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
            # vlim1 = [min(min(inmap_outs.loc[:,vals[0]]),min(inmap_outs.loc[:,vals[1]]))
            #          ,max(max(inmap_outs.loc[:,vals[0]]),max(inmap_outs.loc[:,vals[1]]))]
            vlim1 = [ min(inmap_outs.loc[:, vals[0]].min(), inmap_outs.loc[:, vals[1]].min()),
                    max(inmap_outs.loc[:, vals[0]].max(), inmap_outs.loc[:, vals[1]].max())]
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
        for ax in np.atleast_1d(axs):
            ax.set_xticks([])
            ax.set_yticks([])
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