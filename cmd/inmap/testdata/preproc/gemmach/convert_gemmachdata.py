# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 09:03:10 2023
Convert InMAP NETCDF4 to NETCDF3 

@author: trodge01
"""
#'''
import xarray as xr
import pdb
import os
inpath = 'D:/ECCC_TestData/Test2/BASEGM_2015_017/' #'D:/ECCC_TestData/Test2/RDPS_QC/'#'#'D:\ECCC_TestData\Test2\BASEGM_2015_017/'
#'D:\ECCC_TestData\Test2\GEOPHY_VF'
outpath = inpath+'processed/'#'E:/GitHub/InMAP_Canada/cmd/inmap/testdata/preproc/gemmach/rdpsqc/'
#codetime = time.time()
pdb.set_trace()
for filename in os.listdir(inpath):
    if int(filename[:8]) < 20190102:
        yr,mth,day,hr = filename[:4],filename[4:6],filename[6:8],filename[-5:-3]
        #newname = 'rdpsqctest'+'_'+str(yr)+'-'+str(mth)+'-'+str(day)+'_'+str(hr)+'_00_00.nc'
        newname = 'gemmach'+'_'+str(yr)+'-'+str(mth)+'-'+str(day)+'_'+str(hr)+'_00_00.nc'
        if os.path.exists(outpath+newname):
            continue
        #print(filename)
        xr1 = xr.open_dataset(inpath+filename, decode_coords="all")
        #xr1 = xr1.isel(rlat = slice(0,4),rlon = slice(0,4))
        #Set naming convention. Currently I am just matching the WRF one for convenience.
        
        #newname = 'gemmach'+'_'+str(yr)+'-'+str(mth)+'-'+str(day)+'_'+str(hr)+'_00_00.nc'
        xr1.to_netcdf(outpath+newname,format='NETCDF3_64BIT')
        xr1.close()

#'''
'''
for filename in os.listdir(inpath):
    if int(filename[:8]) < 20190102:
        #print(filename)
        xr1 = xr.open_dataset(inpath+filename, decode_coords="all")
        xr1 = xr1.coarsen(rlat=299,rlon=364,level1=5,level2=5,level3=5,level4=5,boundary="trim").mean()
        yr,mth,day,hr = filename[:4],filename[4:6],filename[6:8],filename[-5:-3]
        #Set naming convention. Currently I am just matching the WRF one for convenience.
        newname = 'gemmachtest'+'_'+str(yr)+'-'+str(mth)+'-'+str(day)+'_'+str(hr)+'_00_00.nc'
        xr1.to_netcdf(outpath+newname,format='NETCDF3_CLASSIC')
        
for filename in os.listdir(outpath):
    if filename[:3] != 'gem':
        yr,mth,day,hr = filename[:4],filename[4:6],filename[6:8],filename[-5:-3]
        newname = 'gemmachtest'+'_'+str(yr)+'-'+str(mth)+'-'+str(day)+'_'+str(hr)+'_00_00.nc'
        os.rename(outpath+filename,outpath+newname)
for filename in ['Gem_geophy.nc']:
    #if int(filename[:8]) < 20190102:
        #print(filename)
    xr1 = xr.open_dataset(inpath+filename, decode_coords="all")
    xr1 = xr1.coarsen(rlat=299,rlon=364,boundary="trim").mean()
    yr,mth,day,hr = filename[:4],filename[4:6],filename[6:8],filename[-5:-3]
    #Set naming convention. Currently I am just matching the WRF one for convenience.
    newname = 'gemmachtest'+'_'+str(yr)+'-'+str(mth)+'-'+str(day)+'_'+str(hr)+'_00_00.nc'
    xr1.to_netcdf(outpath+newname,format='NETCDF3_CLASSIC')

'''
'''

inpath = 'D:\ECCC_TestData\Test2\BASEGM_2015_017/'
outpath = 'E:/GitHub/InMAP_Canada/cmd/inmap/testdata/preproc/gemmach/'
#codetime = time.time()
pdb.set_trace()
for filename in os.listdir(inpath):
    if int(filename[:8]) <= 20190102:
        #print(filename)
        xr1 = xr.open_dataset(inpath+filename, decode_coords="all")
        xr1 = xr1.coarsen(rlat=299,rlon=364,level1=8,level2=8,level3=8,level4=8,boundary="trim").mean()
        yr,mth,day,hr = filename[:4],filename[4:6],filename[6:8],filename[-5:-3]
        #Set naming convention. Currently I am just matching the WRF one for convenience.
        newname = 'gemmachtest'+'_'+str(yr)+'-'+str(mth)+'-'+str(day)+'_'+str(hr)+'_00_00.nc'
        xr1.to_netcdf(outpath+newname,format='NETCDF3_CLASSIC')
        
for filename in os.listdir(outpath):
    if filename[:3] != 'gem':
        yr,mth,day,hr = filename[:4],filename[4:6],filename[6:8],filename[-5:-3]
        newname = 'gemmachtest'+'_'+str(yr)+'-'+str(mth)+'-'+str(day)+'_'+str(hr)+'_00_00.nc'
        os.rename(outpath+filename,outpath+newname)
        
        
        
#coorddict = {}
coorddf = pd.DataFrame(index=xr1.variables)
coorddf.loc[:,'coords'] = np.nan
#coords = 
pdb.set_trace()
for ind,var in enumerate(xr1.variables):    
    try:
        coorddf.loc[var,'coords'] = xr1[var].dims[1]
    except IndexError:
        coorddf.loc[var,'coords'] = None
    #coorddict[var] = xr1[var].dims
#coords = pd.concat(coorddict)
'''
'''

fpath = 'E:/GitHub/InMAP_Canada/cmd/inmap/testdata/preproc/gemmach/'
f1 = 'gemmach_2019-01-01_00_00_00.nc'
f2 = 'gemmach_2019-01-01_01_00_00.nc'
xr1 = xr.open_dataset(fpath+f1, decode_coords="all")
xr1 = xr1.coarsen(rlat=299,rlon=364,level1=8,level2=8,level3=8,level4=8,boundary="trim").mean()
xr1.to_netcdf(fpath+'gemmachtest_2019-01-01_00_00_00.nc')
xr2 = xr.open_dataset(fpath+f2, decode_coords="all")
xr2 = xr2.coarsen(rlat=299,rlon=364,level1=8,level2=8,level3=8,level4=8,boundary="trim").mean()
xr2.to_netcdf(fpath+'gemmachtest_2019-01-01_01_00_00.nc')
fpath = 'E:/GitHub/InMAP_Canada/cmd/inmap/testdata/preproc/'
f3 = 'wrfout_d01_2005-01-01_00_00_00.nc'
xr = xr.open_dataset(fpath+f3, decode_coords="all")
'''