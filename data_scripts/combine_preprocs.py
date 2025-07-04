#This script will take chunked inmap preprocessor outputs and combine them 
#into a single, averaged file. It also removes NaN values and fixes the 
#Dz variable, which were causing problems running the model.
import xarray as xr 
import os
import numpy as np
xr.set_options(keep_attrs=True,display_max_rows=25)

#Path where the preprocessed files are located
inpath = '/home/tfmrodge/scratch/GEMMACH_data/data/Inmap_outputs/Preproc/monthly/'
#Path where you will save the output, combined preprocessor file
outpath = '/home/tfmrodge/scratch/GEMMACH_data/data/Inmap_outputs/Preproc/'
#Name that you want for the final output file
newname = 'inmapData_GEMMACH_BASEGM_2015_017_complete.nc'
#newname = '/home/tfmrodge/scratch/GEMMACH_data/data/Inmap_outputs/Preproc/inmapData_GEMMACH_BASEGM_2015_017_test.nc'
#Number of days in each file. 
daydict = {
    'Jan1-15': 15,'Jan16-31': 16,
    'Feb1-14': 14,'Feb15-28': 14,  
    'Mar01-15': 15,'Mar16-31': 16,
    'Apr1-15': 15,'Apr16-30': 15,
    'May1-15': 15,'May16-31': 15,
    'Jun01-15': 15,'Jun16-30': 15,
    'Jul01-15': 15,'Jul16-31': 16,
    'Aug01-15': 15,'Aug16-31': 16,
    'Sep1-15': 15,'Sep16-24': 9,'Sep25-30':6,
    'Oct01-15': 15,'Oct16-31': 16,
    'Nov01-15': 15,'Nov16-30': 15,
    'Dec01-15': 15,'Dec16-31': 16}
#daydict = {0:15,1:16}
daydenom = 365#31#365#334 #Number of days to average by
firstfile = True
for filename in os.listdir(inpath): #  fnames: # 
    print('loading '+filename)
    xr1 = xr.open_dataset(inpath+filename, decode_coords="all")
    #xrdict[filename] = xr1
    if firstfile == True:
        xrout = xr.zeros_like(xr1)
        firstfile = False
    #This will find the number of days in the file
    timechunk = filename[filename.rfind('_')+1:-3]
    ndays = daydict[timechunk]
    varbs = xr1.variables#['Dz','pS']
    for var in varbs:
        #Set nans properly
        if xr1[var].min()<-1e5:
            xr1[var] = xr1[var].where(xr1[var]!=xr1[var].min()) 
        if xr1[var].max()>1e5:
            xr1[var] = xr1[var].where(xr1[var]!=xr1[var].max()) 
        #Calculate the average in the overall dataset
        xrout[var] += xr1[var]*ndays/daydenom
        print(np.array([xrout[var].min(),xrout[var].mean(),xrout[var].max()]))
    print('Done '+timechunk)
#Save an intermediate file in case something goes wrong
xrout.to_netcdf(outpath+'inmapData_GEMMACH_BASEGM_2015_017.nc',format= "NETCDF3_64BIT")
#xrout = xr.open_dataset(inpath+filename, decode_coords="all")
print('Saved intermediate xrout to file')
#Set nans - output file has very large and small values, we are going to set to nan then fill.
print('Then we fill \'em up!')
xrout = xrout.ffill(dim='x')
xrout = xrout.ffill(dim='y')
xrout = xrout.ffill(dim='xStagger')
xrout = xrout.ffill(dim='yStagger')
xrout = xrout.bfill(dim='x')
xrout = xrout.bfill(dim='y')
xrout = xrout.bfill(dim='xStagger')
xrout = xrout.bfill(dim='yStagger')
print('done fillin')
#Set all dZ 0s to the average for the layer.
for zz in xrout.z:
    xrout['Dz'].loc[dict(z=zz)]= xrout.sel(z=zz)['Dz'].where(xrout.sel(z=zz)['Dz']!=0,xrout['Dz'].sel(z=zz).mean())
print('done dZ')
#Export
xrout.to_netcdf(outpath+newname,format='NETCDF3_64BIT')
print('done')
