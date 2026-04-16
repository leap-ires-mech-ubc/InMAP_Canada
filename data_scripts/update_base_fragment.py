def update_base(outfile,gemmfile=None,crs=None):
    '''Function to map GEMMACH output data onto INMAP-GEMMACH outputs. 
    outfile (str): Filepath to the output shapefile you want to load
    gemmfile (str): Gemmach output data you want to map onto a shapefile
    crs (None or str): Crs you want to convert to. If None, everything will be in the shapefile CRS
    '''
    #First, convert to the given CRS
    if crs != None:
        outgdf = gpd.read_file(outfile).to_crs(crs)
    else:
        outgdf = gpd.read_file(outfile)
        crs = outgdf.crs
    #Next, load the raster data

    if calcdelta:
        for vals in [['BasePM25','TotalPM25'],['BasePNO3','PNO3'],['BasePNH4',
                        'PNH4'],['BasePSO4','PSO4'],['BaseSOA','SOA'],['BasePrimPM25','PrimPM25']]:
            delname = 'delta_'+vals[0][4:]
            if vals[0] == 'BasePrimPM25':
                #Need to calculate out primary PM2.5 in the base scenario
                inmap_outs.loc[:,vals[0]] = (inmap_outs.BasePM25) - (inmap_outs.loc[:,['BasePNO3','BasePNH4','BasePSO4','BaseSOA']].sum(axis=1))
            inmap_outs.loc[:,delname] = inmap_outs.loc[:,vals[1]] - inmap_outs.loc[:,vals[0]] 
    


    #Next, load the basefile if it is present. If basefile is loaded, all values are differences
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
    #Clip to clipped geometry, if requested
    if type(clipped) != bool:
        outgdf = gpd.clip(outgdf,clipped.geometry)
    return outgdf