import pandas as pd
import xarray as xr


def get_lat_lons_from_txt_file(loc, ptloc):
    ## read text file with points
    textpts_fname = '../data/{0}/latlon_{1}.txt'.format(loc, ptloc)
    df = pd.read_csv(textpts_fname, header=None, sep=' ', names=['latitude', 'longitude'], engine='python')
    df['longitude'] = df['longitude']*-1
    df = df
    lons = df['longitude'].values
    lats = df['latitude'].values
    
    return lons, lats

def calc_ivt_vars(fname, forecast, lons, lats):
    ## load the forecast data
    ds = xr.open_dataset(fname)
    ds = ds.assign_coords({"lon": (((ds.lon + 180) % 360) - 180)}) # Convert DataArray longitude coordinates from 0-359 to -180-179
    ds = ds.sel(lon=slice(-180., -1)) # keep only western hemisphere
    if forecast == 'ECMWF':
        ds = ds.rename({'forecast_time': 'forecast_hour'}) # need to rename this to match GEFS
        ds = ds.assign_coords({"forecast_hour": (ds.forecast_hour*3)}) # convert forecast time to forecast hour
    elif forecast == 'W-WRF':
        ds = ds.rename({'ensembles': 'ensemble'}) # need to rename this to match GEFS/ECMWF
    
    ## Calculate IVT for the maps
    ds1 = ds.mean(['forecast_hour', 'ensemble'])
    ds1 = ds1.where(ds.IVT >= 250)
        
    # subset ds to the select points
    x = xr.DataArray(lons, dims=['location'])
    y = xr.DataArray(lats, dims=['location'])
    ds = ds.sel(lon=x, lat=y, method='nearest')
        
    ## Calculate probability and duration IVT >= threshold
    thresholds = [150., 250., 500., 750.]
    probability_lst = []
    duration_lst = []
    for i, thres in enumerate(thresholds):
        data_size = ds.IVT.ensemble.shape
        tmp = ds.where(ds.IVT >= thres) # find where IVT exceeds threshold
        
        # sum the number of ensembles where IVT exceeds threshold
        probability = tmp.count(dim='ensemble') / data_size
        probability_lst.append(probability.IVT)
        
        # sum the number of time steps where IVT exceeds threshold
        duration = tmp.count(dim='forecast_hour')*3 # three hour time intervals
        duration_lst.append(duration.IVT)
        
    # merge duration and probability datasets
    duration_ds = xr.concat(duration_lst, pd.Index(thresholds, name="threshold"))
    prob_ds = xr.concat(probability_lst, pd.Index(thresholds, name="threshold"))
    
    ## Calculate Vectors
    uvec = ds.uIVT.mean('ensemble') # get the ensemble mean uIVT
    vvec = ds.vIVT.mean('ensemble') # get the ensemble mean vIVT
    ensemble_mean = ds.IVT.mean(dim='ensemble') # get the ensemble mean IVT

    # normalize vectors
    u = uvec / ensemble_mean
    v = vvec / ensemble_mean
    
    control = ds.IVT.sel(ensemble=0)
    
    ## place into single dataset
    x1 = ensemble_mean.rename('ensemble_mean')
    x2 = control.rename('control')
    x3 = v.rename('v')
    x4 = u.rename('u')
    x5 = duration_ds.rename('duration')
    x6 = prob_ds.rename('probability')
    
    x_lst = [x1, x2, x3, x4, x5]
    final_ds = xr.merge(x_lst)

    return final_ds, ds1