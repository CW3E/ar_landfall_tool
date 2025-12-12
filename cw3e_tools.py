#!/usr/bin/python3
"""
Filename:    cw3e_tools.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: a collection of useful constants, colors, useful plotting and loading etc.
"""

import os
import re
import shutil
import subprocess
import glob
import xarray as xr
import pandas as pd
import datetime
import numpy as np
import cartopy.crs as ccrs
import cmocean.cm as cmo
from PIL import Image
from matplotlib import font_manager as fm
import matplotlib.pyplot as plt
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor

ivt_colors = {'250': (255./255.0, 174./255.0, 0./255.0), # orange
              '500': (236./255.0, 0./255.0, 7./255.0), # red
              '750': (86./255.0, 0./255.0, 137./255.0), # purple
              '1000': (77./255.0, 77./255.0, 77./255.0) # grey
              }

def plot_terrain(ax, ext):
    fname = '/data/projects/operations/data/ETOPO1_Bed_c_gmt4.grd'
    datacrs = ccrs.PlateCarree()
    grid = xr.open_dataset(fname)
    grid = grid.where(grid.z > 0) # mask below sea level
    grid = grid.sel(x=slice(ext[0], ext[1]), y=slice(ext[2], ext[3]))
    cs = ax.pcolormesh(grid.x, grid.y, grid.z,
                        cmap=cmo.gray_r, transform=datacrs, alpha=0.7)

    return ax

def set_cw3e_font(current_dpi, scaling_factor):
    fm.fontManager.addfont('/data/projects/operations/LandfallTools/ar_landfall_tool/utils/fonts/helvetica.ttc')

    plt.rcParams.update({
                    'font.family' : 'Helvetica',
                    'figure.dpi': current_dpi,
                    # 'font.size': 8 * scaling_factor, #changes axes tick label
                    # 'axes.labelsize': 8 * scaling_factor,
                    # 'axes.titlesize': 8 * scaling_factor,
                    # 'xtick.labelsize': 8 * scaling_factor,#do nothing
                    # 'ytick.labelsize': 8 * scaling_factor, #do nothing
                    # 'legend.fontsize': 5 * scaling_factor,
                    # 'lines.linewidth': 0.7 * scaling_factor,
                    # 'axes.linewidth': 0.2 * scaling_factor,
                    # 'legend.fontsize': 12 * scaling_factor,
                    # 'xtick.major.width': 0.8 * scaling_factor,
                    # 'ytick.major.width': 0.8 * scaling_factor,
                    # 'xtick.minor.width': 0.6 * scaling_factor,
                    # 'ytick.minor.width': 0.6 * scaling_factor,
                    # 'lines.markersize': 6 * scaling_factor
                })

def plot_cw3e_logo(ax, orientation):
    ## location of CW3E logo
    if orientation == 'horizontal':
        im = '/data/projects/operations/data/CW3E_Logo_Suite/1-Horzontal-PRIMARY_LOGO/Digital/JPG-RGB/CW3E-Logo-Horizontal-FullColor-RGB.jpg'
    else:
        im = '/data/projects/operations/data/CW3E_Logo_Suite/5-Vertical-Acronym_Only/Digital/PNG/CW3E-Logo-Vertical-Acronym-FullColor.png'
    img = np.asarray(Image.open(im))
    ax.imshow(img)
    ax.axis('off')
    return ax

def get_every_other_vector(x):
    '''
    stagger matrix setting values to diagonal
    based on https://www.w3resource.com/python-exercises/numpy/basic/numpy-basic-exercise-30.php

    Parameters
    ----------
    x : 2-D array

    Returns
    -------
    x : 2-D array
    same array as input but with the values staggered
    [[ 1.  0.  1.  0.]
     [ 0.  1.  0.  1.]
     [ 1.  0.  1.  0.]
     [ 0.  1.  0.  1.]]
    '''
    x[::2, 1::2] = 0
    x[1::2, ::2] = 0

    return x

def myround(x, base=5):
    return base * round(x/base)

class load_datasets:
    '''
    Copies or downloads the latest 7-d QPF for the current model initialization time
    Loads 7-d QPF
    Loads and calculates IVT vars for plotting

    Parameters
    ----------
    forecast : str
        name of the forecast product
        this can be 'ECMWF', 'GFS', 'WWRF'
    loc : str
        name of the location for the .txt file with latitude and longitude of locations for analysis
        this can be 'US-west', 'AK', or 'SAK'
    ptloc : str
        name of the transect of .txt file with latitude and longitude of locations for analysis
        this file should have no header, with latitude, then longitude, separated by a space
        this can be 'coast', 'foothills', or 'inland'

    Returns
    -------
    grb : grb file
        grb file downloaded for current run of figures

    '''
    def __init__(self, model, loc, ptloc, init_date):
        self.path_to_out = '/data/projects/operations/LandfallTools/ar_landfall_tool/data/'
        self.forecast = model
        self.model_init_date = init_date
        self.ptloc = ptloc
        self.loc = loc
        
        if self.forecast == 'GEFS':
            self.fpath = '/data/projects/derived_products/GEFS_IVT/data/'
            self.fname = f'{self.fpath}GEFS_IVT_{self.model_init_date}.nc'
            self.ensemble_name = 'GEFS'
            self.datasize_min = 15.
            
            
        elif self.forecast == 'ECMWF' or self.forecast == 'ECMWF-GEFS':
            self.fpath = '/data/projects/derived_products/ECMWF_IVT/Ensemble/'
            self.fname = f'{self.fpath}IVT_EC_{self.model_init_date}.nc'
            self.ensemble_name = 'ECMWF'
            self.datasize_min = 25.
            
        elif self.forecast == 'W-WRF':
            self.fpath = '/data/downloaded/WWRF-NRT/2025-2026/Ensemble_IVT/'
            self.fname = f'{self.fpath}IVT_WWRF_{self.model_init_date}.nc'
            self.ensemble_name = 'West-WRF'
            self.datasize_min = 50.
        else:
            print('Forecast product not available! Please choose either GEFS, ECMWF, ECMWF-GEFS, or W-WRF.')

    def download_QPF_dataset(self):
        dt_init = datetime.strptime(self.model_init_date, "%Y%m%d%H")
        date = dt_init.strftime('%Y%m%d') # model init date
        hr = dt_init.strftime('%H') # model init hour

        if self.forecast == 'GEFS':
            ## download from NOMADS
            print(date, hr)
            subprocess.check_call(["download_QPF.sh", date, hr], shell=True) # downloads the latest QPF data
        else:
            mmdyhr_init = dt_init.strftime('%m%d%H') # month-day-hr init date
            date2 = dt_init + datetime.timedelta(days=7)
            date2 = date2.strftime('%m%d%H') # valid date
            fpath = '/data/downloaded/Forecasts/ECMWF/NRT_data/{0}{1}/'.format(date, hr)
            fname = 'S1D{0}00{1}001'.format(mmdyhr_init, date2)
            shutil.copy(fpath+fname, self.path_to_out+'precip_ECMWF') # copy file over to data folder


    def load_prec_QPF_dataset(self):

        if self.forecast == 'GEFS':
            try:
                ## this method directly opens data from NOMADS
                print(self.model_init_date)
                dt_init = datetime.strptime(self.model_init_date, "%Y%m%d%H")
                date = dt_init.strftime('%Y%m%d') # model init date
                hr = dt_init.strftime('%H') # model init hour
                url = 'https://nomads.ncep.noaa.gov/dods/gfs_0p25/gfs{0}/gfs_0p25_{1}z'.format(date, hr)
                ds = xr.open_dataset(url, decode_times=False)
                ds = ds.isel(time=7*8) # get 7-day QPF - the variable is already cumulative
                prec = ds['apcpsfc']/25.4 # convert from mm to inches
            except OSError:
                try:
                    ## This method uses the downloaded data
                    self.download_QPF_dataset()
                    ds = xr.open_dataset('precip_GFS.grb', engine='cfgrib', backend_kwargs={"indexpath": ""})
                    ds = ds.rename({'longitude': 'lon', 'latitude': 'lat'})
                    prec = ds['tp']/25.4 # convert from mm to inches
                except OSError:
                    ## build a fake precip dataset of 0s
                    lats = np.arange(-90., 90.25, .25)
                    lons = np.arange(-180., 180.25, .25)
                    var_dict = {'tp': (['lat', 'lon'], np.zeros((len(lats), len(lons))))}
                    prec = xr.Dataset(var_dict, coords={'lat': (['lat'], lats),
                                                        'lon': (['lon'], lons)})

        else:
            self.download_QPF_dataset()
            var_lst = ['u10','lsm','msl','d2m','z','t2m','stl1', 'stl2', 'stl3', 'stl4', 'swvl4','swvl2', 'swvl3','sst','sp','v10','sd','skt', 'swvl1','siconc','tcwv','tcw']
            ds = xr.open_dataset(self.path_to_out+'precip_ECMWF', drop_variables=var_lst, engine='cfgrib', backend_kwargs={'filter_by_keys': {'typeOfLevel': 'surface'}})
            prec = ds['tp']*39.3701 # convert from m to inches
            prec = prec.rename({'longitude': 'lon', 'latitude': 'lat'}) # need to rename this to match GEFS

        return prec

    def get_lat_lons_from_txt_file(self):
        ## read text file with points
        textpts_fname = self.path_to_out+'{0}/latlon_{1}.txt'.format(self.loc, self.ptloc)
        df = pd.read_csv(textpts_fname, header=None, sep=' ', names=['latitude', 'longitude'], engine='python')
        df['longitude'] = df['longitude']*-1
        df = df
        self.lons = df['longitude'].values
        self.lats = df['latitude'].values

    def calc_ivt_vars(self):

        def background_ivt_calculation(ds):
            ds1 = ds.sel(forecast_hour=slice(0, 24*7)).astype('float64').mean(['forecast_hour', 'ensemble'])
            ds1 = ds1.astype('float32')
            ds1 = ds1.where(ds1.IVT >= 250)
            return ds1

        # Load and preprocess dataset
        print('Reading IVT data ...')
        ds = xr.open_dataset(self.fname)
        ds = ds.astype('float16')
        for coord in ds.coords:
            if np.issubdtype(ds.coords[coord].dtype, np.floating):
                ds.coords[coord] = ds.coords[coord].astype('float16')
        ds = ds.assign_coords(lon=((ds.lon + 180) % 360 - 180)).sortby('lon')
        
        ## updates specific to model name
        if self.forecast == 'ECMWF':
            ds = ds.rename({'forecast_time': 'forecast_hour'})
            dt_init = datetime.strptime(self.model_init_date, "%Y%m%d%H")
            if int(dt_init.strftime('%Y')) > 2020:
                ds['forecast_hour'] = ds.forecast_hour * 3
            
        elif self.forecast == 'W-WRF':
            ds = ds.rename({'ensembles': 'ensemble'})

        self.get_lat_lons_from_txt_file()
        x = xr.DataArray(self.lons, dims=['location'])
        y = xr.DataArray(self.lats, dims=['location'])
        ds = ds.sel(lon=x, lat=y, method='nearest')
        
        executor = concurrent.futures.ThreadPoolExecutor()
        future = executor.submit(background_ivt_calculation, ds)

        ## Calculate probability and duration IVT >= threshold
        thresholds = [100, 150, 250, 500, 750, 1000]
        # Precompute data_size and valid_mask once
        ivt = ds.IVT  # shape: (forecast_hour, ensemble, location)
        ens_size = ivt.sizes['ensemble']

        data_size = ivt.count(dim='ensemble')  # (forecast_hour, location)
        valid_mask = data_size >= self.datasize_min

        # Prepare result lists
        probability_lst = []
        duration_lst = []

        def process_threshold(thres):
            mask = ivt >= thres

            count_ens = mask.sum(dim='ensemble') / ens_size
            count_ens = count_ens.where(valid_mask)

            count_time = mask.sum(dim='forecast_hour') * 3
            count_time = count_time.where(valid_mask)

            return count_ens, count_time

        with ThreadPoolExecutor() as executor:
            results = list(executor.map(process_threshold, thresholds))

        probability_lst, duration_lst = zip(*results)

        # merge duration and probability datasets
        duration_ds = xr.concat(duration_lst, pd.Index(thresholds, name="threshold"))
        prob_ds = xr.concat(probability_lst, pd.Index(thresholds, name="threshold"))

        ## Calculate Vectors
        print('Computing ensemble mean...')
        # get the ensemble mean vIVT and uIVT
        uvec = (
            ds.uIVT.astype('float64')
            .where(data_size >= self.datasize_min)
            .mean(dim='ensemble')
        )
        vvec = (
            ds.vIVT.astype('float64')
            .where(data_size >= self.datasize_min)
            .mean(dim='ensemble')
        )
        # get the ensemble mean IVT where there are enough ensembles
        ensemble_mean = (
            ds.IVT.astype('float64')
            .where(data_size >= self.datasize_min)
            .mean(dim='ensemble')
        )

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

        x_lst = [x1, x2, x3, x4, x5, x6]
        final_ds = xr.merge(x_lst)

        ##Add attribute information
        final_ds = final_ds.assign_attrs(model_init_date=self.model_init_date)
        ds1 = future.result()
        return final_ds, ds1