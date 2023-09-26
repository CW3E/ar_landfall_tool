#!/usr/bin/python3
"""
Filename:    cw3e_tools.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: a collection of useful constants, colors, useful plotting and loading etc.
"""

import os
import shutil
import subprocess
import xarray as xr
import datetime
import numpy as np
import cartopy.crs as ccrs
import cmocean.cm as cmo
from PIL import Image
import requests

ivt_colors = {'250': (255./255.0, 174./255.0, 0./255.0), # orange
              '500': (236./255.0, 0./255.0, 7./255.0), # red 
              '750': (86./255.0, 0./255.0, 137./255.0), # purple
              '1000': (77./255.0, 77./255.0, 77./255.0) # grey
              }

def plot_terrain(ax, ext):
    fname = '/work/bkawzenuk_work/Maps/data/ETOPO1_Bed_c_gmt4.grd'
    datacrs = ccrs.PlateCarree()
    grid = xr.open_dataset(fname)
    grid = grid.where(grid.z > 0) # mask below sea level
    grid = grid.sel(x=slice(ext[0], ext[1]), y=slice(ext[2], ext[3]))
    cs = ax.pcolormesh(grid.x, grid.y, grid.z,
                        cmap=cmo.gray_r, transform=datacrs, alpha=0.7)
    
    return ax

def plot_cw3e_logo(ax, orientation):
    ## location of CW3E logo
    if orientation == 'horizontal':
        im = '/common/CW3E_Logo_Suite/1-Horzontal-PRIMARY_LOGO/Digital/JPG-RGB/CW3E-Logo-Horizontal-FullColor-RGB.jpg'
    else:
        im = '/common/CW3E_Logo_Suite/2-Vertical/Digital/JPG-RGB/CW3E-Logo-Vertical-FullColor-RGB.jpg'
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

def load_prec_QPF_dataset(forecast, model_init_date, date_string):
    date = model_init_date.strftime('%Y%m%d') # model init date
    hr = model_init_date.strftime('%H') # model init hour
    
    if forecast == 'GEFS':
        ## download from NOMADS
        subprocess.check_call(["../preprocess/download_QPF.sh", date, hr], shell=True) # downloads the latest QPF data
        ds = xr.open_dataset('../preprocess/precip.grb', engine='cfgrib')
        ds = ds.rename({'longitude': 'lon', 'latitude': 'lat'})
        prec = ds['tp']/25.4 # convert from mm to inches
    else:
        mmdyhr_init = model_init_date.strftime('%m%d%H') # month-day-hr init date
        date1 = datetime.datetime.strptime(date_string, '%Y%m%d%H')
        date2 = date1 + datetime.timedelta(days=7) 
        date2 = date2.strftime('%m%d%H') # valid date
        fpath = '/data/downloaded/Forecasts/ECMWF/NRT_data/{0}{1}/'.format(date, hr)
        fname = 'S1D{0}00{1}001'.format(mmdyhr_init, date2)
        var_lst = ['u10','lsm','msl','d2m','z','t2m','stl1', 'stl2', 'stl3', 'stl4', 'swvl4','swvl2', 'swvl3','sst','sp','v10','sd','skt', 'swvl1','siconc','tcwv','tcw']
        shutil.copy(fpath+fname, '../data/precip') # copy file over to data folder
        ds = xr.open_dataset('../data/precip', drop_variables=var_lst, engine='cfgrib', backend_kwargs={'filter_by_keys': {'typeOfLevel': 'surface'}})
        prec = ds['tp']*39.3701 # convert from m to inches
        prec = prec.rename({'longitude': 'lon', 'latitude': 'lat'}) # need to rename this to match GEFS
    return prec