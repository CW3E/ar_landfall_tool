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

class download_QPF_prec_dataset:
    '''
    Copies or downloads the latest 7-d QPF for the current model initialization time
    
    Parameters
    ----------
    forecast : str
        name of the forecast product - options include GEFS or ECMWF
  
    Returns
    -------
    grb : grb file
        grb file downloaded for current run of figures
    
    '''
    def __init__(self, forecast):
        path_to_data = '/data/downloaded/SCRATCH/cw3eit_scratch/'
        self.forecast = forecast
        if forecast == 'GEFS':
            self.fpath = path_to_data + 'GEFS/FullFiles/'
            self.ensemble_name = 'GEFS'
        elif forecast == 'ECMWF':
            self.fpath = path_to_data + 'ECMWF/archive/' # will need to adjust when operational
            self.ensemble_name = 'ECMWF'
        else:
            print('Forecast product not available! Please choose either GEFS or ECMWF.')
                  
        ## find the most recent file in the currect directory
        list_of_files = glob.glob(self.fpath+'*.nc')
        self.fname = max(list_of_files, key=os.path.getctime)
        # pull the initialization date from the filename
        regex = re.compile(r'\d+')
        date_string = regex.findall(self.fname)
        self.date_string = date_string[1]
        self.model_init_date = datetime.datetime.strptime(self.date_string, '%Y%m%d%H')
        
    def download_dataset(self):
        date = self.model_init_date.strftime('%Y%m%d') # model init date
        hr = self.model_init_date.strftime('%H') # model init hour

        if self.forecast == 'GEFS':
            ## download from NOMADS
            print(date, hr)
            subprocess.check_call(["download_QPF.sh", date, hr], shell=True) # downloads the latest QPF data
        else:
            mmdyhr_init = self.model_init_date.strftime('%m%d%H') # month-day-hr init date
            date1 = datetime.datetime.strptime(self.date_string, '%Y%m%d%H')
            date2 = date1 + datetime.timedelta(days=7) 
            date2 = date2.strftime('%m%d%H') # valid date
            fpath = '/data/downloaded/Forecasts/ECMWF/NRT_data/{0}{1}/'.format(date, hr)
            fname = 'S1D{0}00{1}001'.format(mmdyhr_init, date2)
            shutil.copy(fpath+fname, 'precip_ECMWF') # copy file over to data folder

    
