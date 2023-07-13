"""
Filename:    ar_tools.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Functions for CW3E AR tools
"""
# Standard Python modules
import os, sys
import glob
import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime
import re
import textwrap
from PIL import Image

# matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colorbar import Colorbar # different way to handle colorbar
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
from matplotlib import dates as mdates

# cartopy
import cartopy.crs as ccrs
from cartopy.mpl.geoaxes import GeoAxes
import cartopy.feature as cfeature

# other
import cmocean.cm as cmo

# import personal modules
import nclcmaps as nclc
import cw3ecmaps as cw3e


def plot_terrain(ax, ext):
    fname = '/work/bkawzenuk_work/Maps/data/ETOPO1_Bed_c_gmt4.grd'
    datacrs = ccrs.PlateCarree()
    grid = xr.open_dataset(fname)
    grid = grid.where(grid.z > 0) # mask below sea level
    grid = grid.sel(x=slice(ext[0], ext[1]), y=slice(ext[2], ext[3]))
    cs = ax.pcolormesh(grid.x, grid.y, grid.z,
                        cmap=cmo.gray_r, transform=datacrs, alpha=0.7)
    
    return ax

def plot_cw3e_logo(ax):
    ## location of CW3E logo
    im = '/common/CW3E_Logo_Suite/1-Horzontal-PRIMARY_LOGO/Digital/JPG-RGB/CW3E-Logo-Horizontal-FullColor-RGB.jpg'
    img = np.asarray(Image.open(im))
    ## TODO make logo bigger
    # ax.imshow(img, extent=[-0.25, 1.25, -0.5, 1.5], aspect='auto')
    ax.imshow(img)
    ax.axis('off')
    # ax.set_aspect(0.6)
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

def load_prec_QPF_dataset(forecast, model_init_date):

    if forecast == 'GEFS':
        date = model_init_date.strftime('%Y%m%d') # model init date
        hr = model_init_date.strftime('%H') # model init hour
        url = 'https://nomads.ncep.noaa.gov/dods/gfs_0p25/gfs{0}/gfs_0p25_{1}z'.format(date, hr)
        ds = xr.open_dataset(url, decode_times=False)
        # subprocess.check_call(["../preprocess/download_QPF.sh", date, hr], shell=True) # downloads the latest QPF data
        ds = ds.isel(time=7*8) # get 7-day QPF - the variable is already cumulative
        prec = ds['apcpsfc']/25.4 # convert from mm to inches
    else:
        print('No precipitation data yet')
        
    return prec

class landfall_tool_contour:
    '''
    Returns a .png file with Cordeira AR Landfall Tool (contour) figure with input locations, chosen Forecast product, and IVT threshold
    
    Parameters
    ----------
    ptloc : str
        name of the .txt file with latitude and longitude of locations for analysis
        this file should have no header, with latitude, then longitude, separated by a space
    forecast : str
        name of the forecast product - options include GEFS, ECMWF, TODO: ECMWF - GEFS and West-WRF
    threshold : int
        threshold for IVT probabilty in kg m-1 s-1 - options include 150, 250, 500, 750
  
    Returns
    -------
    fig : figure
        png file of the figure
    
    '''
    
    def __init__(self, loc, ptloc, forecast='GEFS', threshold=250, orientation='latitude'):
        path_to_data = '/data/downloaded/SCRATCH/cw3eit_scratch/'
        if forecast == 'GEFS':
            fpath = path_to_data + 'GEFS/FullFiles/'
            self.ensemble_name = 'GEFSv12'
        elif forecast == 'ECMWF':
            fpath = path_to_data + 'ECMWF/'
            self.ensemble_name = 'ECMWF EPS'
        else:
            print('Forecast product not available! Please choose either GEFS or ECMWF.')
                  
        ## find the most recent file in the currect directory
        list_of_files = glob.glob(fpath+'*.nc') 
        self.fname = max(list_of_files, key=os.path.getctime)
        # pull the initialization date from the filename
        regex = re.compile(r'\d+')
        date_string = regex.findall(self.fname)
        self.date_string = date_string[1]
        self.model_init_date = datetime.strptime(self.date_string, '%Y%m%d%H')

        ## read text file with points
        self.loc = loc
        self.ptloc = ptloc
        textpts_fname = '../data/{0}/latlon_{1}.txt'.format(self.loc, self.ptloc)
        df = pd.read_csv(textpts_fname, header=None, sep=' ', names=['latitude', 'longitude'], engine='python')
        df['longitude'] = df['longitude']*-1
        self.df = df
        self.threshold = threshold
                  
        ## format dicts for plots
        self.kw_ticklabels = {'size': 6, 'color': 'gray', 'fontweight': 'light'}
        self.kw_grid = {'linewidth': .5, 'color': 'k', 'linestyle': '--', 'alpha': 0.1}
        self.kw_ticks = {'length': 4, 'width': 0.5, 'pad': 2, 'color': 'black'}
        self.IVT_units = 'kg m$^{-1}$ s$^{-1}$'
        
        ## info for plot orientation
        self.orientation = orientation

    def get_date_information(self):

        hr = self.model_init_date.strftime('%H')
        weekday = self.model_init_date.strftime('%a')
        day = self.model_init_date.strftime('%-d')
        month = self.model_init_date.strftime('%b')
        year = self.model_init_date.strftime('%Y')
        self.xlbl = "<-------- Forecast Day from {0}Z on {1} {2} {3} {4} --------".format(hr, weekday, day, month, year)
        self.title = 'Model Run: {0}Z {1} {2} {3} {4}'.format(hr, weekday, day, month, year)
                  
        ## create datetime labels for the x-axis
        date_lst = pd.date_range(self.model_init_date, periods=17, freq='1D')
        xtck_lbl = []
        for i, x in enumerate(date_lst):
            t = pd.to_datetime(str(x))
            xtck_lbl.append(t.strftime('%m/%d'))
        
        self.xtck_lbl = xtck_lbl
        
    def calc_ivt_probability(self):
        ## load the forecast data
        ds = xr.open_dataset(self.fname)
        ds = ds.assign_coords({"lon": (((ds.lon + 180) % 360) - 180)}) # Convert DataArray longitude coordinates from 0-359 to -180-179
        ds = ds.sel(lon=slice(-180., -1)) # keep only western hemisphere
                  
        # subset ds to the select points
        self.lons = self.df['longitude'].values
        self.lats = self.df['latitude'].values
        x = xr.DataArray(self.lons, dims=['location'])
        y = xr.DataArray(self.lats, dims=['location'])

        ds = ds.sel(lon=x, lat=y, method='nearest')

        ## calculate probability IVT >= threshold
        data_size = ds.IVT.ensemble.shape

        # sum the number of ensembles where IVT exceeds threshold
        self.probability = (ds.IVT.where(ds.IVT >= self.threshold).count(dim='ensemble')) / data_size
    
        
    def plot_probability_latitude(self, ax):
        
        self.calc_ivt_probability() # run the calculation
        
        # Contour Filled
        x = self.probability.forecast_hour / 24 # convert forecast hour to forecast day
        y = self.probability.lat
        data = np.flipud(np.rot90(self.probability.values)) # rotate data 90 degrees and flip up down
        cmap, norm, bnds = cw3e.cmap('ivt_probability')
        self.cflevs = bnds
        self.cf = ax.pcolormesh(x, y, data, levels=self.cflevs, cmap=cmap, norm=norm, rasterized=True, extend='neither')
        ax.invert_xaxis() # invert x-axis so that time reads from right to left
                  
        # get tick and label information
        self.get_date_information()

        # apply ytick parameters
        ticks_loc = ax.get_yticks().tolist()
        ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
        ax.set_yticklabels([u"{:0.0f}\N{DEGREE SIGN}N".format(x) for x in ticks_loc], **self.kw_ticklabels)
        # apply xtick parameters
        ax.xaxis.set_major_locator(mticker.MaxNLocator(16))
        ticks_loc = ax.get_xticks().tolist()
        ax.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
        # ax.set_xticklabels([u"{:0.0f}".format(x) for x in ticks_loc], **self.kw_ticklabels) # labels are days since forecast initialization
        ax.set_xticklabels(["{0}".format(x) for x in self.xtck_lbl], **self.kw_ticklabels) # labels are month/day

        # apply gridlines
        ax.minorticks_on()
        ax.grid(visible=None, which='both', axis='y', **self.kw_grid)
        ax.grid(visible=None, which='major', axis='x', **self.kw_grid)
        ax.tick_params(axis='x', which='minor', bottom=False)
        ax.tick_params(axis='x', which='major', **self.kw_ticks)
        ax.tick_params(axis='y', which='major', direction='out', **self.kw_ticks)

        ## labels and subtitles
        ax.set_ylabel("Latitude along West Coast", fontsize=8)
        ax.set_xlabel(self.xlbl, fontsize=8)
        ax.set_title(self.title, loc='right', fontsize=8)
        ax.set_title('16-d {0} Prob of IVT > {1} {2}'.format(self.ensemble_name, self.threshold, self.IVT_units), loc='left', fontsize=8)
        
        return ax
    
    def plot_probability_longitude(self, ax):
        
        self.calc_ivt_probability() # run the calculation
        
        # Contour Filled
        y = self.probability.forecast_hour / 24 # convert forecast hour to forecast day
        x = self.probability.lon
        data = self.probability.values
        cmap, norm, bnds = cw3e.cmap('ivt_probability')
        self.cflevs = bnds
        self.cf = ax.pcolormesh(x, y, data, levels=self.cflevs, cmap=cmap, norm=norm, rasterized=True, extend='neither')
        ax.invert_xaxis() # invert x-axis so that time reads from right to left
                  
        # get tick and label information
        self.get_date_information()

        # apply xtick parameters
        ticks_loc = ax.get_xticks().tolist()
        ax.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
        ax.set_xticklabels([u"{:0.0f}\N{DEGREE SIGN}W".format(x) for x in ticks_loc], **self.kw_ticklabels)
        # apply ytick parameters
        ax.yaxis.set_major_locator(mticker.MaxNLocator(16))
        ticks_loc = ax.get_yticks().tolist()
        ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
        # ax.set_yticklabels([u"{:0.0f}".format(x) for x in ticks_loc], **self.kw_ticklabels) # labels are days since forecast initialization
        ax.set_yticklabels(["{0}".format(x) for x in self.xtck_lbl], **self.kw_ticklabels) # labels are month/day

        # apply gridlines
        ax.minorticks_on()
        ax.grid(visible=None, which='both', axis='x', **self.kw_grid)
        ax.grid(visible=None, which='major', axis='y', **self.kw_grid)
        ax.tick_params(axis='y', which='minor', left=False)
        ax.tick_params(axis='y', which='major', **self.kw_ticks)
        ax.tick_params(axis='x', which='major', direction='out', **self.kw_ticks)

        ## labels and subtitles
        # ax.set_xlabel("Longitude along West Coast", fontsize=8)
        ax.set_ylabel(self.xlbl, fontsize=8)
        ax.set_title(self.title, loc='right', fontsize=8)
        ax.set_title('16-d {0} Prob of IVT $\geq$ {1} {2}'.format(self.ensemble_name, self.threshold, self.IVT_units), loc='left', fontsize=8)
        
        plt.gca().invert_yaxis()
        plt.gca().invert_xaxis()
        
        return ax
    
                  
    def plot_map(self, ax, mapcrs, datacrs):
        ## Set up extent of plot at gridlines
        if self.orientation == 'latitude':          
            lonmin = self.lons.min()-2.
            lonmax = self.lons.max()+2.
            latmin = self.lats.min()
            latmax = self.lats.max()
            
        elif self.orientation == 'longitude':
            lonmin = self.lons.min()
            lonmax = self.lons.max()
            latmin = self.lats.min()-2.
            latmax = self.lats.max()+2.
            
        ext = [lonmin, lonmax, latmin, latmax] # extent of map plot
        du = 5 # how frequent for ticks
        dx = np.arange(lonmin,lonmax+du,du)
        dy = np.arange(latmin,latmax+du,du)
                  
        ax.set_extent(ext, crs=datacrs)
        
        ## Add elevation contours
        ax = plot_terrain(ax, ext)

        # Add map features (continents and country borders)
        # ax.add_feature(cfeature.LAND, facecolor='0.9')      
        ax.add_feature(cfeature.BORDERS, edgecolor='0.4', linewidth=0.4)
        ax.add_feature(cfeature.STATES, edgecolor='0.2', linewidth=0.2)
        ax.add_feature(cfeature.OCEAN, edgecolor='0.4', facecolor='lightskyblue', linewidth=0.2)

        # apply tick parameters   
        gl = ax.gridlines(crs=datacrs, draw_labels=True, **self.kw_grid)
        gl.top_labels = False
        gl.left_labels = True
        gl.right_labels = False
        gl.bottom_labels = True
        gl.xlocator = mticker.FixedLocator(dx)
        gl.ylocator = mticker.FixedLocator(dy)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = self.kw_ticklabels
        gl.ylabel_style = self.kw_ticklabels

        gl.xlines = True
        gl.ylines = True
         
        ## plot point locations
        for i, (x, y) in enumerate(zip(self.lons, self.lats)):
            ax.plot(x, y, 'ko', markersize=2, transform=datacrs)

        ## plot labels
        if self.loc == 'US-west':
            label = 'Forecasts support FIRO/CA-AR Program Intended for research purposes only'
            ax.set_title(textwrap.fill(label, 36), loc='right', fontsize=5)
        else:
            label = 'Forecasts support NSF Coastlines and People Program: #2052972  Intended for research purposes only'
            ax.set_title(textwrap.fill(label, 61), loc='right', fontsize=5)
        
        return ax
    
    def create_figure(self):
        fname = '../figs/{0}/landfall-contour_{1}_{2}_{3}_{4}'.format(self.loc, self.ptloc, self.ensemble_name, self.threshold, self.date_string)
        fmt = 'png'
        
        if self.orientation == 'latitude':
            fig = plt.figure(figsize=(8.5, 4))
            fig.dpi = 300
            nrows = 2
            ncols = 2
            ## Use gridspec to set up a plot with a series of subplots that is
            ## n-rows by n-columns
            gs = GridSpec(nrows, ncols, height_ratios=[1, 0.05], width_ratios = [2, 1], wspace=0.01, hspace=0.4)
            ## use gs[rows index, columns index] to access grids         

            ## Add probability plot         
            ax = fig.add_subplot(gs[0, 0])
            self.plot_probability_latitude(ax)

            ## Add color bar
            cbax = plt.subplot(gs[1,0]) # colorbar axis
            cb = Colorbar(ax = cbax, mappable = self.cf, orientation = 'horizontal', ticklocation = 'bottom', ticks=self.cflevs[::2])
            cb.ax.set_xticklabels(["{:.0%}".format(i) for i in cb.get_ticks()], **self.kw_ticklabels)  # horizontally oriented colorbar

            # Set up projection information for map
            mapcrs = ccrs.PlateCarree()
            datacrs = ccrs.PlateCarree()
            ## Add map
            ax = fig.add_subplot(gs[0, 1], projection=mapcrs)
            self.plot_map(ax, mapcrs, datacrs)

            ## Add CW3E logo
            ax = fig.add_subplot(gs[1, 1])
            ax = plot_cw3e_logo(ax)
            
        
        elif self.orientation == 'longitude':
            fig = plt.figure(figsize=(8, 8))
            fig.dpi = 300
            nrows = 3
            ncols = 1
            ## Use gridspec to set up a plot with a series of subplots that is
            ## n-rows by n-columns
            gs = GridSpec(nrows, ncols, height_ratios=[0.7, 1, 0.05], width_ratios = [1], hspace=0.1)
            ## use gs[rows index, columns index] to access grids         

            ## Add probability plot         
            ax = fig.add_subplot(gs[1, 0])
            self.plot_probability_longitude(ax)

            ## Add color bar
            cbax = plt.subplot(gs[2,0]) # colorbar axis
            cb = Colorbar(ax = cbax, mappable = self.cf, orientation = 'horizontal', ticklocation = 'bottom', ticks=self.cflevs[::2])
            cb.ax.set_xticklabels(["{:.0%}".format(i) for i in cb.get_ticks()], **self.kw_ticklabels)  # horizontally oriented colorbar

            # Set up projection information for map
            mapcrs = ccrs.PlateCarree()
            datacrs = ccrs.PlateCarree()
            ## Add map
            ax = fig.add_subplot(gs[0, 0], projection=mapcrs)
            self.plot_map(ax, mapcrs, datacrs)
            
            ## Add CW3E logo
            ax = fig.add_subplot(gs[1, 1])
            ax = plot_cw3e_logo(ax)
            
        
        fig.savefig('%s.%s' %(fname, fmt), bbox_inches='tight', dpi=fig.dpi)
        plt.show()
        
        
class landfall_tool_vector:
    '''
    Returns a .png file with Cordeira vector landfall tool figure with input locations, chosen Forecast product, and IVT threshold
    
    Parameters
    ----------
    ptloc : str
        name of the .txt file with latitude and longitude of locations for analysis
        this file should have no header, with latitude, then longitude, separated by a space
    forecast : str
        name of the forecast product - options include GEFS, ECMWF, TODO: ECMWF - GEFS and West-WRF
    threshold : int
        threshold for IVT probabilty in kg m-1 s-1 - options include 150, 250, 500, 750
    orientation: str
        orientation of the probability (e.g., latitude or longitude)
  
    Returns
    -------
    fig : figure
        png file of the vector landfall tool
    
    '''
    
    def __init__(self, loc, ptloc, forecast='GEFS', threshold=250, orientation='latitude'):
        path_to_data = '/data/downloaded/SCRATCH/cw3eit_scratch/'
        self.forecast = forecast
        if self.forecast == 'GEFS':
            fpath = path_to_data + 'GEFS/FullFiles/'
            self.ensemble_name = 'GEFSv12'
        elif self.forecast == 'ECMWF':
            fpath = path_to_data + 'ECMWF/'
            self.ensemble_name = 'ECMWF EPS'
        else:
            print('Forecast product not available! Please choose either GEFS or ECMWF.')
                  
        ## find the most recent file in the currect directory
        list_of_files = glob.glob(fpath+'*.nc') 
        self.fname = max(list_of_files, key=os.path.getctime)
        # pull the initialization date from the filename
        regex = re.compile(r'\d+')
        date_string = regex.findall(self.fname)
        self.date_string = date_string[1]
        self.model_init_date = datetime.strptime(self.date_string, '%Y%m%d%H')

        ## read text file with points
        self.loc = loc
        self.ptloc = ptloc
        textpts_fname = '../data/{0}/latlon_{1}.txt'.format(self.loc, self.ptloc)
        df = pd.read_csv(textpts_fname, header=None, sep=' ', names=['latitude', 'longitude'], engine='python')
        df['longitude'] = df['longitude']*-1
        self.df = df
        self.threshold = threshold
                  
        ## format dicts for plots
        self.kw_ticklabels = {'fontsize': 6, 'color': 'gray', 'fontweight': 'light'}
        self.kw_grid = {'linewidth': .5, 'color': 'k', 'linestyle': '--', 'alpha': 0.1}
        self.kw_ticks = {'length': 4, 'width': 0.5, 'pad': 2, 'color': 'black'}
        self.kw_quiver = {'headlength': 6, 'headaxislength': 4.5, 'headwidth': 4.5}
        self.IVT_units = 'kg m$^{-1}$ s$^{-1}$'
        
        ## info for plot orientation and ivt shading
        self.orientation = orientation

    def get_date_information(self):

        hr = self.model_init_date.strftime('%H')
        weekday = self.model_init_date.strftime('%a')
        day = self.model_init_date.strftime('%-d')
        month = self.model_init_date.strftime('%b')
        year = self.model_init_date.strftime('%Y')
        self.xlbl = "<-------- Forecast Day from {0}Z on {1} {2} {3} {4} --------".format(hr, weekday, day, month, year)
        self.title = 'Model Run: {0}Z {1} {2} {3} {4}'.format(hr, weekday, day, month, year)
                  
        ## create datetime labels for the x-axis
        date_lst = pd.date_range(self.model_init_date, periods=17, freq='1D')
        xtck_lbl = []
        for i, x in enumerate(date_lst):
            t = pd.to_datetime(str(x))
            xtck_lbl.append(t.strftime('%m/%d'))
        
        self.xtck_lbl = xtck_lbl
        
    def load_dataset(self, subset=True):
        ## load the forecast data
        ds = xr.open_dataset(self.fname)
        ds = ds.assign_coords({"lon": (((ds.lon + 180) % 360) - 180)}) # Convert DataArray longitude coordinates from 0-359 to -180-179
        ds = ds.sel(lon=slice(-180., -1)) # keep only western hemisphere
        ds = ds.sel(forecast_hour=slice(0, 24*7)) # keep only the first seven days of the forecast
        
        if subset == True: # subset to the latitude points from txt file
            # subset ds to the select points
            self.lons = self.df['longitude'].values
            self.lats = self.df['latitude'].values
            x = xr.DataArray(self.lons, dims=['location'])
            y = xr.DataArray(self.lats, dims=['location'])
        
            ds = ds.sel(lon=x, lat=y, method='nearest')
        
        else:
            ds = ds.mean(['forecast_hour', 'ensemble'])
            ds = ds.where(ds.IVT >= 200)  # is there a threshold set for masking?           
        
        return ds
        
    def calc_ivt_duration(self):
        
        ds = self.load_dataset(subset=True)
        thresholds = [250, 500, 750, 1000]
        duration_lst = []
        for i, thres in enumerate(thresholds):
            tmp = ds.where(ds.IVT >= thres) # find where IVT exceeds threshold
            duration = tmp.count(dim='forecast_hour')*3 # three hour time intervals
            duration_lst.append(duration)
        
        # merge duration datasets
        duration_ds = xr.concat(duration_lst, pd.Index(thresholds, name="duration"))
        
        return duration_ds    
        
    def calc_ivt_probability(self):
        
        ds = self.load_dataset(subset=True)

        ## calculate probability IVT >= threshold
        data_size = ds.IVT.ensemble.shape

        # sum the number of ensembles where IVT exceeds threshold
        self.probability = (ds.IVT.where(ds.IVT >= self.threshold).count(dim='ensemble')) / data_size
        uvec = ds.uIVT.mean('ensemble') # get the ensemble mean uIVT
        vvec = ds.vIVT.mean('ensemble') # get the ensemble mean vIVT
        self.control = ds.IVT.mean('ensemble') # get the ensemble mean IVT
        
        # normalize vectors
        self.u = uvec / self.control 
        self.v = vvec / self.control
        
    def plot_vector_landfall_latitude(self, ax):
        
        self.calc_ivt_probability() # run the calculation
        x = self.probability.forecast_hour / 24 # convert forecast hour to forecast day
        y = self.probability.lat
        data = get_every_other_vector(np.flipud(np.rot90(self.probability.values))) # rotate data 90 degrees and flip up down
        uvec = get_every_other_vector(np.flipud(np.rot90(self.u.values)))
        vvec = get_every_other_vector(np.flipud(np.rot90(self.v.values)))
        ctrl = np.flipud(np.rot90(self.control.values))
        
        # Vectors   
        cmap, norm, bnds = cw3e.cmap('ivt_vector')
        self.cflevs = bnds
        self.cf = ax.quiver(x, y, uvec, vvec, data, cmap=cmap, norm=norm,
                            capstyle='round', units='width', **self.kw_quiver)
        
        ## Contours
        ## add contour lines of control ensemble IVT every 250 kg m-1 s-1
        clevs = np.arange(250, 3500, 250.)
        cs = ax.contour(x, y, ctrl, levels=clevs, colors='k', linewidths=0.5)
        
        ax.invert_xaxis() # invert x-axis so that time reads from right to left
                  
        # get tick and label information
        self.get_date_information()

        # apply ytick parameters
        nticks = 8
        ticks_loc = ax.get_yticks().tolist()
        ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
        ax.set_yticklabels([u"{:0.0f}\N{DEGREE SIGN}N".format(x) for x in ticks_loc], **self.kw_ticklabels)
        # apply xtick parameters
        ax.xaxis.set_major_locator(mticker.MaxNLocator(nticks))
        ticks_loc = ax.get_xticks().tolist()
        ax.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
        # ax.set_xticklabels([u"{:0.0f}".format(x) for x in ticks_loc], **self.kw_ticklabels) # labels are days since forecast initialization
        ax.set_xticklabels(["{0}".format(x) for x in self.xtck_lbl[:nticks]], **self.kw_ticklabels) # labels are month/day
        
        # labels are days since forecast initialization
        for i, x in enumerate(ticks_loc):
            xnew = mdates.date2num(x)
            ax.annotate(u"{:0.0f}".format(x), # this is the text
                       (xnew,0), # these are the coordinates to position the label
                        textcoords="offset points", # how to position the text
                        xytext=(0,0), # distance from text to points (x,y)
                        ha='center', # horizontal alignment can be left, right or center
                        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="k", lw=0.5, alpha=0.8),
                        # xycoords=transform,
                        zorder=200,
                        fontsize=6)

        # apply gridlines
        ax.minorticks_on()
        ax.grid(visible=None, which='both', axis='y', **self.kw_grid)
        ax.grid(visible=None, which='major', axis='x', **self.kw_grid)
        ax.tick_params(axis='x', which='minor', bottom=False)
        ax.tick_params(axis='x', which='major', **self.kw_ticks)
        ax.tick_params(axis='y', which='major', direction='out', **self.kw_ticks)

        ## labels and subtitles
        ax.set_ylabel("Latitude along West Coast", fontsize=5)
        ax.set_xlabel(self.xlbl, fontsize=5)
        ax.set_title(self.title, loc='right', fontsize=5)
        ax.set_title('7-d {0} Ens. Mean IVT colored by Prob of IVT $\geq$ {1} {2}'.format(self.ensemble_name, self.threshold, self.IVT_units), loc='left', fontsize=5)
        
        return ax
    
    def plot_vector_landfall_longitude(self, ax):
        
        self.calc_ivt_probability() # run the calculation
        x = self.probability.lon
        y = self.probability.forecast_hour / 24 # convert forecast hour to forecast day
        data = self.probability.values
        uvec = self.u.values
        vvec = self.v.values
        ctrl = self.control.values
        
        ## Quiver
        cmap, norm, bnds = cw3e.cmap('ivt_vector')
        self.cflevs = bnds
        self.cf = ax.quiver(x, y, u, v, data, cmap=cmap, norm=norm,
                            capstyle='round', units='width', facecolor='white')
        
        ## Contours
        ## add contour lines of control ensemble IVT every 250 kg m-1 s-1
        clevs = np.arange(250, 5000, 250)
        cs = ax.contour(x, y, ctrl, levels=clevs, colors='k', linewidths=0.5)
            
        ax.invert_xaxis() # invert x-axis so that time reads from right to left
                  
        # get tick and label information
        self.get_date_information()

        # apply xtick parameters
        nticks = 8
        ticks_loc = ax.get_xticks().tolist()
        ax.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
        ax.set_xticklabels([u"{:0.0f}\N{DEGREE SIGN}W".format(x) for x in ticks_loc], **self.kw_ticklabels)
        # apply ytick parameters
        ax.yaxis.set_major_locator(mticker.MaxNLocator(nticks))
        ticks_loc = ax.get_yticks().tolist()
        ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
        # ax.set_yticklabels([u"{:0.0f}".format(x) for x in ticks_loc], **self.kw_ticklabels) # labels are days since forecast initialization
        ax.set_yticklabels(["{0}".format(x) for x in self.xtck_lbl], **self.kw_ticklabels) # labels are month/day

        # apply gridlines
        ax.minorticks_on()
        ax.grid(visible=None, which='both', axis='x', **self.kw_grid)
        ax.grid(visible=None, which='major', axis='y', **self.kw_grid)
        ax.tick_params(axis='y', which='minor', left=False)
        ax.tick_params(axis='y', which='major', **self.kw_ticks)
        ax.tick_params(axis='x', which='major', direction='out', **self.kw_ticks)

        ## labels and subtitles
        # ax.set_xlabel("Longitude along West Coast", fontsize=8)
        ax.set_ylabel(self.xlbl, fontsize=5)
        ax.set_title(self.title, loc='right', fontsize=5)
        ax.set_title('7-d {0} Prob of IVT > {1} {2}'.format(self.ensemble_name, self.threshold, self.IVT_units), loc='left', fontsize=5)
        
        plt.gca().invert_yaxis()
        plt.gca().invert_xaxis()
        
        return ax
    
                  
    def plot_map(self, ax, mapcrs, datacrs):
        ## Set up extent of plot at gridlines
        if self.orientation == 'latitude':          
            lonmin = self.lons.min()-2.
            lonmax = self.lons.max()+2.
            latmin = self.lats.min()
            latmax = self.lats.max()
            
        elif self.orientation == 'longitude':
            lonmin = self.lons.min()
            lonmax = self.lons.max()
            latmin = self.lats.min()-2.
            latmax = self.lats.max()+2.
            
        ext = [lonmin, lonmax, latmin, latmax] # extent of map plot
        du = 5 # how frequent for ticks
        dx = np.arange(lonmin,lonmax+du,du)
        dy = np.arange(latmin,latmax+du,du)
                  
        ax.set_extent(ext, crs=datacrs)
        
        ## Add elevation contours
        ax = plot_terrain(ax, ext)
        
        ## Add control ensemble time-averaged IVT
        ds = self.load_dataset(subset=False)
        
        Q = ax.quiver(ds.lon, ds.lat, ds.uIVT, ds.vIVT, transform=datacrs, 
                      color='k', regrid_shape=15,
                      angles='xy', scale_units='xy', scale=250, units='xy', **self.kw_quiver)
        
        # quiver key
        qk = ax.quiverkey(Q, 0.5, -0.15, 500, '500 kg m$^{-1}$ s$^{-1}$', labelpos='E',
                          coordinates='axes', fontproperties={'size': 5.0})
        
        
        ## Add 7-D QPF
        self.prec = load_prec_QPF_dataset(self.forecast, self.model_init_date)
        cmap, norm, bnds = cw3e.cmap('brian_qpf')
        self.qpflevs = bnds
        self.qpf = ax.contourf(self.prec.lon, self.prec.lat, self.prec.values,
                         cmap=cmap, norm=norm, levels=self.qpflevs, alpha=0.8, transform=datacrs)
        

        # Add map features (continents and country borders)
        ax.add_feature(cfeature.LAND, facecolor='0.9')      
        ax.add_feature(cfeature.BORDERS, edgecolor='0.4', linewidth=0.4)
        ax.add_feature(cfeature.STATES, edgecolor='0.2', linewidth=0.2)
        # ax.add_feature(cfeature.OCEAN, edgecolor='0.4', facecolor='lightskyblue', linewidth=0.2)

        # apply tick parameters   
        gl = ax.gridlines(crs=datacrs, draw_labels=True, **self.kw_grid)
        gl.top_labels = False
        gl.left_labels = False
        gl.right_labels = False
        gl.bottom_labels = True
        gl.xlocator = mticker.FixedLocator(dx)
        gl.ylocator = mticker.FixedLocator(dy)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = self.kw_ticklabels
        gl.ylabel_style = self.kw_ticklabels

        gl.xlines = True
        gl.ylines = True
        
        ax.set_yticks(dy, crs=ccrs.PlateCarree()) # add external yticks
        plt.yticks(color='w', size=1) # hack: make the ytick labels white so they don't show up
         
        ## plot point locations
        for i, (x, y) in enumerate(zip(self.lons, self.lats)):
            ax.plot(x, y, 'ko', markersize=2, transform=datacrs)

        ## plot labels
        if self.loc == 'US-west':
            label = 'Forecasts support FIRO/CA-AR Program Intended for research purposes only'
            ax.set_title(textwrap.fill(label, 36), loc='right', fontsize=5)
        else:
            label = 'Forecasts support NSF Coastlines and People Program: #2052972  Intended for research purposes only'
            ax.set_title(textwrap.fill(label, 61), loc='right', fontsize=5)
        
        return ax
    
    def plot_duration(self, ax):
    
        duration_ds = self.calc_ivt_duration()
        colors = ['tab:blue', 'green', 'tab:red', 'black']
        thresholds = [250, 500, 750, 1000]
        y = duration_ds.lat.values
        custom_lines = [] # for the custom legend
        legend_txt = []
        
        for i, thres in enumerate(thresholds):
            tmp = duration_ds.sel(duration=thres)
            custom_lines.append(Line2D([0], [0], color=colors[i], lw=0.6))
            legend_txt.append('$\geq$ {0}'.format(thres))
            # add ensemble mean as thicker line
            x = tmp.IVT.mean('ensemble').values 
            ax.plot(x, y, color=colors[i], linewidth=1)

            for j, ens in enumerate(range(len(duration_ds.ensemble))):
                x = tmp.sel(ensemble=ens).IVT.values
                ax.plot(x, y, color=colors[i], linewidth=0.25, alpha=0.4)

                ## xtick and ytick labels, locations
                ax.set_xlim(0, 144)
                ax.set_ylim(min(y), max(y))
                ax.set_xticks(np.arange(0, 168, 24))
                plt.xticks(**self.kw_ticklabels)
                plt.yticks(color='w', size=1)

                # apply gridlines and minor ticks
                ax.minorticks_on()
                ax.grid(visible=None, which='both', axis='y', **self.kw_grid)
                ax.grid(visible=None, which='major', axis='x', **self.kw_grid)
                ax.tick_params(axis='x', which='minor', bottom=True)
                ax.tick_params(axis='x', which='major', **self.kw_ticks)
                ax.tick_params(axis='y', which='major', direction='out', **self.kw_ticks)

        ## add legend
        ax.legend(custom_lines, legend_txt, loc='upper right', fontsize=5)

        ax.set_xlabel('duration (hours)', fontsize=5)
        ax.set_title('# hours with IVT $\geq$ threshold', loc='left', fontsize=5)
        
        return ax
    
    def create_figure(self):
        fname = '../figs/{0}/landfall-vector_{1}_{2}_{3}_{4}'.format(self.loc, self.ptloc, self.ensemble_name, self.threshold, self.date_string)
        fmt = 'png'
        
        if self.orientation == 'latitude':
            fig = plt.figure(figsize=(8, 3))
            fig.dpi = 300
            nrows = 2
            ncols = 3
            ## Use gridspec to set up a plot with a series of subplots that is
            ## n-rows by n-columns
            gs = GridSpec(nrows, ncols, height_ratios=[1, 0.05], width_ratios = [2, 0.5, 0.75], wspace=0.05, hspace=0.4)
            ## use gs[rows index, columns index] to access grids         

            ## Add probability plot         
            ax = fig.add_subplot(gs[0, 0])
            self.plot_vector_landfall_latitude(ax)

            ## Add color bar
            cbax = plt.subplot(gs[1,0]) # colorbar axis
            cb = Colorbar(ax = cbax, mappable = self.cf, orientation = 'horizontal', ticklocation = 'bottom', ticks=self.cflevs[::2])
            cb.ax.set_xticklabels(["{:.0%}".format(i) for i in cb.get_ticks()], **self.kw_ticklabels)  # horizontally oriented colorbar
            
            ## Add duration plot
            ax = fig.add_subplot(gs[0, 1])
            self.plot_duration(ax)

            # Set up projection information for map
            mapcrs = ccrs.PlateCarree()
            datacrs = ccrs.PlateCarree()
            ## Add map
            ax = fig.add_subplot(gs[0, 2], projection=mapcrs)
            self.plot_map(ax, mapcrs, datacrs)
            
            ## Add color bar for QPF
            cbax = plt.subplot(gs[1,-1]) # colorbar axis
            cb = Colorbar(ax = cbax, mappable = self.qpf, orientation = 'horizontal', ticklocation = 'bottom', ticks=self.qpflevs)
            cb.ax.set_xticklabels(["{0}".format(i) for i in cb.get_ticks()], **self.kw_ticklabels)  # horizontally oriented colorbar
            cb.set_label('QPF (in.)', fontsize=6)

            # ## Add CW3E logo
            # ax = fig.add_subplot(gs[1, 1:])
            # ax = plot_cw3e_logo(ax)
            
        
        elif self.orientation == 'longitude':
            fig = plt.figure(figsize=(8, 10))
            fig.dpi = 300
            nrows = 4
            ncols = 1
            ## Use gridspec to set up a plot with a series of subplots that is
            ## n-rows by n-columns
            gs = GridSpec(nrows, ncols, height_ratios=[0.7, 0.5, 1, 0.05], width_ratios = [1], hspace=0.1)
            ## use gs[rows index, columns index] to access grids         

            # Set up projection information for map
            mapcrs = ccrs.PlateCarree()
            datacrs = ccrs.PlateCarree()
            ## Add map
            ax = fig.add_subplot(gs[0, 0], projection=mapcrs)
            self.plot_map(ax, mapcrs, datacrs)
            
            ## Add duration plot
            ax = fig.add_subplot(gs[1, 0])
            self.plot_duration(ax)
            
            ## Add probability plot         
            ax = fig.add_subplot(gs[2, 0])
            self.plot_vector_landfall_longitude(ax)

            ## Add color bar
            cbax = plt.subplot(gs[3,0]) # colorbar axis
            cb = Colorbar(ax = cbax, mappable = self.cf, orientation = 'horizontal', ticklocation = 'bottom', ticks=self.cflevs[::2])
            cb.ax.set_xticklabels(["{:.0%}".format(i) for i in cb.get_ticks()], **self.kw_ticklabels)  # horizontally oriented colorbar

            ## Add CW3E logo
            ax = fig.add_subplot(gs[1, 1])
            ax = plot_cw3e_logo(ax)
        
        fig.savefig('%s.%s' %(fname, fmt), bbox_inches='tight', dpi=fig.dpi)
        plt.show()