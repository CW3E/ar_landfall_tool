#!/usr/bin/python3
"""
Filename:    ar_landfall_tool_contour.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Functions for CW3E AR Landfall Tool (contour)
"""
# Standard Python modules
import os, sys
import glob
import shutil
import numpy as np
import pandas as pd
import xarray as xr
import datetime
import re
import textwrap
from PIL import Image

# matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm, colors as clr
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec
from matplotlib.colorbar import Colorbar # different way to handle colorbar
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
from matplotlib import dates as mdates
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# cartopy
import cartopy.crs as ccrs
from cartopy.mpl.geoaxes import GeoAxes
import cartopy.feature as cfeature

# other
import cmocean.cm as cmo

# import personal modules
import cw3ecmaps as cw3e
from cw3e_tools import ivt_colors, plot_terrain, plot_cw3e_logo, get_every_other_vector, myround, load_prec_QPF_dataset


class landfall_tool_contour:
    '''
    Returns a .png file with Cordeira AR Landfall Tool (contour) figure with input locations, chosen Forecast product, and IVT threshold
    
    Parameters
    ----------
    ptloc : str
        name of the .txt file with latitude and longitude of locations for analysis
        this file should have no header, with latitude, then longitude, separated by a space
    forecast : str
        name of the forecast product - options include GEFS, ECMWF, ECMWF-GEFS, and West-WRF
    threshold : int
        threshold for IVT probabilty in kg m-1 s-1 - options include 150, 250, 500, 750
  
    Returns
    -------
    fig : figure
        png file of the figure
    
    '''
    
    def __init__(self, loc, ptloc, forecast='GEFS', threshold=250, orientation='latitude'):
        path_to_data = '/data/downloaded/SCRATCH/cw3eit_scratch/'
        self.forecast = forecast
        if forecast == 'GEFS':
            fpath = path_to_data + 'GEFS/FullFiles/'
            self.ensemble_name = 'GEFS'
        elif forecast == 'ECMWF':
            fpath = path_to_data + 'ECMWF/archive/' 
            self.ensemble_name = 'ECMWF'
        elif forecast == 'W-WRF':
            fpath = '/data/downloaded/WWRF-NRT/2023-2024/Ensemble_IVT/' 
            self.ensemble_name = 'West-WRF'
        elif forecast == 'ECMWF-GEFS':
            fpath = path_to_data + 'ECMWF/archive/' 
            self.ensemble_name = 'ECMWF-GEFS'
        else:
            print('Forecast product not available! Please choose either GEFS, ECMWF, ECMWF-GEFS, or W-WRF.')
                  
        ## find the most recent file in the currect directory
        list_of_files = glob.glob(fpath+'*.nc')
        self.fname = max(list_of_files, key=os.path.getctime)
        # pull the initialization date from the filename
        regex = re.compile(r'\d+')
        date_string = regex.findall(self.fname)
        if self.forecast == 'W-WRF':
            self.date_string = date_string[2]
        else:
            self.date_string = date_string[1]
        self.model_init_date = datetime.datetime.strptime(self.date_string, '%Y%m%d%H')
        
        # pull the name of the matching GEFS for ECMWF-GEFS
        self.fname2 = path_to_data + 'GEFS/FullFiles/IVT_Full_{0}.nc'.format(self.date_string)

        ## read text file with points
        self.loc = loc
        self.ptloc = ptloc
        textpts_fname = '../data/{0}/latlon_{1}.txt'.format(self.loc, self.ptloc)
        df = pd.read_csv(textpts_fname, header=None, sep=' ', names=['latitude', 'longitude'], engine='python')
        df['longitude'] = df['longitude']*-1
        self.df = df
        self.lons = self.df['longitude'].values
        self.lats = self.df['latitude'].values
        self.threshold = threshold
        
        ## info for plot orientation
        self.orientation = orientation
                  
        ## format dicts for plots
        if self.orientation == 'longitude':
            self.fontsize = 10
        else:
            self.fontsize = 12
            
        self.kw_ticklabels = {'size': self.fontsize-2, 'color': 'dimgray', 'weight': 'light'}
        self.kw_grid = {'linewidth': .5, 'color': 'k', 'linestyle': '--', 'alpha': 0.1}
        self.kw_ticks = {'length': 4, 'width': 0.5, 'pad': 2, 'color': 'black',
                         'labelsize': self.fontsize-2, 'labelcolor': 'dimgray'}
        self.IVT_units = 'kg m$^{-1}$ s$^{-1}$'
        
        
        
        
        if self.loc == 'US-west_old':
            self.grant_info = 'FIRO/CA-AR Program'
        elif self.loc == 'US-west':
            self.grant_info = 'FIRO/CA-AR Program and NSF #2052972'
        else:
            self.grant_info = 'NSF #2052972'
            
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
        
    def get_shared_axis_map_ticks(self):
        '''
        Returns
        -------
        ext : list
            list with [lonmin, lonmax, latmin, latmax]
        dx : array
            array with tick locations in the x-direction
        dy : array
            array with tick locations in the y-direction
        '''
        if self.orientation == 'latitude':
            # this extends the domain of the plot 2 degrees in the longitude direction
            londx=2. 
            lonmin = self.lons.min()-londx
            lonmax = self.lons.max()+londx
            # this extends the domain 0.25 degrees in the latitude direction
            # to show the entire point location
            latdy = 0.25
            latmin = self.lats.min()-latdy
            latmax = self.lats.max()+latdy

        elif self.orientation == 'longitude':
            # this extends the domain 0.25 degrees in the longitude direction
            # to show the entire point location
            londx = 0.25
            lonmin = self.lons.min()-londx
            lonmax = self.lons.max()+londx
            # this extends the domain of the plot 2 degrees in the latitude direction
            latdy = 2.
            latmin = self.lats.min()-latdy
            latmax = self.lats.max()+latdy

        self.ext = [lonmin, lonmax, latmin, latmax] # extent of map plot

        ## set nice intervals for dx based on diff(lonmax, lonmin)
        if (lonmax - lonmin) < 14:
            dux = 2
        else: 
            dux = 5

        ## set nice intervals for dy based on diff(latmax, latmin)
        if (latmax - latmin) < 14:
            duy = 2
        else:
            duy = 5

        self.dx = np.arange(myround(lonmin+londx, dux),(myround((lonmax-londx)+dux, dux)),dux)
        self.dy = np.arange(myround(latmin+latdy, duy),(myround((latmax-latdy)+duy, duy)),duy)
        
    def calc_ivt_probability(self, fname, forecast):
        ## load the forecast data
        ds = xr.open_dataset(fname)
        ds = ds.assign_coords({"lon": (((ds.lon + 180) % 360) - 180)}) # Convert DataArray longitude coordinates from 0-359 to -180-179
        ds = ds.sel(lon=slice(-180., -1)) # keep only western hemisphere
        if forecast == 'ECMWF':
            ds = ds.rename({'forecast_time': 'forecast_hour'}) # need to rename this to match GEFS
            ds = ds.assign_coords({"forecast_hour": (ds.forecast_hour*3)}) # convert forecast time to forecast hour
        elif forecast == 'W-WRF':
            ds = ds.rename({'ensembles': 'ensemble'}) # need to rename this to match GEFS/ECMWF
        
        self.ndays = str(int(ds.forecast_hour.values.max()/24.)) # get max number of forecast days
        
        # subset ds to the select points
        x = xr.DataArray(self.lons, dims=['location'])
        y = xr.DataArray(self.lats, dims=['location'])
        ds = ds.sel(lon=x, lat=y, method='nearest')

        ## calculate probability IVT >= threshold
        data_size = ds.IVT.ensemble.shape

        # sum the number of ensembles where IVT exceeds threshold
        probability = (ds.IVT.where(ds.IVT >= self.threshold).count(dim='ensemble')) / data_size
        
        return probability
            
    def plot_probability_latitude(self, ax):
        
        if self.forecast == 'ECMWF-GEFS':
            ECMWF = self.calc_ivt_probability(self.fname, forecast='ECMWF') # run the calculation for ECMWF
            GFS = self.calc_ivt_probability(self.fname2, forecast='GEFS') # run the calculation for GFS
            self.probability = ECMWF - GFS # find the difference
            cmap = plt.cm.get_cmap('BrBG')
            self.cflevs = np.arange(-0.55, 0.6, 0.05)
            norm = mcolors.BoundaryNorm(self.cflevs, cmap.N)
        else:
            self.probability = self.calc_ivt_probability(self.fname, forecast=self.forecast) # run the calculation for model
            cmap, norm, bnds = cw3e.cmap('ivt_probability')
            self.cflevs = bnds
        
        # self.calc_ivt_probability() # run the calculation
        
        # Contour Filled
        x = self.probability.forecast_hour.values / 24 # convert forecast hour to forecast day
        y = self.probability.lat.values
        data = np.flipud(np.rot90(self.probability.values)) # rotate data 90 degrees and flip up down

        self.cf = ax.pcolormesh(x, y, data, cmap=cmap, norm=norm, rasterized=True)
        ax.invert_xaxis() # invert x-axis so that time reads from right to left

        # apply ytick parameters (latitude labels)
        ax.yaxis.set_major_locator(mticker.FixedLocator(self.dy))
        ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)
        for tick in ax.get_yticklabels():
            tick.set_fontweight('light')

        # apply xtick parameters
        positions = np.arange(0, 17, 1)
        ax.xaxis.set_major_locator(mticker.FixedLocator(positions))
        ax.xaxis.set_major_formatter(mticker.FixedFormatter(self.xtck_lbl))
        for tick in ax.get_xticklabels():
            tick.set_fontweight('light')
            
        # labels are days since forecast initialization
        for i, x in enumerate(positions):
            ax.annotate(u"{:0.0f}".format(x), # this is the text
                       (x,y.min()), # these are the coordinates to position the label
                        textcoords="offset points", # how to position the text
                        xytext=(0,0), # distance from text to points (x,y)
                        ha='center', # horizontal alignment can be left, right or center
                        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="k", lw=0.5, alpha=0.8),
                        xycoords='data',
                        zorder=200,
                        fontsize=self.fontsize)
            
        ax.set_xlim(16.125, 0)

        # apply gridlines
        ax.minorticks_on()
        ax.grid(visible=None, which='both', axis='y', **self.kw_grid)
        ax.grid(visible=None, which='major', axis='x', **self.kw_grid)
        ax.tick_params(axis='x', which='minor', bottom=False)
        ax.tick_params(axis='x', which='major', **self.kw_ticks)
        ax.tick_params(axis='y', which='major', direction='out', **self.kw_ticks)

        ## labels and subtitles
        ax.set_ylabel("Latitude along West Coast", fontsize=self.fontsize)
        ax.set_xlabel(self.xlbl, fontsize=self.fontsize)
        ax.set_title('{0}-d {1} Probability of IVT $\geq$ {2} {3}'.format(self.ndays, self.ensemble_name, self.threshold, self.IVT_units), loc='left', fontsize=self.fontsize)
        
        return ax
    
    def plot_probability_longitude(self, ax):
        
        if self.forecast == 'ECMWF-GEFS':
            ECMWF = self.calc_ivt_probability(self.fname, forecast='ECMWF') # run the calculation for ECMWF
            GFS = self.calc_ivt_probability(self.fname2, forecast='GEFS') # run the calculation for GFS
            self.probability = ECMWF - GFS # find the difference
            cmap = plt.cm.get_cmap('BrBG')
            self.cflevs = np.arange(-0.55, 0.6, 0.05)
            norm = mcolors.BoundaryNorm(self.cflevs, cmap.N)
        else:
            self.probability = self.calc_ivt_probability(self.fname, forecast=self.forecast) # run the calculation for model
            cmap, norm, bnds = cw3e.cmap('ivt_probability')
            self.cflevs = bnds
        
        # Contour Filled
        y = self.probability.forecast_hour / 24 # convert forecast hour to forecast day
        x = self.probability.lon
        data = self.probability.values
        
        self.cf = ax.pcolormesh(x, y, data, cmap=cmap, norm=norm, rasterized=True)
        ax.invert_xaxis() # invert x-axis so that time reads from right to left

        # apply xtick parameters (longitude labels)
        ax.xaxis.set_major_locator(mticker.FixedLocator(self.dx))
        ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
        for tick in ax.get_xticklabels():
            tick.set_fontweight('light')
        
        # apply ytick parameters
        positions = np.arange(0, 17, 1)
        ax.yaxis.set_major_locator(mticker.FixedLocator(positions))
        ax.yaxis.set_major_formatter(mticker.FixedFormatter(self.xtck_lbl))
        for tick in ax.get_yticklabels():
            tick.set_fontweight('light')
        
        # labels are days since forecast initialization
        for i, y in enumerate(positions):
            ax.annotate(u"{:0.0f}".format(y), # this is the text
                       (x.min(),y), # these are the coordinates to position the label
                        textcoords="offset points", # how to position the text
                        xytext=(-3,-3), # distance from text to points (x,y)
                        ha='left', # horizontal alignment can be left, right or center
                        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="k", lw=0.5, alpha=0.8),
                        xycoords='data',
                        zorder=200,
                        fontsize=self.fontsize)
        ax.set_ylim(0, 16.125)
        
        # apply gridlines
        ax.minorticks_on()
        ax.grid(visible=None, which='both', axis='x', **self.kw_grid)
        ax.grid(visible=None, which='major', axis='y', **self.kw_grid)
        ax.tick_params(axis='y', which='minor', left=False)
        ax.tick_params(axis='y', which='major', **self.kw_ticks)
        ax.tick_params(axis='x', which='major', direction='out', **self.kw_ticks)
        
        ## labels and subtitles
        ax.set_ylabel(self.xlbl, fontsize=self.fontsize)
        # ax.set_xlabel('Probability of IVT $\geq$ {0} {1}'.format(self.threshold, self.IVT_units), fontsize=self.fontsize)
        
        plt.gca().invert_yaxis()
        plt.gca().invert_xaxis()
        
        return ax
    
                  
    def plot_map(self, ax, mapcrs, datacrs):
        ## Set up extent         
        ax.set_extent(self.ext, crs=datacrs)
        
        ## Add elevation contours
        ax = plot_terrain(ax, self.ext)

        # Add map features (continents and country borders)
        # ax.add_feature(cfeature.LAND, facecolor='0.9')      
        ax.add_feature(cfeature.BORDERS, edgecolor='0.4', linewidth=0.4)
        ax.add_feature(cfeature.STATES, edgecolor='0.2', linewidth=0.2)
        ax.add_feature(cfeature.OCEAN, edgecolor='0.4', facecolor='lightskyblue', linewidth=0.2)

        # add gridlines
        gl = ax.gridlines(crs=datacrs, draw_labels=True, **self.kw_grid)
        gl.top_labels = False
        gl.right_labels = False
        if self.orientation == 'latitude':
            gl.left_labels = False
            gl.bottom_labels = True
            mk_size = 3
        else:
            gl.left_labels = True
            gl.bottom_labels = False
            mk_size = 5
        gl.xlocator = mticker.FixedLocator(self.dx)
        gl.ylocator = mticker.FixedLocator(self.dy)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = self.kw_ticklabels
        gl.ylabel_style = self.kw_ticklabels

        gl.xlines = True
        gl.ylines = True
        
        # apply tick parameters
        ax.set_xticks(self.dx, crs=datacrs)
        ax.set_yticks(self.dy, crs=datacrs)
        plt.yticks(color='w', size=1) # hack: make the ytick labels white so they don't show up
        plt.xticks(color='w', size=1) # hack: make the ytick labels white so they don't show up
         
        ## plot point locations
        for i, (x, y) in enumerate(zip(self.lons, self.lats)):
            ax.plot(x, y, 'ko', markersize=mk_size, transform=datacrs)
            
        if self.orientation == 'longitude':
            txt = '{0}-d {1} Probability of IVT $\geq$ {2} {3}'.format(self.ndays, self.ensemble_name, self.threshold, self.IVT_units) 
            ax.set_title(txt, loc='left', fontsize=self.fontsize)
        
        ax.set_title(self.title, loc='right', fontsize=self.fontsize)
        ax.set_extent(self.ext, crs=datacrs)
        ax.set_aspect('auto')
        return ax
    
    def create_figure(self):
        fname1 = '../figs/{0}/{1}_LandfallTool_{2}_{3}_current'.format(self.loc, self.forecast, self.threshold, self.ptloc)
        fname2 = '../figs/{0}/{1}_LandfallTool_{2}_{3}_{4}'.format(self.loc, self.forecast, self.threshold, self.ptloc, self.date_string)
        fmt = 'png'
        
        # get tick and label information
        self.get_date_information()
        self.get_shared_axis_map_ticks()
        
        if self.orientation == 'latitude':
            fig = plt.figure(figsize=(13., 6))
            fig.dpi = 300
            nrows = 4
            ncols = 2
            ## Use gridspec to set up a plot with a series of subplots that is
            ## n-rows by n-columns
            gs = GridSpec(nrows, ncols, height_ratios=[1, 0.05, 0.05, 0.05], width_ratios = [2.25, 0.75], wspace=0.05, hspace=0.2)
            ## use gs[rows index, columns index] to access grids         

            ## Add probability plot         
            ax = fig.add_subplot(gs[0, 0])
            self.plot_probability_latitude(ax)

            ## Add color bar
            cbax = plt.subplot(gs[2,0]) # colorbar axis
            cb = Colorbar(ax = cbax, mappable = self.cf, orientation = 'horizontal', ticklocation = 'bottom', ticks=self.cflevs[1:-1:2])
            cb.ax.set_xticklabels(["{:.0%}".format(i) for i in cb.get_ticks()], **self.kw_ticklabels)  # horizontally oriented colorbar

            # Set up projection information for map
            mapcrs = ccrs.PlateCarree()
            datacrs = ccrs.PlateCarree()
            ## Add map
            ax = fig.add_subplot(gs[0, 1], projection=mapcrs)
            self.plot_map(ax, mapcrs, datacrs)
            
            ## labels and subtitles
            ax = fig.add_subplot(gs[3, 0])
            ax.axis('off')
            title = 'Forecasts support {0} | Intended for research purposes only'.format(self.grant_info)
            ax.annotate(title, # this is the text
                       (0, 0.1), # these are the coordinates to position the label
                        textcoords="offset points", # how to position the text
                        xytext=(0,0), # distance from text to points (x,y)
                        ha='left', # horizontal alignment can be left, right or center
                        fontsize=self.fontsize-2)

            ## Add CW3E logo
            ax = fig.add_subplot(gs[1:, 1])
            ax = plot_cw3e_logo(ax, orientation='horizontal')
            
        
        elif self.orientation == 'longitude':
            fig = plt.figure(figsize=(9, 12))
            fig.dpi = 300
            nrows = 5
            ncols = 2
            ## Use gridspec to set up a plot with a series of subplots that is
            ## n-rows by n-columns
            gs = GridSpec(nrows, ncols, height_ratios=[0.7, 1, 0.02, 0.05, 0.05], width_ratios = [0.75, 0.25], hspace=0.06, wspace=0.02)
            ## use gs[rows index, columns index] to access grids
            ## Add probability plot         
            ax = fig.add_subplot(gs[1, :])
            self.plot_probability_longitude(ax)

            ## Add color bar
            cbax = plt.subplot(gs[3,0]) # colorbar axis
            cb = Colorbar(ax = cbax, mappable = self.cf, orientation = 'horizontal', ticklocation = 'bottom', ticks=self.cflevs[1:-1:2])
            cb.ax.set_xticklabels(["{:.0%}".format(i) for i in cb.get_ticks()], **self.kw_ticklabels)  # horizontally oriented colorbar
       
            # Set up projection information for map
            mapcrs = ccrs.PlateCarree()
            datacrs = ccrs.PlateCarree()
            ## Add map
            ax = fig.add_subplot(gs[0, :], projection=mapcrs)
            self.plot_map(ax, mapcrs, datacrs)

                 
            ## labels and subtitles
            ax = fig.add_subplot(gs[4, 0])
            ax.axis('off')
            title = 'Forecasts support {0} | Intended for research purposes only'.format(self.grant_info)
            ax.annotate(title, # this is the text
                       (0, 0.), # these are the coordinates to position the label
                        textcoords="offset points", # how to position the text
                        xytext=(0,0), # distance from text to points (x,y)
                        ha='left', # horizontal alignment can be left, right or center
                        fontsize=self.fontsize-2)

            ## Add CW3E logo
            ax = fig.add_subplot(gs[3:, 1])
            ax = plot_cw3e_logo(ax, orientation='horizontal')
            
        
        fig.savefig('%s.%s' %(fname1, fmt), bbox_inches='tight', dpi=fig.dpi) # save generic "current"
        fig.savefig('%s.%s' %(fname2, fmt), bbox_inches='tight', dpi=fig.dpi) # save with date/time
        # close figure
        plt.close(plt.gcf())
        
   