#!/usr/bin/python3
"""
Filename:    ar_landfall_tool_IVT_mag.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Functions for CW3E AR Landfall Tool (IVT Magnitude)
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
from cw3e_tools import ivt_colors, plot_terrain, plot_cw3e_logo, get_every_other_vector, myround


class landfall_tool_IVT_magnitude:
    '''
    Returns a .png file with Cordeira AR Landfall Tool (IVT magnitude) figure with input locations, chosen Forecast product, for Control or Ensemble Mean   
    Parameters
    ----------
    ptloc : str
        name of the .txt file with latitude and longitude of locations for analysis
        this file should have no header, with latitude, then longitude, separated by a space
    forecast : str
        name of the forecast product - options include GEFS, ECMWF, ECMWF-GEFS, and West-WRF
    mag_type : str
        choose the magnitude type to calculate - ensemble or control (ensemble = 0)
  
    Returns
    -------
    fig : figure
        png file of the figure
    
    '''
    
    def __init__(self, ds_pt, loc, ptloc, forecast='GEFS', mag_type='control', orientation='latitude'):
        self.path_to_out = '/home/cw3eit/ARPortal/gefs/scripts/ar_landfall_tool/figs/'
        self.date_string = ds_pt.model_init_date.strftime('%Y%m%d%H')
        self.model_init_date = datetime.datetime.strptime(self.date_string, '%Y%m%d%H')
        self.loc = loc
        self.ptloc = ptloc
        self.lons = ds_pt.lon.values
        self.lats = ds_pt.lat.values
        self.orientation = orientation
        self.forecast = forecast
        self.mag_type = mag_type
        self.magnitude = ds_pt[self.mag_type]
        if self.mag_type == 'ensemble_mean':
            self.name = 'Ensemble Mean'
        else:
            self.name = 'Control'
                  
        ## format dicts for plots
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
            
    def plot_magnitude_latitude(self, ax):
        
        if self.forecast == 'ECMWF-GEFS':
            cmap = plt.cm.get_cmap('BrBG')
            self.cflevs = np.arange(-275, 300, 25)
            norm = mcolors.BoundaryNorm(self.cflevs, cmap.N)
        else:
            cmap, norm, bnds = cw3e.cmap('ivt_magnitude')
            self.cflevs = bnds
        
        # Contour Filled
        x = self.magnitude.forecast_hour.values / 24 # convert forecast hour to forecast day
        y = self.magnitude.lat.values
        data = np.flipud(np.rot90(self.magnitude.values)) # rotate data 90 degrees and flip up down
        
        
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
        ax.set_title('16-d {0} {1} IVT ({2})'.format(self.forecast, self.name, self.IVT_units), loc='left', fontsize=self.fontsize)
        
        return ax
    
    def plot_magnitude_longitude(self, ax):
        
        if self.forecast == 'ECMWF-GEFS':
            cmap = plt.cm.get_cmap('BrBG')
            self.cflevs = np.arange(-275, 300, 25)
            norm = mcolors.BoundaryNorm(self.cflevs, cmap.N)
        else:
            cmap, norm, bnds = cw3e.cmap('ivt_magnitude')
            self.cflevs = bnds
        
        # Contour Filled
        y = self.magnitude.forecast_hour / 24 # convert forecast hour to forecast day
        x = self.magnitude.lon
        data = self.magnitude.values
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
            txt = '16-d {0} {1} IVT ({2})'.format(self.forecast, self.name, self.IVT_units) 
            ax.set_title(txt, loc='left', fontsize=self.fontsize)
        
        ax.set_title(self.title, loc='right', fontsize=self.fontsize)
        ax.set_extent(self.ext, crs=datacrs)
        ax.set_aspect('auto')
        return ax
    
    def create_figure(self):
        # GEFS_LandfallTool_Control_coast_current.png
        fname1 = self.path_to_out+'{0}/{1}_LandfallTool_{2}_{3}_current'.format(self.loc, self.forecast, self.mag_type, self.ptloc)
        fname2 = self.path_to_out+'{0}/{1}_LandfallTool_{2}_{3}_{4}'.format(self.loc, self.forecast, self.mag_type, self.ptloc, self.date_string)
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
            self.plot_magnitude_latitude(ax)

            ## Add color bar
            cbax = plt.subplot(gs[2,0]) # colorbar axis
            cb = Colorbar(ax = cbax, mappable = self.cf, orientation = 'horizontal', ticklocation = 'bottom', ticks=self.cflevs[1:-2:2])
            cb.ax.set_xticklabels(["{0}".format(i) for i in cb.get_ticks()], **self.kw_ticklabels)  # horizontally oriented colorbar

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
            self.plot_magnitude_longitude(ax)

            ## Add color bar
            cbax = plt.subplot(gs[3,0]) # colorbar axis
            cb = Colorbar(ax = cbax, mappable = self.cf, orientation = 'horizontal', ticklocation = 'bottom', ticks=self.cflevs[1:-2:2])
            cb.ax.set_xticklabels(["{0}".format(i) for i in cb.get_ticks()], **self.kw_ticklabels)  # horizontally oriented colorbar
       
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