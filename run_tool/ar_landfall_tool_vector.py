#!/usr/bin/python3
"""
Filename:    ar_landfall_tool_vector.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Functions for CW3E AR Landfall Tool (vector)
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
        if forecast == 'GEFS':
            fpath = path_to_data + 'GEFS/FullFiles/'
            self.ensemble_name = 'GEFS'
        elif forecast == 'ECMWF':
            fpath = path_to_data + 'ECMWF/archive/' # will need to adjust when operational
            self.ensemble_name = 'ECMWF'
        else:
            print('Forecast product not available! Please choose either GEFS or ECMWF.')
                  
        ## find the most recent file in the currect directory
        list_of_files = glob.glob(fpath+'*.nc')
        self.fname = max(list_of_files, key=os.path.getctime)
        # pull the initialization date from the filename
        regex = re.compile(r'\d+')
        date_string = regex.findall(self.fname)
        self.date_string = date_string[1]
        self.model_init_date = datetime.datetime.strptime(self.date_string, '%Y%m%d%H')

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
                  
        ## format dicts for plots
        self.fontsize = 12
        self.kw_ticklabels = {'size': self.fontsize-2, 'color': 'dimgray', 'weight': 'light'}
        self.kw_grid = {'linewidth': .5, 'color': 'k', 'linestyle': '--', 'alpha': 0.1}
        self.kw_ticks = {'length': 4, 'width': 0.5, 'pad': 2, 'color': 'black',
                         'labelsize': self.fontsize-2, 'labelcolor': 'dimgray'}
        self.kw_quiver = {'headlength': 6, 'headaxislength': 4.5, 'headwidth': 4.5}
        self.IVT_units = 'kg m$^{-1}$ s$^{-1}$'
        
        ## info for plot orientation
        self.orientation = orientation
        
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
        date_lst = pd.date_range(self.model_init_date, periods=8, freq='1D')
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
        
    def load_dataset(self, subset=True):
        ## load the forecast data
        ds = xr.open_dataset(self.fname)
        ds = ds.assign_coords({"lon": (((ds.lon + 180) % 360) - 180)}) # Convert DataArray longitude coordinates from 0-359 to -180-179
        ds = ds.sel(lon=slice(-180., -1)) # keep only western hemisphere
        if self.forecast == 'ECMWF':
            ds = ds.rename({'forecast_time': 'forecast_hour'}) # need to rename this to match GEFS
            ds = ds.assign_coords({"forecast_hour": (ds.forecast_hour*3)}) # convert forecast time to forecast hour
        
        self.ndays = str(int(ds.forecast_hour.values.max()/24.)) # get max number of forecast days
        
        
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
        
        
    def plot_duration_cbar(self, cbax):
        # create custom colorbar for duration plot
        upper = 1000 # the upper limit for the colorbar
        lower = 250 # the lower limit for the colorbar
        N = 4 # the number of discrete intervals
        deltac = (upper-lower)/(2*(N-1))
        cmap, norm, cbar_tcks = cw3e.cmap('ivt_duration')
        norm = clr.Normalize() # this alters the state of the Normalize object
        duration_cbar = cm.ScalarMappable(norm=norm, cmap=cmap)
        duration_cbar.set_array([lower-deltac,upper+deltac])
        if self.orientation == 'latitude':
            cb = Colorbar(ax = cbax, mappable = duration_cbar, orientation = 'horizontal', ticklocation = 'bottom', ticks=[250, 500, 750, 1000])
            cb.set_label('IVT $\geq$ threshold', fontsize=self.fontsize)
            cb.ax.set_xticklabels(["{0}".format(i) for i in cb.get_ticks()], **self.kw_ticklabels)
        else:
            cb = Colorbar(ax = cbax, mappable = duration_cbar, orientation = 'vertical', ticklocation = 'right', ticks=[250, 500, 750, 1000])
            cb.set_label('IVT $\geq$ threshold', fontsize=self.fontsize)
            cb.ax.set_yticklabels(["{0}".format(i) for i in cb.get_ticks()], **self.kw_ticklabels)
        
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
        cs = ax.contour(x, y, ctrl, levels=clevs, colors='k', linewidths=1.0)
        
        ax.invert_xaxis() # invert x-axis so that time reads from right to left
                  
        # apply ytick parameters (latitude labels)
        ax.yaxis.set_major_locator(mticker.FixedLocator(self.dy))
        ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)
        for tick in ax.get_yticklabels():
            tick.set_fontweight('light')

        # apply xtick parameters
        positions = np.arange(0, 8, 1)
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
                        fontsize=self.fontsize-2)

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
        ax.set_title(self.title, loc='right', fontsize=self.fontsize)
        ax.set_title('(a) 7-d {0} Ens. Mean IVT'.format(self.ensemble_name), loc='left', fontsize=self.fontsize)
        
        return ax
    
    def plot_vector_landfall_longitude(self, ax):
        
        self.calc_ivt_probability() # run the calculation
        x = self.probability.lon
        y = self.probability.forecast_hour / 24 # convert forecast hour to forecast day
        data = get_every_other_vector(self.probability.values)
        uvec = get_every_other_vector(self.u.values)
        vvec = get_every_other_vector(self.v.values)
        ctrl = self.control.values
        
        ## Quiver
        cmap, norm, bnds = cw3e.cmap('ivt_vector')
        self.cflevs = bnds
        self.cf = ax.quiver(x, y, uvec, vvec, data, cmap=cmap, norm=norm,
                            capstyle='round', units='width', **self.kw_quiver)
        
        ## Contours
        ## add contour lines of control ensemble IVT every 250 kg m-1 s-1
        clevs = np.arange(250, 5000, 250)
        cs = ax.contour(x, y, ctrl, levels=clevs, colors='k', linewidths=1.0)
            
        ax.invert_xaxis() # invert x-axis so that time reads from right to left
                  
        # apply xtick parameters (longitude labels)
        ax.xaxis.set_major_locator(mticker.FixedLocator(self.dx))
        ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
        for tick in ax.get_xticklabels():
            tick.set_fontweight('light')
        
        # apply ytick parameters
        positions = np.arange(0, 8, 1)
        ax.yaxis.set_major_locator(mticker.FixedLocator(positions))
        ax.yaxis.set_major_formatter(mticker.FixedFormatter(self.xtck_lbl))
        for tick in ax.get_yticklabels():
            tick.set_fontweight('light')
        
        # labels are days since forecast initialization
        for i, y in enumerate(positions):
            ax.annotate(u"{:0.0f}".format(y), # this is the text
                       (x.min(),y), # these are the coordinates to position the label
                        textcoords="offset points", # how to position the text
                        xytext=(0,-3), # distance from text to points (x,y)
                        ha='center', # horizontal alignment can be left, right or center
                        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="k", lw=0.5, alpha=0.8),
                        xycoords='data',
                        zorder=200,
                        fontsize=self.fontsize-2)
            
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
        ax.set_extent(self.ext, crs=datacrs)
        
        ## Add elevation contours
        ax = plot_terrain(ax, self.ext)
        
        ## Add 7-D QPF
        self.prec = load_prec_QPF_dataset(self.forecast, self.model_init_date, self.date_string)
        cmap, norm, bnds = cw3e.cmap('brian_qpf')
        self.qpflevs = bnds
        self.qpf = ax.contourf(self.prec.lon, self.prec.lat, self.prec.values,
                         cmap=cmap, norm=norm, levels=self.qpflevs, alpha=0.8, transform=datacrs)
        
        ## Add control ensemble time-averaged IVT
        ds = self.load_dataset(subset=False)
        
        Q = ax.quiver(ds.lon, ds.lat, ds.uIVT, ds.vIVT, transform=datacrs, 
                      color='k', regrid_shape=15,
                      angles='xy', scale_units='xy', scale=250, units='xy')
        

        # Add map features (continents and country borders)
        # ax.add_feature(cfeature.LAND, facecolor='0.9')      
        ax.add_feature(cfeature.BORDERS, edgecolor='0.4', linewidth=0.4)
        ax.add_feature(cfeature.STATES, edgecolor='0.2', linewidth=0.2)
        # ax.add_feature(cfeature.OCEAN, edgecolor='0.4', facecolor='lightskyblue', linewidth=0.2)

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

        ## subtitles
        ax.set_title('(c)', loc='left', fontsize=self.fontsize)
    
            
        if self.orientation == 'longitude':
            txt = '7-d {0} Ensemble Mean IVT'.format(self.ensemble_name) 
            ax.set_title(txt, loc='left', fontsize=self.fontsize)
            ax.set_title(self.title, loc='right', fontsize=self.fontsize)
            qk_x = 0.82
            qk_y = 0.95
            
        else:
            qk_x = 0.45
            qk_y = -0.1
              
        # quiver key
        qk = ax.quiverkey(Q, qk_x, qk_y, 500, '500 kg m$^{-1}$ s$^{-1}$', labelpos='E',
                          coordinates='axes', fontproperties={'size': self.fontsize-2})
        
        t = qk.text.set_backgroundcolor('w')
            
        ax.set_extent(self.ext, crs=datacrs)
        ax.set_aspect('auto')
        return ax
    
    def plot_duration_latitude(self, ax):
    
        duration_ds = self.calc_ivt_duration()
        colors = [ivt_colors['250'], ivt_colors['500'], ivt_colors['750'], ivt_colors['1000']]
        thresholds = [250., 500., 750., 1000.]
        y = duration_ds.lat.values
        self.custom_lines = [] # for the custom legend
        self.legend_txt = []
        
        for i, thres in enumerate(thresholds):
            tmp = duration_ds.sel(duration=thres)
            self.custom_lines.append(Line2D([0], [0], color=colors[i], lw=0.8))
            self.legend_txt.append('$\geq$ {0}'.format(thres))
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
                
                # get all the labels of this axis
                labels = ax.get_xticklabels()
                # remove the first and the last labels
                labels[0] = labels[-1] = ""
                # set these new labels
                ax.set_xticklabels(labels)

                # apply gridlines and minor ticks
                ax.minorticks_on()
                ax.grid(visible=None, which='both', axis='y', **self.kw_grid)
                ax.grid(visible=None, which='major', axis='x', **self.kw_grid)
                ax.tick_params(axis='x', which='minor', bottom=True)
                ax.tick_params(axis='x', which='major', **self.kw_ticks)
                ax.tick_params(axis='y', which='major', direction='out', **self.kw_ticks)
                plt.yticks(color='w', size=1)

        txt = 'duration (hours)'
        ax.set_xlabel(textwrap.fill(txt, 25), fontsize=self.fontsize)
        
        ## subtitles
        ax.set_title('(b)', loc='left', fontsize=self.fontsize)
        
        return ax
    
    def plot_duration_longitude(self, ax):
    
        duration_ds = self.calc_ivt_duration()
        colors = [ivt_colors['250'], ivt_colors['500'], ivt_colors['750'], ivt_colors['1000']]
        thresholds = [250., 500., 750., 1000.]
        x = duration_ds.lon.values
        self.custom_lines = [] # for the custom legend
        self.legend_txt = []
        
        for i, thres in enumerate(thresholds):
            tmp = duration_ds.sel(duration=thres)
            self.custom_lines.append(Line2D([0], [0], color=colors[i], lw=0.6))
            self.legend_txt.append('$\geq$ {0}'.format(thres))
            # add ensemble mean as thicker line
            y = tmp.IVT.mean('ensemble').values 
            ax.plot(x, y, color=colors[i], linewidth=1)

            for j, ens in enumerate(range(len(duration_ds.ensemble))):
                y = tmp.sel(ensemble=ens).IVT.values
                ax.plot(x, y, color=colors[i], linewidth=0.25, alpha=0.4)

                ## xtick and ytick labels, locations
                ax.set_ylim(0, 144)
                ax.set_xlim(min(x), max(x))
                ax.set_yticks(np.arange(0, 168, 24))
                plt.yticks(**self.kw_ticklabels)
                
                # apply gridlines and minor ticks
                ax.minorticks_on()
                ax.grid(visible=None, which='both', axis='x', **self.kw_grid)
                ax.grid(visible=None, which='major', axis='y', **self.kw_grid)
                ax.tick_params(axis='y', which='minor', bottom=True)
                ax.tick_params(axis='y', which='major', **self.kw_ticks)
                ax.tick_params(axis='x', which='major', direction='out', **self.kw_ticks)
                plt.xticks(color='w', size=1)

        txt = 'duration (hours)'
        ax.set_ylabel(textwrap.fill(txt, 25), fontsize=self.fontsize)
        
        return ax
    
    def create_figure(self):
        fname1 = '../figs/{0}/{1}_LandfallTool_Vectors_{2}_{3}_current'.format(self.loc, self.forecast, self.threshold, self.ptloc)
        fname2 = '../figs/{0}/{1}_LandfallTool_Vectors_{2}_{3}_{4}'.format(self.loc, self.forecast, self.threshold, self.ptloc, self.date_string)
        fmt = 'png'
        
        # get tick and label information
        self.get_date_information()
        self.get_shared_axis_map_ticks()
        
        ## (a) (b) (c) label location info
        x=0.01
        y=0.973
        
        if self.orientation == 'latitude':
            fig = plt.figure(figsize=(12, 6))
            fig.dpi = 300
            nrows = 3
            ncols = 3
            ## Use gridspec to set up a plot with a series of subplots that is
            ## n-rows by n-columns
            gs = GridSpec(nrows, ncols, height_ratios=[1, 0.05, 0.05], width_ratios = [2, 0.5, 0.75], wspace=0.05, hspace=0.5)
            ## use gs[rows index, columns index] to access grids         

            ## Add probability plot         
            ax = fig.add_subplot(gs[0, 0])
            self.plot_vector_landfall_latitude(ax)


            ## Add color bar
            cbax = plt.subplot(gs[1,0]) # colorbar axis
            cb = Colorbar(ax = cbax, mappable = self.cf, orientation = 'horizontal', ticklocation = 'bottom', ticks=self.cflevs[::2])
            cb.ax.set_xticklabels(["{:.0%}".format(i) for i in cb.get_ticks()], **self.kw_ticklabels)  # horizontally oriented colorbar
            cb.set_label('Probability of IVT $\geq$ {0} {1}'.format(self.threshold, self.IVT_units), fontsize=self.fontsize)
            
            ## Add duration plot
            ax = fig.add_subplot(gs[0, 1])
            self.plot_duration_latitude(ax)
            
            ## add cmap for duration plot
            cbax = fig.add_subplot(gs[1, 1])
            # create custom colorbar for duration plot
            self.plot_duration_cbar(cbax)

            # Set up projection information for map
            mapcrs = ccrs.PlateCarree()
            datacrs = ccrs.PlateCarree()
            ## Add map
            ax = fig.add_subplot(gs[0, 2], projection=mapcrs)
            self.plot_map(ax, mapcrs, datacrs)
            
            ## Add color bar for QPF
            cbax = plt.subplot(gs[1,2]) # colorbar axis
            # cbax = inset_axes(ax, width="3%", height="45%", loc='lower left')
            cb = Colorbar(ax = cbax, mappable = self.qpf, orientation = 'horizontal', ticklocation = 'bottom', ticks=self.qpflevs)
            # horizontally oriented colorbar
            cb.ax.set_xticklabels(["{0}".format(i) for i in cb.get_ticks()], **self.kw_ticklabels)  
            cb.set_label('QPF (in.)', fontsize=self.fontsize)
            # cbax.set_title('QPF (in.)', fontsize=self.fontsize-2)

            ## Add CW3E logo
            in_ax = inset_axes(ax, width="40%", height="25%", loc='lower left')
            in_ax = plot_cw3e_logo(in_ax, orientation='vertical')
            
            ## Add grant information
            ax = fig.add_subplot(gs[2, :])
            ax.axis('off')
            title = 'Forecasts support {0} | Intended for research purposes only'.format(self.grant_info)
            ax.annotate(title, # this is the text
                       (0, 0.), # these are the coordinates to position the label
                        textcoords="offset points", # how to position the text
                        xytext=(0,0), # distance from text to points (x,y)
                        ha='left', # horizontal alignment can be left, right or center
                        fontsize=self.fontsize-2)
        
        elif self.orientation == 'longitude':
            fig = plt.figure(figsize=(10, 15))
            fig.dpi = 300
            nrows = 6
            ncols = 3
            ## Use gridspec to set up a plot with a series of subplots that is
            ## n-rows by n-columns
            gs = GridSpec(nrows, ncols, height_ratios=[0.7, 0.5, 1, 0.02, 0.05, 0.05], width_ratios = [0.75, 0.25, 0.03], hspace=0.06, wspace=0.01)
            ## use gs[rows index, columns index] to access grids         

            # Set up projection information for map
            mapcrs = ccrs.PlateCarree()
            datacrs = ccrs.PlateCarree()
            ## Add map
            ax = fig.add_subplot(gs[0, 0:2], projection=mapcrs)
            self.plot_map(ax, mapcrs, datacrs)
            ax.text(x, y, '(a)', ha='left', va='top', transform=ax.transAxes, fontsize=12., backgroundcolor='white', zorder=101)

            
            # add QPF colorbar
            cbax = plt.subplot(gs[0,2]) # colorbar axis
            # cbax = inset_axes(ax, width="30%", height="3%", loc='upper right')
            cb = Colorbar(ax = cbax, mappable = self.qpf, orientation = 'vertical', ticklocation = 'right', ticks=self.qpflevs)
            cb.ax.set_yticklabels(["{0}".format(i) for i in cb.get_ticks()], **self.kw_ticklabels)  # horizontally oriented colorbar
            cb.set_label('QPF (in.)', fontsize=self.fontsize)
            
            ## Add duration plot
            ax = fig.add_subplot(gs[1, 0:2])
            self.plot_duration_longitude(ax)
            ax.text(x, y, '(b)', ha='left', va='top', transform=ax.transAxes, fontsize=12., backgroundcolor='white', zorder=101)
            
            ## add legend for duration plot
            cbax = fig.add_subplot(gs[1, 2])
            # create custom colorbar for duration plot
            self.plot_duration_cbar(cbax)
            
            ## Add probability plot         
            ax = fig.add_subplot(gs[2, 0:2])
            self.plot_vector_landfall_longitude(ax)
            ax.text(x, y, '(c)', ha='left', va='top', transform=ax.transAxes, fontsize=12., backgroundcolor='white', zorder=101)

            ## Add color bar
            cbax = plt.subplot(gs[2,2]) # colorbar axis
            cb = Colorbar(ax = cbax, mappable = self.cf, orientation = 'vertical', ticklocation = 'right', ticks=self.cflevs[::2])
            cb.ax.set_yticklabels(["{:.0%}".format(i) for i in cb.get_ticks()], **self.kw_ticklabels)  # horizontally oriented colorbar
            cb.set_label('Probability of IVT $\geq$ {0} {1}'.format(self.threshold, self.IVT_units), fontsize=self.fontsize)
            
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
            ax = fig.add_subplot(gs[4:, 1:])
            ax = plot_cw3e_logo(ax, orientation='horizontal')
        
        fig.savefig('%s.%s' %(fname1, fmt), bbox_inches='tight', dpi=fig.dpi)
        fig.savefig('%s.%s' %(fname2, fmt), bbox_inches='tight', dpi=fig.dpi)
        # close figure
        plt.close(plt.gcf())