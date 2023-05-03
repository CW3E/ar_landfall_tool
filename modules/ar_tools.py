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

# matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colorbar import Colorbar # different way to handle colorbar
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker

# cartopy
import cartopy.crs as ccrs
from cartopy.mpl.geoaxes import GeoAxes
import cartopy.feature as cfeature

# import personal modules
import nclcmaps as nclc

class waterfall_ivt_probability:
    '''
    Returns a .png file with Cordeira "waterfall" figure with input locations, chosen Forecast product, and IVT threshold
    
    Parameters
    ----------
    filename : str
        name of the .txt file with latitude and longitude of locations for analysis
        this file should have no header, with latitude, then longitude, separated by a space
    forecast : str
        name of the forecast product - options include GEFS, ECMWF, TODO: ECMWF - GEFS and West-WRF
    threshold : int
        threshold for IVT probabilty in kg m-1 s-1 - options include 150, 250, 500, 750
  
    Returns
    -------
    fig : figure
        png file of the waterfall figure
    
    '''
    
    def __init__(self, textpts_fname, forecast='GEFS', threshold=250):
        path_to_data = '/data/downloaded/SCRATCH/cw3eit_scratch/'
        if forecast == 'GEFS':
            fpath = path_to_data + 'GEFS/FullFiles/'
            self.ensemble_name = 'GEFSv12'
        elif forecast == 'ECMWF':
            fpath = path_to_data + 'GEFS/FullFiles/'
            self.ensemble_name = 'GEFSv12'
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
        df = pd.read_csv(textpts_fname, header=None, sep=' ', names=['latitude', 'longitude'], engine='python')
        df['longitude'] = df['longitude']*-1
        self.df = df
        self.threshold = threshold
                  
        ## format dicts for plots
        self.kw_ticklabels = {'size': 6, 'color': 'gray', 'fontweight': 'light'}
        self.kw_grid = {'linewidth': .5, 'color': 'k', 'linestyle': '--', 'alpha': 0.1}
        self.kw_ticks = {'length': 4, 'width': 0.5, 'pad': 2, 'color': 'black'}
        self.IVT_units = 'kg m$^{-1}$ s$^{-1}$'

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

        ## calculate probability IVT > threshold
        data_size = ds.IVT.ensemble.shape

        # sum the number of ensembles where IVT exceeds threshold
        self.probability = (ds.IVT.where(ds.IVT > self.threshold).count(dim='ensemble')) / data_size
    
        
    def plot_waterfall(self, ax):
        
        self.calc_ivt_probability() # run the calculation
        
        # Contour Filled
        x = self.probability.forecast_hour / 24 # convert forecast hour to forecast day
        y = self.probability.lat
        data = np.flipud(np.rot90(self.probability.values)) # rotate data 90 degrees and flip up down
        self.cflevs = [0., 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.] # levels for IVT probability
        cmap = nclc.cmap('WhiteBlueGreenYellowRed') # cmap for IVT probability
        self.cf = ax.contourf(x, y, data, levels=self.cflevs, cmap=cmap, extend='neither')
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
    
                  
    def plot_map(self, ax, mapcrs, datacrs):
                  
        ## Set up extent of plot at gridlines
        lonmin = self.lons.min()-5.
        lonmax = self.lons.max()+5.
        latmin = self.lats.min()
        latmax = self.lats.max()

        ext = [lonmin, lonmax, latmin, latmax] # extent of map plot
        du = 5 # how frequent for ticks ##TODO: create a more flexible option
        dx = np.arange(lonmin,lonmax+du,du)
        dy = np.arange(latmin,latmax+du,du)
                  
        ax.set_extent(ext, crs=datacrs)

        # Add map features (continents and country borders)
        ax.add_feature(cfeature.LAND, facecolor='0.9')      
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
                  
        ## TODO: Add elevation contours
        # cf = ax.contourf(lons, lats, elev, transform=datacrs,
        #                  levels=cflevs, cmap=cmap, alpha=0.9, extend='max')
         
        ## plot point locations
        for i, (x, y) in enumerate(zip(self.lons, self.lats)):
            ax.plot(x, y, 'ko', markersize=2, transform=datacrs)

        ## plot labels
        label = 'Forecasts support FIRO/CA-AR Program Intended for research purposes only'
        ax.set_title(textwrap.fill(label, 36), loc='right', fontsize=5)
        
        return ax
    
    def create_figure(self):

        fig = plt.figure(figsize=(10, 3))
        fig.dpi = 300
        fname = '../figs/waterfall_{0}_{1}'.format(self.ensemble_name, self.date_string)
        fmt = 'png'

        nrows = 2
        ncols = 2

        ## Use gridspec to set up a plot with a series of subplots that is
        ## n-rows by n-columns
        gs = GridSpec(nrows, ncols, height_ratios=[1, 0.05], width_ratios = [2, 1], wspace=0.05, hspace=0.4)
        ## use gs[rows index, columns index] to access grids         
        
        ## Add waterfall plot         
        ax = fig.add_subplot(gs[0, 0])
        self.plot_waterfall(ax)
                  
        ## Add color bar
        cbax = plt.subplot(gs[1,0]) # colorbar axis
        cb = Colorbar(ax = cbax, mappable = self.cf, orientation = 'horizontal', ticklocation = 'bottom', ticks=self.cflevs[::2])
                  
        
        # Set up projection information for map
        mapcrs = ccrs.PlateCarree()
        datacrs = ccrs.PlateCarree()
        ## Add map
        ax = fig.add_subplot(gs[0, 1], projection=mapcrs)
        self.plot_map(ax, mapcrs, datacrs)
                  
        ## TODO: Add CW3E logo
        # im = 'path/to/cw3e-logo'
        # ax = fig.add_subplot(gs[1, 1])
        # ax.imshow(im)
        # ax.axis('off')
        
        fig.savefig('%s.%s' %(fname, fmt), bbox_inches='tight', dpi=fig.dpi)
        plt.show()