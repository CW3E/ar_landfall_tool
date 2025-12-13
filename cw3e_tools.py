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
# from concurrent.futures import ThreadPoolExecutor
from typing import Tuple, Optional, List, Sequence, Dict


from profiler import StepProfiler

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


class LoadDatasets:
    """
    Refactored loader that caches on a per-(forecast, init_date) basis to avoid
    repeatedly opening and preprocessing large .nc files.

    Usage:
        loader = LoadDatasets('ECMWF', 'US-west', 'coast', '2025121100')
        ds_pt, ds_full = loader.calc_ivt_vars()

    Public methods:
        - calc_ivt_vars() -> (ds_pt, ds_full)
        - load_prec_QPF_dataset() -> precipitation dataset (cached)
    """

    # Class-level caches (shared across instances)
    _cached_ds_model = {}         # key: (forecast, init_date) -> preprocessed full ds
    _cached_vector_mean = {}      # key: (forecast, init_date) -> precomputed ensemble mean vector ds (optional)
    _cached_prec = {}             # key: (forecast, init_date, forecast) -> prec dataset

    def __init__(self, model: str, loc: str, ptloc: str, init_date: str):

        self.forecast = model
        self.loc = loc
        self.ptloc = ptloc
        self.model_init_date = init_date  # format "YYYYmmddHH"

        # paths and model-specific settings
        self._setup_paths_and_settings()

        # ensure the (big) dataset is loaded and preprocessed once
        key = (self.forecast, self.model_init_date)
        if key not in LoadDatasets._cached_ds_model:
            LoadDatasets._cached_ds_model[key] = self.read_ivt_data()
        self.ds_full: xr.Dataset = LoadDatasets._cached_ds_model[key]

    # --------------------------
    # Internal helpers
    # --------------------------
    def _setup_paths_and_settings(self):
        """Set file paths and model-specific attributes."""
        self.path_to_out = '/data/projects/operations/LandfallTools/ar_landfall_tool/data/'

        if self.forecast == 'GEFS':
            self.fpath = '/data/projects/derived_products/GEFS_IVT/data/'
            self.fname = os.path.join(self.fpath, f'GEFS_IVT_{self.model_init_date}.nc')
            self.ensemble_name = 'GEFS'
            self.datasize_min = 15.0

        elif self.forecast in ('ECMWF', 'ECMWF-GEFS'):
            self.fpath = '/data/projects/derived_products/ECMWF_IVT/Ensemble/'
            self.fname = os.path.join(self.fpath, f'IVT_EC_{self.model_init_date}.nc')
            self.ensemble_name = 'ECMWF'
            self.datasize_min = 25.0

        elif self.forecast == 'W-WRF':
            self.fpath = '/data/downloaded/WWRF-NRT/2025-2026/Ensemble_IVT/'
            self.fname = os.path.join(self.fpath, f'IVT_WWRF_{self.model_init_date}.nc')
            self.ensemble_name = 'West-WRF'
            self.datasize_min = 50.0

        else:
            raise ValueError("Forecast product not available! Choose GEFS, ECMWF, ECMWF-GEFS, or W-WRF.")

    # --------------------------
    # Dataset IO & preprocessing
    # --------------------------
    def read_ivt_data(self) -> xr.Dataset:
        """
        Open and preprocess the IVT dataset once per (model, init_date).
        Operations include: open, cast coords, normalize lon to -180..180,
        rename dims/coords to a consistent set, and any model-specific fixes.
        """
        print(f"Opening and preprocessing IVT dataset for {self.forecast} {self.model_init_date}: {self.fname}")

        # try opening the dataset
        try:
            ds = xr.open_dataset(self.fname)
        except Exception as e:
            # raise a helpful error so upstream code can handle missing files
            raise OSError(f"Unable to open IVT file {self.fname}: {e}")

        # normalize longitudes to [-180, 180)
        if 'lon' in ds.coords or 'longitude' in ds.coords:
            # prefer 'lon' coordinate name but support both
            lon_name = 'lon' if 'lon' in ds.coords else 'longitude'
            ds = ds.assign_coords({lon_name: ((ds[lon_name] + 180) % 360 - 180)}).sortby(lon_name)
            # if needed, rename to 'lon'/'lat' canonical names
            if 'longitude' in ds.coords and 'lon' not in ds.coords:
                ds = ds.rename({'longitude': 'lon'})
            if 'latitude' in ds.coords and 'lat' not in ds.coords:
                ds = ds.rename({'latitude': 'lat'})

        # Model-specific adjustments
        if self.forecast == 'ECMWF':
            if 'forecast_time' in ds:
                ds = ds.rename({'forecast_time': 'forecast_hour'})
            # historic: for years > 2020 forecast_hour was stored differently
            try:
                dt_init = datetime.datetime.strptime(self.model_init_date, "%Y%m%d%H")
                if int(dt_init.strftime('%Y')) > 2020 and 'forecast_hour' in ds:
                    ds['forecast_hour'] = ds['forecast_hour'] * 3
            except Exception:
                pass

        if self.forecast == 'W-WRF':
            if 'ensembles' in ds:
                ds = ds.rename({'ensembles': 'ensemble'})

        # Guarantee standard coordinate names: forecast_hour, ensemble, location, lat, lon
        # (if certain names don't exist downstream code should handle gracefully)
        return ds

    def download_QPF_dataset(self):
        """
        Keep your original QPF download/copy logic here.
        This function remains available for load_prec_QPF_dataset fallback usage.
        """
        dt_init = datetime.datetime.strptime(self.model_init_date, "%Y%m%d%H")
        date = dt_init.strftime('%Y%m%d')  # model init date
        hr = dt_init.strftime('%H')  # model init hour

        if self.forecast == 'GEFS':
            print(date, hr)
            # shell=True kept for legacy usage in original code; consider removing for safety later
            subprocess.check_call(["download_QPF.sh", date, hr], shell=True)
        else:
            mmdyhr_init = dt_init.strftime('%m%d%H')
            date2 = dt_init + datetime.timedelta(days=7)
            date2 = date2.strftime('%m%d%H')
            fpath = f'/data/downloaded/Forecasts/ECMWF/NRT_data/{date}{hr}/'
            fname = f'S1D{mmdyhr_init}00{date2}001'
            shutil.copy(fpath + fname, os.path.join(self.path_to_out, 'precip_ECMWF'))

    def load_prec_QPF_dataset(self) -> xr.Dataset:
        """
        Load precipitation (QPF) dataset. Cached per (forecast, init_date).
        """
        key = (self.forecast, self.model_init_date)
        if key in LoadDatasets._cached_prec:
            return LoadDatasets._cached_prec[key]

        # GEFS path: try remote open first, then local fallback
        if self.forecast == 'GEFS':
            try:
                dt_init = datetime.datetime.strptime(self.model_init_date, "%Y%m%d%H")
                date = dt_init.strftime('%Y%m%d')
                hr = dt_init.strftime('%H')
                url = f'https://nomads.ncep.noaa.gov/dods/gfs_0p25/gfs{date}/gfs_0p25_{hr}z'
                ds = xr.open_dataset(url, decode_times=False)
                ds = ds.isel(time=7 * 8)  # 7-day QPF
                prec = ds['apcpsfc'] / 25.4  # mm -> inches
            except OSError:
                try:
                    self.download_QPF_dataset()
                    ds = xr.open_dataset('precip_GFS.grb', engine='cfgrib', backend_kwargs={"indexpath": ""})
                    ds = ds.rename({'longitude': 'lon', 'latitude': 'lat'})
                    prec = ds['tp'] / 25.4
                except Exception:
                    # fallback: zeros grid
                    lats = np.arange(-90., 90.25, 0.25)
                    lons = np.arange(-180., 180.25, 0.25)
                    var_dict = {'tp': (['lat', 'lon'], np.zeros((len(lats), len(lons))))}
                    prec = xr.Dataset(var_dict, coords={'lat': (['lat'], lats), 'lon': (['lon'], lons)})
        else:
            # ECMWF/WWRF path (unchanged logic, but cached)
            self.download_QPF_dataset()
            var_lst = ['u10', 'lsm', 'msl', 'd2m', 'z', 't2m', 'stl1', 'stl2', 'stl3', 'stl4',
                       'swvl4', 'swvl2', 'swvl3', 'sst', 'sp', 'v10', 'sd', 'skt', 'swvl1',
                       'siconc', 'tcwv', 'tcw']
            ds = xr.open_dataset(
                os.path.join(self.path_to_out, 'precip_ECMWF'),
                drop_variables=var_lst,
                engine='cfgrib',
                backend_kwargs={'filter_by_keys': {'typeOfLevel': 'surface'}}
            )
            prec = ds['tp'] * 39.3701  # meters -> inches
            prec = prec.rename({'longitude': 'lon', 'latitude': 'lat'})

        LoadDatasets._cached_prec[key] = prec
        return prec

    # --------------------------
    # IVT computations
    # --------------------------
    def calc_ivt_mean_for_vector_plots(self, ds: Optional[xr.Dataset] = None) -> xr.Dataset:
        """
        Compute ensemble-mean uIVT/vIVT/IVT mean needed for vector plots for the first 7 days.
        This can be cached per (model, init_date) because it does not depend on ptloc.
        """
        key = (self.forecast, self.model_init_date)
        if key in LoadDatasets._cached_vector_mean:
            return LoadDatasets._cached_vector_mean[key]

        if ds is None:
            ds = self.ds_full

        # slice first 7 days (0..24*7), compute mean across forecast_hour & ensemble
        ds_small = ds.sel(forecast_hour=slice(0, 24 * 7)).astype('float64')
        ensemble_mean = ds_small.mean(dim=['forecast_hour', 'ensemble']).astype('float32')
        LoadDatasets._cached_vector_mean[key] = ensemble_mean
        return ensemble_mean

    def compute_intermediate_products(
        self,
        ds: Optional[xr.Dataset] = None,
        thresholds: Optional[Sequence[int]] = None,
        duration_multiplier: int = 3,
        chunking: Optional[Dict[str,int]] = None,
        compute: bool = True,
        compressor=None
    ) -> str:
        """
        Compute intermediate IVT products on the full grid (lazy/dask).

        Args
        ----
        ds : xarray.Dataset, optional
            Full IVT dataset. If None, uses self.ds_full (must be set).
        thresholds : sequence of int
            IVT thresholds to compute (default [100,150,250,500,750,1000]).
        duration_multiplier : int
            Multiplier to convert counts to hours (default 3 for 3-hour timesteps).
        chunking : dict
            Chunk sizes for dask (e.g. {'ensemble': -1, 'lat': 200, 'lon':200}).
            If None, a sensible default is chosen.
        compute : bool
            If True, triggers `.to_zarr(..., compute=True)` which executes the dask graph.
            If False, the zarr metadata is written lazily (requires client/workers to compute later).
        compressor : numcodecs compressor instance or None
            Optional compressor for Zarr store (e.g., zarr.Blosc(cname='zstd', clevel=3)).
        Returns
        -------
        """

        # Defaults
        if ds is None:
            if not hasattr(self, 'ds_full') or self.ds_full is None:
                raise ValueError("Full dataset not provided and self.ds_full is not set.")
            ds = self.ds_full

        if thresholds is None:
            thresholds = [100, 150, 250, 500, 750, 1000]

#         if out_zarr_path is None:
#             out_zarr_path = f"data/tmp/ivt_intermediate_{self.forecast}_{self.model_init_date}.zarr"
#         print('Chunking data...')
#         # sensible chunking defaults (tweak for your machine)
#         if chunking is None:
#             chunking = {}
#             # try to preserve existing dims where available
#             if 'ensemble' in ds.dims:
#                 chunking['ensemble'] = -1   # keep ensemble as one chunk (good for mean reductions)
#             if 'forecast_hour' in ds.dims:
#                 chunking['forecast_hour'] = min(24*7, ds.sizes.get('forecast_hour', 24))
#             if 'lat' in ds.dims:
#                 chunking['lat'] = 200
#             if 'lon' in ds.dims:
#                 chunking['lon'] = 200

#         # Re-chunk dataset for efficient reductions
#         ds = ds.chunk(chunking)
        
        print('Computing thresholds and duration...')
        # Compute data_size (number of non-missing ensembles) per (forecast_hour, lat, lon)
        ivt = ds['IVT']
        if 'ensemble' not in ivt.dims:
            raise KeyError("IVT must have 'ensemble' dimension for these computations.")
        data_size = ivt.count(dim='ensemble')  # (forecast_hour, lat, lon)
        valid_mask = data_size >= self.datasize_min  # boolean (forecast_hour, lat, lon)

        # Build threshold DataArray that will broadcast to (threshold, forecast_hour, ensemble, lat, lon)
        thr = np.asarray(thresholds)
        thr_da = xr.DataArray(thr, dims=['threshold'])

        # Broadcast thresholds across ivt by adding a threshold dim (xarray will broadcast)
        # mask shape: (threshold, forecast_hour, ensemble, lat, lon)
        mask = ivt >= thr_da

        # Probability: fraction of ensembles >= threshold -> mean over ensemble axis
        # result dims: (threshold, forecast_hour, lat, lon)
        probability = mask.mean(dim='ensemble')  # lazy

        # Duration: count of forecast_hour where condition true -> sum over forecast_hour then * multiplier
        # First sum over forecast_hour: dims (threshold, ensemble, lat, lon) -> then mean over ensemble
        # We want duration per (threshold, lat, lon) aggregated across forecast_hour.
        duration = mask.sum(dim='forecast_hour') * duration_multiplier  # dims (threshold, ensemble, lat, lon)
        # compute duration as mean across ensemble after valid_mask is applied:
        duration = duration.mean(dim='ensemble')  # dims: (threshold, lat, lon)
        
        valid_loc = data_size.max(dim='forecast_hour') >= self.datasize_min  # dims (lat, lon)
        # Broadcast valid_loc to probability and duration (threshold, lat, lon)
        probability = probability.where(valid_loc)
        duration = duration.where(valid_loc)
        
        # Explicit threshold coordinate (CRITICAL)
        threshold_coord = xr.DataArray(
            thresholds,
            dims="threshold",
            name="threshold"
        )

        probability = probability.assign_coords(threshold=threshold_coord)
        duration = duration.assign_coords(threshold=threshold_coord)

        print('Calculating ensemble means...')
        # Ensemble means (reduce FIRST)
        ensemble_mean_ivt = ds['IVT'].mean(dim='ensemble')
        ensemble_mean_u   = ds['uIVT'].mean(dim='ensemble')
        ensemble_mean_v   = ds['vIVT'].mean(dim='ensemble')

        # Spatial validity mask
        valid_loc = data_size.max(dim='forecast_hour') >= self.datasize_min

        ensemble_mean_ivt = ensemble_mean_ivt.where(valid_loc)
        ensemble_mean_u   = ensemble_mean_u.where(valid_loc)
        ensemble_mean_v   = ensemble_mean_v.where(valid_loc)

        # Normalize vectors safely
        eps = 1e-6
        u_norm = xr.where(ensemble_mean_ivt > eps,
                          ensemble_mean_u / ensemble_mean_ivt,
                          np.nan)
        v_norm = xr.where(ensemble_mean_ivt > eps,
                          ensemble_mean_v / ensemble_mean_ivt,
                          np.nan)

        # Control member
        control = ds['IVT'].isel(ensemble=0)


        # Assemble intermediate dataset
        print('Assemble intermediate dataset...')
        # Align dims: probability dims (threshold, forecast_hour, lat, lon) but we applied mean(dim='ensemble') so it is (threshold, forecast_hour, lat, lon)
        # To pack into intermediate dataset consistent with your previous structure:
        # We'll aggregate probability across forecast_hour by taking max or keep forecast_hour dimension? Your prior code produced prob shape (forecast_hour, location)
        # To keep generality, we'll keep forecast_hour dimension in probability (i.e., probability per forecast_hour)
        # But earlier under mask.mean(dim='ensemble') probability preserves forecast_hour. We applied valid_loc (no forecast_hour dim); that's ok.
        # For duration we have (threshold, lat, lon). We'll add 'location' later when subsetting.

        intermediate = xr.Dataset({
            "probability": probability,               # dims: (threshold, forecast_hour, lat, lon)
            "duration": duration,                     # dims: (threshold, lat, lon)
            "ensemble_mean": ensemble_mean_ivt,       # dims: (forecast_hour?, lat, lon) -> ensemble reduced; here we averaged across ensemble but kept forecast_hour
            "u": u_norm,                              # dims: same as ensemble_mean
            "v": v_norm,                              # dims: same as ensemble_mean
            "control": control                        # dims: (forecast_hour, lat, lon, ensemble removed)
        })

        # Add attributes
        intermediate = intermediate.assign_attrs({
            "source_forecast": self.forecast,
            "model_init_date": self.model_init_date,
            "thresholds": thresholds
        })
        
        self.intermediate = intermediate
        
#         print('Writing to ZARR...')
#         # Optionally set compressor on each variable via encoding when writing, or leave to xarray default
#         # Ensure output dir exists
#         out_dir = os.path.dirname(out_zarr_path)
#         if out_dir and not os.path.exists(out_dir):
#             os.makedirs(out_dir, exist_ok=True)

#         # Write to zarr (lazy) or compute & persist
#         # If compute==False, we still write metadata references but not compute values
#         write_kwargs = {"mode": "w"}

#         # Provide encoding/chunking hints
#         encoding = {}
#         for vname in intermediate.data_vars:
#             encoding[vname] = {"chunks": None}
#             if compressor is not None:
#                 encoding[vname]["compressor"] = compressor

#         # Save to zarr (this triggers computation if compute==True)
#         # Use safe remove if existing
#         if os.path.exists(out_zarr_path):
#             import shutil
#             shutil.rmtree(out_zarr_path)

#         # Use to_zarr - this will execute lazily or compute depending on compute flag and scheduler config
#         intermediate.to_zarr(out_zarr_path, mode="w", compute=compute, consolidated=True)

        return intermediate
    
    def extract_points_from_intermediate(
        self,
        loc: str,
        ptloc: str,
        out_nc_path: Optional[str] = None,
        method: str = "nearest",
        save_nc: bool = True
    ) -> xr.Dataset:
        """
        Open the intermediate data, subset to points defined by (loc, ptloc),
        and optionally save the result as a small NetCDF for plotting.

        Args
        ----
        loc : str
          Location group (used to find lat/lons file path via self.path_to_out)
        ptloc : str
          Transect name (coast, foothills, inland).
        out_nc_path : str, optional
          Path to save extracted netCDF. If None, will be:
          f"intermediate_{self.forecast}_{self.model_init_date}_{loc}_{ptloc}.nc"
        method : str
          Selection method for xr.sel (default 'nearest')
        save_nc : bool
          If True, save extracted dataset to netCDF (small file).

        Returns
        -------
        ds_pt : xarray.Dataset
          Point-subset dataset (contains probability, duration, ensemble_mean, u, v, control)
        """

        # pull intermediate data
        ds = self.intermediate

        # read lat/lon file for the ptloc
        textpts_fname = os.path.join(self.path_to_out, f'{loc}/latlon_{ptloc}.txt')
        df = pd.read_csv(textpts_fname, header=None, sep=' ', names=['latitude', 'longitude'], engine='python')
        df['longitude'] = df['longitude'] * -1
        lons = df['longitude'].values
        lats = df['latitude'].values

        # Build DataArray for selection
        x = xr.DataArray(lons, dims=['location'])
        y = xr.DataArray(lats, dims=['location'])

        # For probability which may have threshold & forecast_hour, use multi-dim selection
        # We'll select lat/lon by nearest
        ds_pt = ds.sel(lon=x, lat=y, method=method)

        # Optional: add coords for the point names or indices
        ds_pt = ds_pt.assign_coords(location=np.arange(len(lons)))
        ds_pt = ds_pt.assign_coords(ptloc=np.array([ptloc]))

        if save_nc:
            # Save small netCDF (no compression by default)
            out_dir = os.path.dirname(out_nc_path)
            if out_dir and not os.path.exists(out_dir):
                os.makedirs(out_dir, exist_ok=True)
            ds_pt.to_netcdf(out_nc_path)

        return ds_pt
