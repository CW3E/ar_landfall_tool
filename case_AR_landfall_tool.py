#!/usr/bin/python3
"""
Filename:    case_AR_landfall_tool.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Run the AR landfall tool for a select case. User to update the filepath and figpath.
/data/projects/case_studies/comp_test/case_study_1995/ForecastTools/data/1995030512/ECMWF_EPS_IVT_1995030512.nc (plus Mach 6 and 7).
"""
import pandas as pd
import os
import sys
import xarray as xr
import numpy as np
from datetime import datetime
import traceback

from utils import clear_tmp_dir
from cw3e_tools import LoadDatasets
from ar_landfall_tool_contour import landfall_tool_contour
from ar_landfall_tool_IVT_mag import landfall_tool_IVT_magnitude

# --- USER UPDATE ---
# Create date range every 24 hours (1 day)
dates = pd.date_range(
    start="1995-03-05 12:00",
    end="1995-03-07 12:00",
    freq="24H"
)

# Format as YYYYMMDDHH
date_lst = dates.strftime("%Y%m%d%H").tolist()
model="ECMWF_archive"

# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------

MODEL_CONFIG = {
    "ECMWF_archive": {
        "locs": ['US-west']*3,
        "oris": ['latitude']*3,
        "ptlocs": ['coast', 'foothills', 'inland',]
    },
}

# ---------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------
def plot_magnitudes(ds_pt, loc, ptloc, model, orientation, path_to_out):
    """Plot control and ensemble mean magnitude figures."""
    for mag_type in ["ensemble_mean"]:
        print("\n--------------------------------------------")
        print(f" Magnitude | {mag_type}")
        print("--------------------------------------------")
        print("Elapsed:", datetime.now() - startTime)
        fig = landfall_tool_IVT_magnitude(
            ds_pt=ds_pt,
            loc=loc,
            ptloc=ptloc,
            forecast=model,
            mag_type=mag_type,
            orientation=orientation,
            path_to_out=path_to_out,
        )
        fig.create_figure()


def threshold_list(ptloc):
    """Thresholds differ for 'intwest'."""
    return [100, 150, 250, 500, 750] if ptloc == "intwest" else \
           [150, 250, 500, 750]

def load_intermediate_data(model, locs, ptlocs, init_date):
    # We temporarily initialize with dummy loc/ptloc; these get updated later
    loader = LoadDatasets(model, locs[0], ptlocs[0], init_date)

    print("Reading IVT dataset once...")
    ds_full = loader.read_ivt_data()         # <-- cached internally & reused everywhere
    print("Elapsed:", datetime.now() - startTime)
    
    print("Computing intermediate products once")
    # compute intermediate products once (lazy dask)
    intermediate = loader.compute_intermediate_products(
        ds=ds_full,
        thresholds=[100,150,250,500,750,1000],
        chunking={'ensemble': -1, 'forecast_hour': 168, 'lat': 200, 'lon': 200}
    )
    print("Elapsed:", datetime.now() - startTime)
    
    return loader, intermediate

# ---------------------------------------------------------------------
# Main Script
# ---------------------------------------------------------------------

startTime = datetime.now()

# -------------------------------
# Inputs passed to this script
# -------------------------------
if model not in MODEL_CONFIG:
    raise ValueError(f"Unknown model: {model}")

cfg = MODEL_CONFIG[model]
locs, oris, ptlocs = cfg["locs"], cfg["oris"], cfg["ptlocs"]

prec = None

# ================================================================
# 0. Remove tmp files
# ================================================================
print('Removing tmp intermediate data files...') 
# Specify the directory and the pattern
tmp_directory = f"/home/dnash/repos/ar_landfall_tool/data/tmp/{model}/"
clear_tmp_dir(tmp_directory)


for init_date in date_lst:

    startTime = datetime.now()

    print("\n===============================================")
    print(f" Running AR Landfall Tool for {model} {init_date}")
    print("===============================================\n")

    path_to_out = f'/home/dnash/repos/ar_landfall_tool/figs/{init_date}/' # the location where you want the figs to save

    # --- load data and create plots --- 
    # for each model and point location, load the data, then calculate each metric
    loader, intermediate = load_intermediate_data(model, locs, ptlocs, init_date)

    # then for each ptloc just extract and save a small netcdf
    print("Extracting ptlocs to save as netcdf..")
    for loc, ptloc in zip(locs, ptlocs):
        loader.extract_points_from_intermediate(
            loc=loc,
            ptloc=ptloc,
            out_nc_path=f"/home/dnash/repos/ar_landfall_tool/data/tmp/{model}/intermediate_{model}_{init_date}_{loc}_{ptloc}.nc",
            save_nc=True
            )
    print("Elapsed:", datetime.now() - startTime)
    # you can now free memory and later load the small per-ptloc netCDF for plotting
    del intermediate

    # ================================================================
    # 2. Load and Plot Intermediate Data
    # ================================================================
    for i, (loc, ori, ptloc) in enumerate(zip(locs, oris, ptlocs)):
        print("\n--------------------------------------------")
        print(f" {i+1}/{len(locs)} :: {model} | {loc} | {ptloc}")
        print("--------------------------------------------")
        print("Elapsed:", datetime.now() - startTime)
    
        try:
    
            ds_pt = xr.open_dataset(f"/home/dnash/repos/ar_landfall_tool/data/tmp/{model}/intermediate_{model}_{init_date}_{loc}_{ptloc}.nc")
    
            # Save or plot results
            # -----------------------------------------
            # Magnitude Plots
            # -----------------------------------------
            plot_magnitudes(ds_pt, loc, ptloc, "ECMWF", ori, path_to_out)
    
            # -----------------------------------------
            # Contour + Vector Plots for thresholds
            # -----------------------------------------
            for thres in threshold_list(ptloc):
                print("\n--------------------------------------------")
                print(f" Contour | {thres}")
                print("--------------------------------------------")
                print("Elapsed:", datetime.now() - startTime)
    
                # Contour plot
                contour = landfall_tool_contour(
                    ds_pt=ds_pt, loc=loc, ptloc=ptloc,
                    forecast="ECMWF", threshold=thres,
                    orientation=ori,
                    path_to_out=path_to_out,
                )
                contour.create_figure()

                    # Clean up before next iteration
            del ds_pt 

        except Exception as e:
            print(f"\nERROR processing {loc}, {ptloc}: {e}")
            traceback.print_exc()
            continue

    # ================================================================
    # 3. Final Cleanup After Workflow Completes
    # ================================================================
    
    print('Removing tmp intermediate data files...') 
    # Specify the directory and the pattern
    tmp_directory = f"/home/dnash/repos/ar_landfall_tool/data/tmp/{model}/"
    clear_tmp_dir(tmp_directory)
    
    print("\n===============================================")
    print(" Workflow Complete")
    print(" Total Time:", datetime.now() - startTime)
    print("===============================================\n")