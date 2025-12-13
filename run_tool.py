#!/usr/bin/env python3
"""
Filename:    run_tool.py
Author:      Deanna Nash, dnash@ucsd.edu

Description:
For GEFS, ECMWF, W-WRF, or 'ECMWF-GEFS':
Output .png files of AR landfall tool plots (contour and vector)
for different IVT thresholds and coastal/foothill/inland points.
"""

import os
import sys
import xarray as xr
from datetime import datetime
import traceback

from utils import clear_tmp_dir
from cw3e_tools import LoadDatasets
from ar_landfall_tool_contour import landfall_tool_contour
from ar_landfall_tool_vector import landfall_tool_vector
from ar_landfall_tool_IVT_mag import landfall_tool_IVT_magnitude


# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------

MODEL_CONFIG = {
    "ECMWF": {
        "locs": ['US-west']*4 + ['SAK']*3 + ['AK']*2,
        "oris": ['latitude']*4 + ['longitude']*3 + ['latitude']*2,
        "ptlocs": ['coast', 'foothills', 'inland', 'intwest',
                   'coast', 'foothills', 'inland',
                   'coast', 'inland']
    },
    "GEFS": {
        "locs": ['US-west']*4 + ['SAK']*3 + ['AK']*2,
        "oris": ['latitude']*4 + ['longitude']*3 + ['latitude']*2,
        "ptlocs": ['coast', 'foothills', 'inland', 'intwest',
                   'coast', 'foothills', 'inland',
                   'coast', 'inland']
    },
    "ECMWF-GEFS": {
        "locs": ['US-west']*4 + ['SAK']*3 + ['AK']*2,
        "oris": ['latitude']*4 + ['longitude']*3 + ['latitude']*2,
        "ptlocs": ['coast', 'foothills', 'inland', 'intwest',
                   'coast', 'foothills', 'inland',
                   'coast', 'inland']
    },
    "W-WRF": {
        "locs": ['US-west']*4 + ['SAK']*2,
        "oris": ['latitude']*4 + ['longitude']*2,
        "ptlocs": ['coast', 'foothills', 'inland', 'intwest',
                   'coast', 'foothills']
    }
}


# ---------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------
def plot_magnitudes(ds_pt, loc, ptloc, model, orientation):
    """Plot control and ensemble mean magnitude figures."""
    for mag_type in ["control", "ensemble_mean"]:
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
            orientation=orientation
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
model = sys.argv[1]        # e.g., "GEFS"
init_date = sys.argv[2]    # e.g., "2025013012"

if model not in MODEL_CONFIG:
    raise ValueError(f"Unknown model: {model}")

cfg = MODEL_CONFIG[model]
locs, oris, ptlocs = cfg["locs"], cfg["oris"], cfg["ptlocs"]

prec = None

startTime = datetime.now()

print("\n===============================================")
print(f" Running AR Landfall Tool for {model} {init_date}")
print("===============================================\n")

# ================================================================
# 0. Remove tmp files
# ================================================================
print('Removing tmp intermediate data files...') 
# Specify the directory and the pattern
tmp_directory = "/data/projects/operations/LandfallTools/ar_landfall_tool/data/tmp/"
clear_tmp_dir(tmp_directory)


# ================================================================
# 1. CREATE ONE LOADER PER MODEL RUN (not per-location)
# ================================================================
if model == "ECMWF-GEFS":
    # Load & compute intermediates separately
    loader_ecmwf, interm_ecmwf = load_intermediate_data(
        "ECMWF", locs, ptlocs, init_date
    )
    loader_gefs, interm_gefs = load_intermediate_data(
        "GEFS", locs, ptlocs, init_date
    )

    # Align and subtract
    interm_ecmwf = interm_ecmwf.drop_vars(["duration", "ensemble"])
    interm_gefs = interm_gefs.drop_vars(["duration", "ensemble"])
    interm_gefs = interm_gefs.sel(forecast_hour=interm_ecmwf.forecast_hour.values)
    print(interm_ecmwf)
    print(interm_gefs)
    interm_ecmwf, interm_gefs = xr.align(interm_ecmwf, interm_gefs, join="exact")
    intermediate = interm_ecmwf - interm_gefs

    # Choose one loader to own the differenced data
    loader = loader_ecmwf
    loader.intermediate = intermediate   # <-- THIS is the key line

else: ## all other model choices
    loader, intermediate = load_intermediate_data(model, locs, ptlocs, init_date)


# Only load precipitation dataset once if the model is GEFS or ECMWF
# These are needed for the vector plots, which we do not compute for W-WRF or ECMWF-GEFS
print("Loading QPF once...")
if model in ("ECMWF", "GEFS"):
    ds_qpf = loader.load_prec_QPF_dataset()  # optional depending on workflow
    print("Elapsed:", datetime.now() - startTime)
    print("Computing IVT ensemble mean for vector plots once...")
    ds_ivt_mean = loader.calc_ivt_mean_for_vector_plots()
    print("Elapsed:", datetime.now() - startTime)

# then for each ptloc just extract and save a small netcdf
print("Extracting ptlocs to save as netcdf..")
for loc, ptloc in zip(locs, ptlocs):
    loader.extract_points_from_intermediate(
        loc=loc,
        ptloc=ptloc,
        out_nc_path=f"data/tmp/intermediate_{model}_{init_date}_{loc}_{ptloc}.nc",
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

        ds_pt = xr.open_dataset(f"data/tmp/intermediate_{model}_{init_date}_{loc}_{ptloc}.nc")

        # Save or plot results
        # -----------------------------------------
        # Magnitude Plots
        # -----------------------------------------
        plot_magnitudes(ds_pt, loc, ptloc, model, ori)

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
                forecast=model, threshold=thres,
                orientation=ori
            )
            contour.create_figure()
            
            
            print("\n--------------------------------------------")
            print(f" Vector | {thres}")
            print("--------------------------------------------")
            print("Elapsed:", datetime.now() - startTime)
            # Vector plot (only for ECMWF/GEFS)
            # we do not compute for W-WRF or ECMWF-GEFS
            if model in ("ECMWF", "GEFS"):
                vector = landfall_tool_vector(
                    ds_pt=ds_pt, ds=ds_ivt_mean, prec=ds_qpf,
                    loc=loc, ptloc=ptloc,
                    forecast=model, threshold=thres,
                    orientation=ori
                )
                vector.create_figure()

        # Clean up before next iteration
        del ds_pt  

    except Exception as e:
        print(f"\nERROR processing {loc}, {ptloc}: {e}")
        traceback.print_exc()
        continue


# ================================================================
# 3. Final Cleanup After Workflow Completes
# ================================================================
if model in ("ECMWF", "GEFS"):
    del ds_ivt_mean
    del ds_qpf

print('Removing tmp intermediate data files...') 
# Specify the directory and the pattern
tmp_directory = "/data/projects/operations/LandfallTools/ar_landfall_tool/data/tmp/"
clear_tmp_dir(tmp_directory)

print("\n===============================================")
print(" Workflow Complete")
print(" Total Time:", datetime.now() - startTime)
print("===============================================\n")
