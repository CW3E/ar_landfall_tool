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
from datetime import datetime
import traceback

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
# 1. CREATE ONE LOADER PER MODEL RUN (not per-location)
# ================================================================
# We temporarily initialize with dummy loc/ptloc; these get updated later
loader = LoadDatasets(model, locs[0], ptlocs[0], init_date)

print("Reading IVT dataset once...")
ds_ivt = loader.read_ivt_data()         # <-- cached internally & reused everywhere

# Only load precipitation dataset once if the model is GEFS or ECMWF
print("Loading QPF once...")
if model in ("ECMWF", "GEFS"):
    ds_qpf = loader.load_prec_QPF_dataset()  # optional depending on workflow

# ================================================================
# 2. MAIN LOOP OVER LOCATIONS USING THE SAME ds_ivt
# ================================================================
for i, (loc, ori, ptloc) in enumerate(zip(locs, oris, ptlocs)):
    print("\n--------------------------------------------")
    print(f" {i+1}/{len(locs)} :: {model} | {loc} | {ptloc}")
    print("--------------------------------------------")
    print("Elapsed:", datetime.now() - startTime)

    try:
        # -----------------------------------------
        # Update location info inside loader
        # -----------------------------------------
        loader.loc = loc
        loader.ptloc = ptloc
        loader.get_lat_lons_from_txt_file()

        # -----------------------------------------
        # Compute point-based probabilities
        # Using the SAME ds_ivt loaded once above
        # -----------------------------------------
        ds_pt = loader.calc_ivt_probability_and_duration_for_points(ds_ivt)

        # Save or plot results
        # -----------------------------------------
        # Magnitude Plots
        # -----------------------------------------
        plot_magnitudes(ds_pt, loc, ptloc, model, ori)

        # -----------------------------------------
        # Contour + Vector Plots for thresholds
        # -----------------------------------------
        for thres in threshold_list(ptloc):

            # Contour plot
            contour = landfall_tool_contour(
                ds_pt=ds_pt, loc=loc, ptloc=ptloc,
                forecast=model, threshold=thres,
                orientation=ori
            )
            contour.create_figure()

            # Vector plot (only for ECMWF/GEFS)
            if model in ("ECMWF", "GEFS"):
                vector = landfall_tool_vector(
                    ds_pt=ds_pt, ds=ds, prec=ds_qpf,
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
print("\nReleasing internal dataset cache...")
loader.release_ivt_dataset()

del ds_ivt
del ds_qpf

print("\n===============================================")
print(" Workflow Complete")
print(" Total Time:", datetime.now() - startTime)
print("===============================================\n")
