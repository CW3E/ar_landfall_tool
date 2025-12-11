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
from concurrent.futures import ThreadPoolExecutor

from cw3e_tools import load_datasets
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

def load_intermediate_data(model, loc, ptloc, init_date):
    """Loads and computes IVT variables for the chosen model."""
    if model == "ECMWF-GEFS":
        s_ecmwf = load_datasets('ECMWF', loc, ptloc, init_date)
        s_gefs = load_datasets('GEFS', loc, ptloc, init_date)

        # Parallel computation
        with ThreadPoolExecutor() as executor:
            fut_ec = executor.submit(s_ecmwf.calc_ivt_vars)
            fut_ge = executor.submit(s_gefs.calc_ivt_vars)

            ds_pt_ec, ds_ec = fut_ec.result()
            ds_pt_ge, ds_ge = fut_ge.result()

        # ECMWF - GEFS
        ds_pt = ds_pt_ec - ds_pt_ge
        ds = ds_ec - ds_ge

        # Add metadata
        ds_pt = ds_pt.assign_attrs(model_init_date=init_date)

    else:
        s = load_datasets(model, loc, ptloc, init_date)
        ds_pt, ds = s.calc_ivt_vars()

    return ds_pt, ds


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

model = sys.argv[1]
init_date = sys.argv[2]

if model not in MODEL_CONFIG:
    raise ValueError(f"Unknown model: {model}")

cfg = MODEL_CONFIG[model]
locs, oris, ptlocs = cfg["locs"], cfg["oris"], cfg["ptlocs"]

prec = None

for i, (loc, ori, ptloc) in enumerate(zip(locs, oris, ptlocs)):
    print(model, loc, ori, ptloc, ":", datetime.now() - startTime)

    # -----------------------------------------
    # Load intermediate data
    # -----------------------------------------
    ds_pt, ds = load_intermediate_data(model, loc, ptloc, init_date)

    # Only load precipitation dataset once
    if i == 0 and model in ("ECMWF", "GEFS"):
        # s exists only for single-model workflow
        s = load_datasets(model, loc, ptloc, init_date)
        prec = s.load_prec_QPF_dataset()

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
                ds_pt=ds_pt, ds=ds, prec=prec,
                loc=loc, ptloc=ptloc,
                forecast=model, threshold=thres,
                orientation=ori
            )
            vector.create_figure()

print("Time to execute:", datetime.now() - startTime)
