"""
Filename:    ar_landfall_tool_AK-south.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: For GFS or ECMWF, output .png files of ar landfall tool plot (contour and vector) for different IVT thresholds and coastal, foothill and inland points
"""

import os, sys
from cw3e_tools import load_datasets
from ar_landfall_tool_contour import landfall_tool_contour
from ar_landfall_tool_vector import landfall_tool_vector
from ar_landfall_tool_IVT_mag import landfall_tool_IVT_magnitude

loc = 'SAK'
ori = 'longitude'
model_lst = ['GEFS', 'ECMWF', 'W-WRF', 'ECMWF-GEFS']
ptloc_lst = ['coast', 'foothills', 'inland']
threshold_lst = [150, 250, 500, 750]

# for each model and point location, load the data, then calculate each metric
for i, (model, ptloc) in enumerate(zip(model_lst, ptloc_lst)):
    
    ################################
    ### Create Intermediate Data ###
    ################################
    ## load the data - this loads and calculates for all metrics and IVT thresholds for the given model and pt location
    s = load_datasets(model, loc, ptloc)
    ds_pt, ds = s.calc_ivt_vars()
    
    if model == 'ECMWF' or model == 'GEFS':
        prec = s.load_prec_QPF_dataset()
    
    ##########################################
    ### Create Plots for Intermediate Data ###
    ##########################################
    ## plot control magnitude plots (this doesn't need to loop through thresholds)
    s = landfall_tool_IVT_magnitude(ds_pt=ds_pt, loc=loc, ptloc=ptloc, forecast=model, mag_type='control', orientation=ori)
    s.create_figure()
    
    ## plot ensemble mean magnitude plots (this doesn't need to loop through thresholds)
    s = landfall_tool_IVT_magnitude(ds_pt=ds_pt, loc=loc, ptloc=ptloc, forecast=model, mag_type='ensemble_mean', orientation=ori)
    s.create_figure()
    
    ## plot vector and contour landfall plots (this will loop through all the thresholds)
    for j, thres in enumerate(threshold_lst):
        s = landfall_tool_contour(ds_pt=ds_pt, loc=loc, ptloc=ptloc, forecast=model, threshold=thres, orientation=ori)
        s.create_figure()
        
        if model == 'ECMWF' or model == 'GEFS':
            s = landfall_tool_vector(ds_pt=ds_pt, ds=ds, prec=prec, loc=loc, ptloc=ptloc, forecast=model, threshold=thres, orientation=ori)
            s.create_figure()