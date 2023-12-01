"""
Filename:    run_tool.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: For GEFS, ECMWF, W-WRF, or 'ECMWF-GEFS': output .png files of ar landfall tool plot (contour and vector) for different IVT thresholds and coastal, foothill and inland points
"""

import os, sys
from cw3e_tools import load_datasets
from ar_landfall_tool_contour import landfall_tool_contour
from ar_landfall_tool_vector import landfall_tool_vector
from ar_landfall_tool_IVT_mag import landfall_tool_IVT_magnitude
from datetime import datetime
# from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['PT Sans']})
# rc('text', usetex=True)


startTime = datetime.now() # get start time of script
model = sys.argv[1]

if (model == 'ECMWF') | (model == 'GEFS') | (model == 'ECMWF-GEFS'):
    loc_lst = ['US-west']*4 + ['SAK']*3 + ['AK']*2
    ori_lst = ['latitude']*4 + ['longitude']*3 + ['latitude']*2
    ptloc_lst = ['coast', 'foothills', 'inland', 'intwest']*2 + ['coast', 'inland']

elif (model == 'W-WRF'):
    loc_lst = ['US-west']*4 + ['SAK']*2
    ori_lst = ['latitude']*4 + ['longitude']*2 
    ptloc_lst = ['coast', 'foothills', 'inland', 'intwest'] + ['coast', 'foothills']

threshold_lst = [150, 250, 500, 750]

# for each model and point location, load the data, then calculate each metric
for i, (loc, ori, ptloc) in enumerate(zip(loc_lst, ori_lst, ptloc_lst)):
    print(model, loc, ori, ptloc)
    
    ################################
    ### Create Intermediate Data ###
    ################################
    ## load the data - this loads and calculates for all metrics and IVT thresholds for the given model and pt location
    if model == 'ECMWF-GEFS':   
        s = load_datasets('ECMWF', loc, ptloc)
        ds_pt_ECMWF, ds_ECMWF = s.calc_ivt_vars()
        model_init_date = ds_pt_ECMWF.model_init_date
        date_string = model_init_date.strftime('%Y%m%d%H')
        path_to_data = '/data/downloaded/SCRATCH/cw3eit_scratch/'
        fname = path_to_data + 'GEFS/FullFiles/IVT_Full_{0}.nc'.format(date_string)

        s = load_datasets('GEFS', loc, ptloc, fname)
        ds_pt_GEFS, ds_GEFS = s.calc_ivt_vars()

        ## subtract ECMWF - GEFS
        ds_pt = ds_pt_ECMWF-ds_pt_GEFS
        ds = ds_ECMWF-ds_GEFS
        
        ##Add attribute information
        ds_pt = ds_pt.assign_attrs(model_init_date=model_init_date)
        
    else:
        s = load_datasets(model, loc, ptloc)
        ds_pt, ds = s.calc_ivt_vars()
    
    if (i == 0 and model == 'ECMWF') or (i == 0 and model == 'GEFS'):
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

print('Time to execute:', datetime.now() - startTime)