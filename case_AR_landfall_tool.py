#!/usr/bin/python3
"""
Filename:    case_AR_landfall_tool.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Run the AR landfall tool for a select case. User to update the filepath and figpath.
"""

from cw3e_tools import load_datasets
from ar_landfall_tool_contour import landfall_tool_contour
from ar_landfall_tool_IVT_mag import landfall_tool_IVT_magnitude

# --- USER UPDATE ---
model = 'ECMWF' # model name (WWRF, ECMWF, or GEFS)
datadir = '/data/projects/case_studies/comp_test/case_study_1997/ForecastTools/data/1996122512/' # the dir where your IVT nc file is
fname = 'ECMWF_EPS_IVT_1996122512.nc' # the name of the IVT .nc file
figpath = '/home/dnash/repos/ar_landfall_tool/figs/' # the location where you want the figs to save

loc_lst = ['US-west']*3 # the domain
ori_lst = ['latitude']*3 # latitude or longitude (longitude only for Alaska)
ptloc_lst = ['coast', 'foothills', 'inland'] # point location - coast, foothills, or inland only

# --- load data and create plots --- 
# for each model and point location, load the data, then calculate each metric
for i, (loc, ori, ptloc) in enumerate(zip(loc_lst, ori_lst, ptloc_lst)):
    print(model, loc, ori, ptloc)
    s = load_datasets(model, loc, ptloc, fname=datadir+fname)
    ds_pt, ds = s.calc_ivt_vars()
    
    ## plot ensemble mean magnitude plots (this doesn't need to loop through thresholds)
    s = landfall_tool_IVT_magnitude(ds_pt=ds_pt, loc=loc, ptloc=ptloc, forecast=model, mag_type='ensemble_mean', orientation=ori,
                                   path_to_out=figpath)
    s.create_figure()
    
    threshold_lst = [150, 250, 500, 750]
    ## plot contour landfall plots (this will loop through all the thresholds)
    for j, thres in enumerate(threshold_lst):
        s = landfall_tool_contour(ds_pt=ds_pt, loc=loc, ptloc=ptloc, forecast=model, threshold=thres, orientation=ori,
                                   path_to_out=figpath)
        s.create_figure()