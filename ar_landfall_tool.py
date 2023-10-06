"""
Filename:    ar_landfall_tool.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: For GFS or ECMWF, output .png files of ar landfall tool plot (contour and vector) for different IVT thresholds and coastal, foothill and inland points
"""

import os, sys
from cw3e_tools import download_QPF_prec_dataset
from ar_landfall_tool_contour import landfall_tool_contour
from ar_landfall_tool_vector import landfall_tool_vector
from ar_landfall_tool_IVT_mag import landfall_tool_IVT_magnitude

## start by downloading/copying QPF dataset for both GEFS and ECMWF
model_lst = ['ECMWF']
for i, model in enumerate(model_lst):
    print('Downloading {0} 7-d QPF data...'.format(model))
    d = download_QPF_prec_dataset(model)
    d.download_dataset()

ptloc_lst = ['coast', 'inland'] + ['coast', 'foothills', 'inland']*2
loc_lst = ['AK']*2 + ['SAK']*3 + ['US-west']*3
ori_lst = ['latitude']*2 + ['longitude']*3 + ['latitude']*3
threshold_lst = [150, 250, 500, 750]
model_lst = ['GEFS', 'ECMWF', 'W-WRF', 'ECMWF-GEFS']
mag_type_lst = ['control', 'ensemble']

for i, (ptloc, loc, ori) in enumerate(zip(ptloc_lst, loc_lst, ori_lst)):
    for j, thres in enumerate(threshold_lst):
        for k, model in enumerate(model_lst[1:2]): # only run the first two models for the vector tool
            print('Running AR Landfall Vector Tool for {0}, {1}, {2}, {3}...'.format(model, thres, ptloc, loc))
            s = landfall_tool_vector(loc=loc, ptloc=ptloc, forecast=model, threshold=thres, orientation=ori)
            s.create_figure()
            
        for k, model in enumerate(model_lst): # run all models             
            s = landfall_tool_contour(loc=loc, ptloc=ptloc, forecast=model, threshold=thres, orientation=ori)
            s.create_figure()
            
            
for i, (ptloc, loc, ori) in enumerate(zip(ptloc_lst, loc_lst, ori_lst)):
    for j, model in enumerate(model_lst):
        for k, mag in enumerate(mag_type_lst):
            s = landfall_tool_IVT_magnitude(loc=loc, ptloc=ptloc, forecast=model, mag_type=mag, orientation=ori)
            s.create_figure()
            
            
            
################################
### Create Intermediate Data ###
################################

model_lst = ['ECMWF', 'GEFS', 'WWRF']
ptloc_lst = ['coast', 'foothills', 'inland']

for i, model in enumerate(model_lst):
    for j, ptloc in enumerate(ptloc_lst):
        
        ################################
        ### Create Intermediate Data ###
        ################################

        load = load_data(model, ptloc)
        
        ## if ECMWF or GEFS, save
        

        ##########################################
        ### Create Plots for Intermediate Data ###
        ##########################################
        
        s = calculate metric
        
    