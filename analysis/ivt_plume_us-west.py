"""
Filename:    ivt_plume_us-west-coast.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: For GFS, output .png files of ar landfall tool plot (contour and vector) for different IVT thresholds and coastal, foothill and inland points
"""

import os, sys
sys.path.append('../modules')
from ar_tools import waterfall_ivt_probability

threshold_lst = [150, 250, 500, 750]
textpts_loc_lst = ['coast', 'foothills', 'inland']

for i, threshold in enumerate(threshold_lst):
    for j, textpts_loc in enumerate(textpts_loc_lst):
        s = landfall_tool_contour(loc='US-west', ptloc=textpts_loc, forecast='GEFS', threshold=threshold, orientation='latitude')
        s.create_figure()
        
        s = landfall_tool_vector(loc='US-west', ptloc=textpts_loc, forecast='GEFS', threshold=threshold, orientation='latitude')
        s.create_figure()