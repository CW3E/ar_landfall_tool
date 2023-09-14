"""
Filename:    ar_landfall_tool.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: For GFS or ECMWF, output .png files of ar landfall tool plot (contour and vector) for different IVT thresholds and coastal, foothill and inland points
"""

import os, sys
sys.path.append('../modules')
from cw3e_tools import landfall_tool_contour, landfall_tool_vector

ptloc_lst = ['coast', 'inland'] + ['coast', 'foothills', 'inland']*2
loc_lst = ['AK']*2 + ['SAK']*3 + ['US-west_new']*3
ori_lst = ['latitude']*2 + ['longitude']*3 + ['latitude']*3
ds = []
for i, (ptloc, loc, ori) in enumerate(zip(ptloc_lst, loc_lst, ori_lst)):
    s = landfall_tool_contour(loc=loc, ptloc=ptloc, forecast='ECMWF', threshold=150, orientation=ori)
    s.create_figure()
    
    s = landfall_tool_vector(loc=loc, ptloc=ptloc, forecast='GEFS', threshold=150, orientation=ori)
    s.create_figure()