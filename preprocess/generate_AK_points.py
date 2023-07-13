"""
Filename:    generate_AK_points.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Functions for generating text file of points for AK AR Landfall Tools
"""
# Standard Python modules
import os, sys
import glob
import numpy as np
import pandas as pd


##################
### West Coast ###
##################

lats_coast = np.arange(54.0, 71.5, 0.5)
lons_coast = [165.5, 165.0, 163.5, 162.5, 160.5, 159.5, 158.5, 158.0, 
              157.5, 157.5, 162.0, 162.0, 162.0, 165.0, 165.5, 166.0, 
              166.0, 165.5, 165.0, 162.5, 161.0, 166.5, 166.5, 168.0, 
              167.0, 165.5, 163.0, 164.0, 165.0, 166.5, 166.0, 163.0, 
              162.5, 160.5, 157.5]

d = {'lat': lats_coast, 'lon': lons_coast}
df = pd.DataFrame(data=d)
np.savetxt(r'../data/AK/latlon_coast.txt', df.values, fmt='%.2f')


lats_inland = np.arange(57.5, 71.5, 0.5)
lons_inland = [157., 156.5, 156.5, 161.5, 161.5, 161.5, 163., 163., 
               164., 164., 164., 162., 160.5, 160.5, 160.5, 166., 
               165.5, 166., 160., 161., 163., 164., 165.5, 163., 
               162., 161.5, 159., 156.5]

d = {'lat': lats_inland, 'lon': lons_inland}
df = pd.DataFrame(data=d)
np.savetxt(r'../data/AK/latlon_inland.txt', df.values, fmt='%.2f')

##############################
### South- Southeast Coast ###
##############################
lons = np.arange(130., 165.5, .5)
lats_coast =[54.0, 54.5, 54.5, 55.0, 55.0, 55.,  55.0, 55.5,  56., 56.,  56.5, 57.,  57.5,
             58.0, 58.5, 58.5, 59.0, 59.0,  59.5, 59.5,  60.0, 60.,  60.0, 60.,  60.,
             60.0, 60.0, 60.0, 60.0, 60.0,  60.5, 60.5,  60.5, 60.5,  61.0, 61.,  60.5,
             60.0, 60.0, 60.0, 59.5, 59.5,  59.5, 59.,  59.5, 60.,  60.0, 59.5,  58.5,
             58.0, 58.0, 57.5, 57.5, 57.0,  57.0, 56.5,  56.5, 56.,  56.0, 55.5,  55.5,
             55.5, 55.5, 55.5, 55.0, 55.0,  55.0, 54.5,  54.5, 54.5, 54.5]

d = {'lat': lats_coast, 'lon': lons}
df = pd.DataFrame(data=d)
np.savetxt(r'../data/SAK/latlon_coast.txt', df.values, fmt='%.2f')

lats_inland = [55.0, 55.5, 55.5, 56.0, 56.0, 56.0, 56.0, 56.5, 57.0, 57.0, 
               57.5, 58.0, 58.5, 59.0, 59.5, 59.5, 60.0, 60.0, 60.5, 60.5, 
               60.5, 61.0, 61.0, 61.0, 61.0, 61.0, 61.0, 61.0, 61.0, 61.0, 
               61.5, 61.5, 61.5, 61.5, 62.0, 62.0, 61.5, 61.5, 61.0, 61.0, 
               60.5, 60.5, 60.5, 61.5, 61.5, 61.0, 61.0, 60.5, 60.0, 59.5, 
               59.0, 58.5, 58.0, 58.0, 57.5, 57.0, 57.0, 56.5, 56.5, 56.0, 
               56.0, 56.0, 56.0, 56.0, 55.5, 55.5, 55.0, 55.0, 55.0, 55.0, 
               54.5]

d = {'lat': lats_inland, 'lon': lons}
df = pd.DataFrame(data=d)
np.savetxt(r'../data/SAK/latlon_inland.txt', df.values, fmt='%.2f')