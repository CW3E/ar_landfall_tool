"""
Filename:    cw3ecmaps.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Functions for cmaps from the CW3E website adapted from from https://github.com/samwisehawkins/nclcmaps
"""

import numpy as np
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap

__all__ = ['cw3ecmaps']

colors = {"ivt_probability": [[10, 180, 255],
                              [12, 253, 248],
                              [7, 228, 170],
                              [3, 203, 87],
                              [6, 182, 6],
                              [95, 208, 3],
                              [176, 231, 2],
                              [255, 255, 3],
                              [255, 225, 3],
                              [255, 199, 2],
                              [255, 170, 2],
                              [255, 115, 1],
                              [255, 57, 1],
                              [255, 1, 0],
                              [205, 0, 37],
                              [150, 0, 86],
                              [91, 0, 133]],
          
        "ivt_vector": [[255., 255., 255.], #0.0-0.1
                       [9., 156., 255.], #0.1-0.2
                       [8., 235., 191.], #0.2-0.3
                       [2., 186., 24.], #0.3-0.4
                       [155., 225., 2.], #0.4-0.5
                       [255, 233, 3], #0.5-0.6
                       [255, 175, 2], #0.6-0.7
                       [255, 69, 1], #0.7-0.8
                       [215, 0, 28], #0.8-0.9
                       [102, 0, 125]], #0.9-1.0
          
          "jay_qpf": [[10, 193, 255],
                      [156, 47, 206],
                      [206, 33, 33],
                      [140, 67, 0],
                      [255, 128, 1],
                      [140, 101, 1],
                      [255, 187, 7],
                      [255, 223, 140],
                      [255, 131, 173]],
          
          "brian_qpf": [[154, 219, 165],
                        [106, 198, 165],
                        [83, 174, 164],
                        [72, 148, 160],
                        [63, 123, 155],
                        [60, 97, 150],
                        [63, 68, 133],
                        [56, 46, 91],
                        [40, 25, 51]],
          
          "test_terrain": [[255, 255, 255],
                      [156, 47, 206],
                      [206, 33, 33],
                      [140, 67, 0],
                      [255, 128, 1],
                      [140, 101, 1],
                      [255, 187, 7],
                      [255, 223, 140],
                      [255, 131, 173]],
          "land_sea": [[38, 63, 120],
                       [255, 255, 255]],
          "ivt_duration": [[255, 174, 0],
                           [255, 78, 0],
                           [236, 0, 7],
                           [86, 0, 137]]
         }

bounds = {"ivt_probability": [0., 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.],
          "ivt_vector": [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.],
          "jay_qpf": [1, 2, 3, 4, 5, 7, 10, 15, 20],
          "brian_qpf": [1, 2, 3, 4, 5, 7, 10, 15, 20],
          "test_terrain": [0., 250., 500., 1000., 1500., 2000., 2500., 3000., 3500.],
          "land_sea": [0., 1.],
          "ivt_duration": [0, 250, 500, 750, 1000]}


def cmap(name):
    data = np.array(colors[name])
    data = data / np.max(data)
    cmap = ListedColormap(data, name=name)
    bnds = bounds[name]
    norm = mcolors.BoundaryNorm(bnds, cmap.N)
    return cmap, norm, bnds