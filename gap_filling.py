# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 15:36:03 2015

Generic imputation routines

@author: ian_mchugh
"""
import numpy as np
import pdb
from scipy.interpolate import griddata

def generic_2d_linear(data_2d):

    """
    Takes a 2d array as input and;
     1) tiles this into a 3 x 3 space (9 repeats of the original 2d array in 3 columns and 3 rows)
     2) removes the missing data (c.missing_value) from the tiled array
     3) does a bi-linear interpolation to replace the the missing data
     4) returns the central tile
     The effect is to replace missing data in the original 2d array with data from a bi-linear
     interpolation, the tiling repeats the original array along its boundaries to avoid problems
     at the array edges.
    """
    
    data_2d_tiled = np.tile(data_2d, (3,3))   
    
    num_x = np.shape(data_2d_tiled)[1]
    flat_x = np.arange(0, num_x)
    num_y = np.shape(data_2d_tiled)[0]
    flat_y = np.arange(0, num_y)

    # Make the regular grid to project the data onto
    coords_x, coords_y = np.meshgrid(flat_x, flat_y)
    
    # Make a flat array of the tiled data
    data_flat = data_2d_tiled.flatten()

    # Define an index that will return all valid data for the array
    index = np.where(~np.isnan(data_flat))
    
    # Generate a 2d array with existing 2d coordinates of the tiled data
    data_coords = np.column_stack([coords_x.flatten(), 
                                   coords_y.flatten()])

    # Do the interpolation
    grid_z = griddata(data_coords[index], data_flat[index], 
                      (coords_x, coords_y), method = 'linear')
    
    # Return the central tile
    return grid_z[num_y / 3: num_y / 3 * 2, num_x / 3: num_x / 3 * 2]
   
