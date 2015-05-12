# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 15:36:03 2015

Generic imputation routines

@author: ian_mchugh
"""
import numpy as np
import pdb
from scipy.interpolate import griddata as griddata_sc
from matplotlib.mlab import griddata as griddata_np

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
    grid_z = griddata_sc(data_coords[index], data_flat[index], 
                      (coords_x, coords_y), method = 'linear')
    
    # Return the central tile
    return grid_z[num_y / 3: num_y / 3 * 2, num_x / 3: num_x / 3 * 2]

def do_2dinterpolation(array_2d):
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
    
    missing_value = -9999 
    eps = 0.0000001    
    
    WasMA = False
    if np.ma.isMA(array_2d):
        WasMA = True
        array_2d = np.ma.filled(array_2d, float(missing_value))
    
    # Tile the 2d array into a 3 by 3 array
    data_2d_tiled = np.tile(array_2d,(3,3))
    
    # Get the dimensions of the tiled array and create coordinates and grid
    num_x = np.shape(data_2d_tiled)[1]
    array_x = np.arange(0, num_x)
    num_y = np.shape(data_2d_tiled)[0]
    array_y = np.arange(0, num_y)
    coords_x, coords_y = np.meshgrid(array_x, array_y)
    
    # Make a flat array of the tiled data
    data_1d = data_2d_tiled.flatten()
    
    # Make a 2d array of the coordinates
    data_coords = np.column_stack([coords_x.flatten(), 
                                   coords_y.flatten()])
    
    # Define an index that will return all valid data for the array
    index = np.where(data_1d!= missing_value)    
    
    # Do the interpolation
    grid_z = griddata_sc(data_coords[index], data_1d[index], 
                      (coords_x, coords_y), method = 'linear')

    # Retrieve the central tile
    array_2d_filled = grid_z[num_y / 3: num_y / 3 * 2, num_x / 3: num_x / 3 * 2]

    # Check something...
    if WasMA:
        array_2d_filled = np.ma.masked_where(abs(array_2d_filled - np.float64(missing_value)) < eps, array_2d_filled)
        array_2d = np.ma.masked_where(abs(array_2d - np.float64(missing_value)) < eps, array_2d)

    # Return the filled array
    return array_2d_filled   
