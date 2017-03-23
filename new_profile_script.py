#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 09:52:14 2017

@author: ian
"""

import numpy as np
import pandas as pd
import copy as cp
import pdb
import matplotlib.pyplot as plt

default_press = 101.325
site_alt = None
output_int = 30
smooth_window = 10

def downsample_data(df, output_freq = 30, smooth_window = 0):

    """
    This function downsamples profile data to the requested output frequency
    (generally 30 minutes);
    - args: 'df' (pandas sdataframe with datetime index)
    - kwargs: 'output_freq' (int, minutes) - the output interval required
              'smooth_window' (int, minutes) - applies a centered running 
              average over the requested interval
    """
    
    df.index = pd.to_datetime(df.Datetime)
    df.drop('Datetime', axis = 1, inplace = True)
    
    if not smooth_window == 0:
        current_freq = pd.infer_freq(df.index)
        upsample_df = df.resample('1T').interpolate()
        smooth_df = upsample_df.rolling(window = smooth_window, 
                                        center = True).mean()
        df = smooth_df.resample(current_freq).pad()
        return df.resample('30T').pad()
    else:
        return df.resample('30T').pad()

def calculate_storage(df, site_alt = None):
    
    pass

def sort_names(name_list, return_levels = True):
    
    str_levels_list = [name.split('_')[1][:-1] for name in name_list]
    num_levels_list = []
    for level in str_levels_list:
        try:
            num_levels_list.append(int(level))
        except ValueError:
            num_levels_list.append(float(level))
    sorted_levels_list = [i[1] for i in sorted(enumerate(num_levels_list), 
                          key = lambda x: x[1])]
    sort_idx = [i[0] for i in sorted(enumerate(num_levels_list), 
                key = lambda x: x[1])]
    sorted_name_list = [name_list[i] for i in sort_idx]
    if return_levels:
        return sorted_name_list, sorted_levels_list
    else:
        return sorted_name_list

def get_layers():
    pass
#------------------------------------------------------------------------------

path_to_file = '/home/ian/Documents/profile.csv'
default_press = 101.325
use_Ta_var = None

this_df = pd.read_csv(path_to_file)
df = downsample_data(this_df, smooth_window = 10)

if not 'press' in df.columns:
    if not site_alt is None:
        pass
    else:
        df['press'] = default_press

# Make name and level lists for CO2 (don't assume order will be read correctly)
CO2_names_list = [var for var in df.columns if 'CO2' in var]
CO2_names_list, CO2_levels_list = sort_names(CO2_names_list)

# Make name and level lists for T (don't assume order will be read correctly)
if not use_Ta_var is None:
    T_names_list = [use_Ta_var] * len(CO2_names_list)
else:
    T_names_list = [var for var in df.columns if 'Tair' in var]
    if len(T_names_list) == len(CO2_names_list):
        T_names_list, T_levels_list = sort_names(T_names_list)
        if not CO2_levels_list == T_levels_list: 
            raise Exception('Heights for CO2 and Temperature in input file ' 
                            'do not match!')
    elif len(T_names_list) == 1:
        T_names_list = T_names_list * len(CO2_names_list)
    else:
        raise Exception('Number of temperature variables does not match the '
                        'number of CO2 variables - when this is the case, '
                        'either a single temperature variable should be '
                        'included or the name of the desired temperature '
                        'variable should be passed as a keyword argument '
                        '(use_Ta_var); exiting!')

# Create list of layer depths and layer names
levels_list = CO2_levels_list
zero_levels_list = cp.copy(CO2_levels_list)
zero_levels_list.insert(0, 0)
layer_depths_list = [zero_levels_list[i] - zero_levels_list[i - 1] 
                     for i in range(1, len(zero_levels_list))]
CO2_layer_names_list = ['CO2_{0}-{1}m'.format(str(zero_levels_list[i - 1]),
                                              str(zero_levels_list[i])) 
                        for i in range(1, len(zero_levels_list))]

# Calculate the layer averages (mean of upper and lower boundary);
# lowest layer is just the value observed at the lowest boundary)
layers_df = pd.DataFrame(index = df.index)
for i in range(len(levels_list)):
    if i == 0:
        level_name = CO2_names_list[i]
        layer_name = CO2_layer_names_list[i]
        layers_df[layer_name] = df[level_name]
    else:
        upper_level_name = CO2_names_list[i]
        lower_level_name = CO2_names_list[i - 1]
        layer_name = CO2_layer_names_list[i]
        layers_df[layer_name] = (df[upper_level_name] + 
                                 df[lower_level_name]) / 2


fig, ax = plt.subplots(1, 1, figsize = (12, 8))
fig.patch.set_facecolor('white')
colour_idx = np.linspace(0, 1, len(CO2_layer_names_list))
ax.tick_params(axis = 'x', labelsize = 14)
ax.tick_params(axis = 'y', labelsize = 14)
ax.set_xlabel('$Date$', fontsize = 18)
ax.set_ylabel('$CO2\/\/[ppm]$', fontsize = 18)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
for i, var in enumerate(CO2_layer_names_list):
    color = plt.cm.cool(colour_idx[i])
    plt.plot(layers_df.index, layers_df[var], label = var, color = color)
plt.legend(loc='upper left', frameon = False)



#
#
#
#
#
#
#
#
#
#
#
#






