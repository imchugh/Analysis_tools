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

class storage(object):
    
    def __init__(self, target_file, CO2_str = 'CO2', Ta_str = 'Tair'):
        
        self.target_file = target_file
        self.df = self.get_data()
        self.CO2_level_names, self.CO2_levels = self.get_level_names(CO2_str)
        self.T_channels = [var for var in self.df.columns if Ta_str in var]
#        self.levels = [int(this_label.split('_')[1][:-1]) for 
#                       this_label in self.CO2_channels]

    def get_data(self):
    
        return pd.read_csv(self.target_file)

    def downsample_data(self, output_freq = 30, smooth_window = 0):
    
        """
        This function downsamples profile data to the requested output frequency
        (generally 30 minutes);
        - args: 'df' (pandas sdataframe with datetime index)
        - kwargs: 'output_freq' (int, minutes) - the output interval required
                  'smooth_window' (int, minutes) - applies a centered running 
                  average over the requested interval
        """
        
        if not isinstance(output_freq, int):
            raise Exception('Output frequency must be an integer (units = minutes)')
        output_freq_string = '{0}T'.format(str(output_freq))

        local_df = cp.copy(self.df)

        local_df.index = pd.to_datetime(local_df.Datetime)
        local_df.drop('Datetime', axis = 1, inplace = True)
        
        if not smooth_window == 0:
            current_freq = pd.infer_freq(local_df.index)
            upsample_df = local_df.resample('1T').interpolate()
            smooth_df = upsample_df.rolling(window = smooth_window, 
                                            center = True).mean()
            downsample_df = smooth_df.resample(current_freq).pad()
            return downsample_df.resample(output_freq_string).pad()
        else:
            return local_df.resample(output_freq_string).pad()

    def get_level_names(self, search_str, return_levels = True):

        name_list = [var for var in self.df.columns if search_str in var]        
        
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

def calculate_CO2_storage(df, site_alt = None):
    
    CO2_layer_names_list = get_names(df)
    Tair_layer_names_list = get_names(df, var = 'Tair')
    if not len(Tair_layer_names_list) == len(CO2_layer_names_list):
        Tair_layer_names_list = Tair_layer_names_list * len(CO2_layer_names_list)
    lookup_dict = dict(zip(CO2_layer_names_list, Tair_layer_names_list))
    
    storage_names = ['Sc{0}'.format(this_var.split('CO2')[1]) 
                     for this_var in CO2_layer_names_list]
    
    storage_depths = get_layer_depths_from_layer_names(CO2_layer_names_list)
    
    storage_df = pd.DataFrame(index = df.index)
    
    if not 'ps' in df.columns:
        df['ps'] = get_standard_press(site_alt)
    
    for i, var in enumerate(CO2_layer_names_list):
        
        molar_density = (df['ps'] * 10**3 / 
                         (8.3143 * (273.15 + df[Tair_layer_names_list[i]])))
        
        molar_density_CO2 = molar_density * df[var] * 10**-6
                                              
        layer_moles_CO2 = molar_density_CO2 * storage_depths[i]

        delta_layer_moles_CO2 = (layer_moles_CO2 - layer_moles_CO2.shift()) * 10**6                                              
        delta_layer_moles_CO2_per_sec = delta_layer_moles_CO2 / 1800                                      
        storage_df[storage_names[i]] = delta_layer_moles_CO2_per_sec
        
    storage_df['Sc_total'] = storage_df[storage_names].sum(axis = 1)
        
   
    return storage_df

def get_standard_press(elevation):
    
    return 101.3


def get_layer_names(levels_list, var = 'CO2'):
    
    zero_levels_list = cp.copy(levels_list)
    zero_levels_list.insert(0, 0)
    return ['{0}_{1}-{2}m'.format(var,
                                  str(zero_levels_list[i - 1]),
                                  str(zero_levels_list[i])) 
            for i in range(1, len(zero_levels_list))]

def get_layer_depths_from_layer_names(layer_names_list):
    
    new_list = [name.replace('-', '_').split('_') for name in layer_names_list]
    lower_plane = [int(i[1]) for i in new_list]    
    upper_plane = [int(i[2][:-1]) for i in new_list]
    return [upper_plane[i] - lower_plane[i] for i in range(len(upper_plane))]

def get_layer_depths_from_levels(levels_list):
    
    zero_levels_list = cp.copy(levels_list)
    zero_levels_list.insert(0, 0)
    return [zero_levels_list[i] - zero_levels_list[i - 1] 
            for i in range(1, len(zero_levels_list))]
    
def get_layer_means(df, var = 'CO2'):
    
    level_names_list = get_names(df, var = var)    
    level_names_list, levels_list = sort_level_names(level_names_list)
    layer_names_list = get_layer_names(levels_list, var = var)
    
    layers_df = pd.DataFrame(index = df.index)
    for i in range(len(levels_list)):   
        if i == 0:
            level_name = level_names_list[i]
            layer_name = layer_names_list[i]
            layers_df[layer_name] = df[level_name]
        else:
            upper_level_name = level_names_list[i]
            lower_level_name = level_names_list[i - 1]
            layer_name = layer_names_list[i]
            layers_df[layer_name] = (df[upper_level_name] + 
                                     df[lower_level_name]) / 2
    
    return layers_df

#------------------------------------------------------------------------------

path_to_file = '/home/ian/Documents/profile.csv'
default_press = 101.325
site_alt = None
output_int = 30
smooth_window = 10
use_Tair_var = False

#this_df = pd.read_csv(path_to_file)
#df = downsample_data(this_df, smooth_window = 10)
#
#if not 'press' in df.columns:
#    if not site_alt is None:
#        pass
#    else:
#        df['press'] = default_press
#
## Make name and level lists for CO2
#CO2_level_names_list = [var for var in df.columns if 'CO2' in var]
#CO2_level_names_list, CO2_levels_list = sort_level_names(CO2_level_names_list)
#
## Make name and level lists for T (include various checks)
#Tair_level_names_list = [var for var in df.columns if 'Tair' in var]
#if use_Tair_var:
#    if not use_Tair_var in Tair_level_names_list:
#        raise KeyError('User-specified air temperature variable not found in '
#                        'data file - exiting!')
#    n_Tair = 1
#else:
#    if len(Tair_level_names_list) == 1:
#        use_Tair_var == Tair_level_names_list[0]
#        n_Tair = 1
#    if len(Tair_level_names_list) == len(CO2_level_names_list):
#        Tair_level_names_list, Tair_levels_list = sort_level_names(Tair_level_names_list)
#        n_Tair = len(CO2_level_names_list)
#    else:
#        raise Exception('Number of temperature variables does not match the '
#                        'number of CO2 variables - when this is the case, '
#                        'either a single temperature variable should be '
#                        'included or the name of the desired temperature '
#                        'variable should be passed as a keyword argument '
#                        '(use_Ta_var); exiting!') 
#
#CO2_layers_df = get_layer_means(df)
#if n_Tair == 1:
#    layers_df = CO2_layers_df
#    layers_df['Tair'] = df[use_Tair_var]
#else:    
#    Tair_layers_df = get_layer_means(df, var = 'Tair')
#    layers_df = CO2_layers_df.join(Tair_layers_df)
#
#CO2_layer_names_list = get_layer_names(CO2_levels_list)
#
#test = calculate_CO2_storage(layers_df)
#
#diurnal_df = test.groupby([lambda x: x.hour, lambda y: y.minute]).mean()
#diurnal_df.index = np.arange(48) / 2.0


#fig, ax = plt.subplots(1, 1, figsize = (12, 8))
#fig.patch.set_facecolor('white')
#colour_idx = np.linspace(0, 1, len(CO2_layer_names_list))
#ax.tick_params(axis = 'x', labelsize = 14)
#ax.tick_params(axis = 'y', labelsize = 14)
#ax.set_xlabel('$Date$', fontsize = 18)
#ax.set_ylabel('$CO2\/\/[ppm]$', fontsize = 18)
#ax.xaxis.set_ticks_position('bottom')
#ax.yaxis.set_ticks_position('left')
#ax.spines['right'].set_visible(False)
#ax.spines['top'].set_visible(False)
#for i, var in enumerate(CO2_layer_names_list):
#    color = plt.cm.cool(colour_idx[i])
#    plt.plot(layers_df.index, layers_df[var], label = var, color = color)
#plt.legend(loc='upper left', frameon = False)

#fig, ax = plt.subplots(1, 1, figsize = (12, 8))
#fig.patch.set_facecolor('white')
#colour_idx = np.linspace(0, 1, len(CO2_layer_names_list))
#ax.tick_params(axis = 'x', labelsize = 14)
#ax.tick_params(axis = 'y', labelsize = 14)
#ax.set_xlabel('$Date$', fontsize = 18)
#ax.set_ylabel('$CO2\/\/[ppm]$', fontsize = 18)
#ax.xaxis.set_ticks_position('bottom')
#ax.yaxis.set_ticks_position('left')
#ax.spines['right'].set_visible(False)
#ax.spines['top'].set_visible(False)
#for i, var in enumerate(CO2_layer_names_list):
#    color = plt.cm.cool(colour_idx[i])
#    plt.plot(layers_df.index, layers_df[var], label = var, color = color)
#plt.legend(loc='upper left', frameon = False)