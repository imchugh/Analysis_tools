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
    """
    Docstring coming soon!
    """    
    def __init__(self, df, CO2_str = 'CO2', Tair_str = 'Tair', use_Tair = None):
        
        self.df = df
        self.use_Tair = use_Tair
        self.CO2_level_names, self.CO2_levels = self.get_levels(CO2_str)
        self.Tair_level_names, self.Tair_levels = self.get_levels(Tair_str)
        self.n_Tair = len(self.Tair_level_names)
        self.Tair_level_names = self.cross_check_Tair()
        self.levels = self.CO2_levels
        self.CO2_layer_names, self.CO2_layers = self.get_layers(CO2_str)

    def get_levels(self, search_str, return_levels = True):

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

    def get_layers(self, layer_str, return_layers = True):
        
        zero_levels_list = cp.copy(self.levels)
        zero_levels_list.insert(0, 0)
        
        names_list = []
        depths_list = []
        for i in range(1, len(zero_levels_list)):
            depths_list.append(zero_levels_list[i] - zero_levels_list[i - 1])
            names_list.append('{0}_{1}-{2}m'.format(layer_str,
                                                    str(zero_levels_list[i - 1]),
                                                    str(zero_levels_list[i])))
        if return_layers:
            return names_list, depths_list
        else:
            return names_list
        
    def cross_check_Tair(self):
        
        if self.use_Tair:
            if isinstance(self.use_Tair, str):
                if self.use_Tair in self.Tair_level_names:
                    return [self.use_Tair] * len(self.CO2_level_names)
                else:
                    raise KeyError('Kwarg \'Tair_var\' not found in dataframe!')
            else:
                raise IOError('Kwarg \'Tair_var\' must be a string!')
        elif len(self.Tair_level_names) == 1:
            return self.Tair_level_names * len(self.CO2_level_names)
        elif len(self.Tair_level_names) == len(self.CO2_level_names):        
            if not self.Tair_levels == self.CO2_levels:
                raise KeyError('Levels for air temperature do not match '
                               'levels for CO2!')
        else:
            raise IOError('Number of air temperature variables in data file '
                          'does not match number of CO2 variables!')
        return self.Tair_level_names
    
def downsample_data(df, output_freq = 30, smooth_window = 0):

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

    local_df = cp.copy(df)
    
    if not smooth_window == 0:
        current_freq = pd.infer_freq(local_df.index)
        upsample_df = local_df.resample('1T').interpolate()
        smooth_df = upsample_df.rolling(window = smooth_window, 
                                        center = True).mean()
        downsample_df = smooth_df.resample(current_freq).pad()
        return downsample_df.resample(output_freq_string).pad()
    else:
        return local_df.resample(output_freq_string).pad()

def get_formatted_data():
    
    file_path = '/home/ian/Documents/profile.csv'
    
    df = pd.read_csv(file_path)
    df.index = pd.to_datetime(df.Datetime)
    return df

def get_layer_means(df, level_names, levels, layer_names):
    
    layers_df = pd.DataFrame(index = df.index)
    for i in range(len(levels)):   
        if i == 0:
            level_name = level_names[i]
            layer_name = layer_names[i]
            layers_df[layer_name] = df[level_name]
        else:
            upper_level_name = level_names[i]
            lower_level_name = level_names[i - 1]
            layer_name = layer_names[i]
            layers_df[layer_name] = (df[upper_level_name] + 
                                     df[lower_level_name]) / 2
    
    return layers_df

def check_ps(df, site_alt):
    if not 'ps' in df.columns:
        df['ps'] = 101.3

def calculate_CO2_storage(df, 
                          CO2_layer_names_list,
                          Tair_layer_names_list,
                          layer_depths_list):
    
    storage_df = pd.DataFrame(index = df.index)

    for i, var in enumerate(CO2_layer_names_list):
        
        molar_density = (df['ps'] * 10**3 / 
                         (8.3143 * (273.15 + df[Tair_layer_names_list[i]])))
        
        molar_density_CO2 = molar_density * df[var] * 10**-6
                                              
        layer_moles_CO2 = molar_density_CO2 * layer_depths_list[i]

        delta_layer_moles_CO2 = (layer_moles_CO2 - layer_moles_CO2.shift()) * 10**6                                              
        delta_layer_moles_CO2_per_sec = delta_layer_moles_CO2 / 1800                                      
        storage_df[storage_names[i]] = delta_layer_moles_CO2_per_sec
        
    storage_df['Sc_total'] = storage_df[storage_names].sum(axis = 1)
        
   
    return storage_df


def main(site_alt = None, use_Tair = None):
    
    raw_data = get_formatted_data()
    
    profile_obj = storage(raw_data, use_Tair = use_Tair)
    
    profile_obj.cross_check_Tair()
        
    downsample_df = downsample_data(profile_obj.df)

    check_ps(downsample_df, site_alt = 100)
    
    layers_df = get_layer_means(downsample_df, 
                                profile_obj.CO2_level_names, 
                                profile_obj.CO2_levels,
                                profile_obj.CO2_layer_names) 

        
    
    return profile_obj



def get_standard_press(elevation):
    
    return 101.3

    

#------------------------------------------------------------------------------


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