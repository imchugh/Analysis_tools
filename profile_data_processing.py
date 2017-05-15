#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 09:52:14 2017

@author: ian
"""

import numpy as np
import pandas as pd
import copy as cp
import os
import pdb
import matplotlib.pyplot as plt
import Tkinter, tkFileDialog

import site_profile_data_processing as spdp
reload(spdp)
#------------------------------------------------------------------------------
class storage(object):
    """
    This class is a container for a set of profile data attributes required 
    for calculation of the storage term, including names and heights of
    CO2 and temperature variables, and layer depths used for vertical scaling 
    of the storage estimates; the only positional argument required is a 
    list of column names (str; note that the column names must...)
    """    
    def __init__(self, col_names, CO2_str = 'CO2', Tair_str = 'Tair', 
                 use_Tair = None):
        
        self.col_names = col_names
        self.CO2_level_names, self.CO2_levels = self.get_levels(CO2_str)
        self.Tair_level_names, self.Tair_levels = self.get_levels(Tair_str)
        self.use_Tair = use_Tair
        if len(self.Tair_level_names) == 1: 
            self.use_Tair = self.Tair_level_names[0]
        self.Tair_level_names = self.cross_check_Tair()
        self.levels = self.CO2_levels
        self.CO2_layer_names, self.CO2_layers = self.get_layers(CO2_str)
        if self.use_Tair is None:
            self.Tair_layer_names, self.Tair_layers = self.get_layers(Tair_str)
        else:
            self.Tair_layer_names = [self.use_Tair] * len(self.CO2_levels)

    def get_levels(self, search_str, return_levels = True):

        name_list = [var for var in self.col_names if search_str in var]        
        
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
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def calculate_CO2_storage(df, profile_obj):
    
    storage_names_list = ['Sc{0}'.format(var.split('CO2')[1]) for var in 
                          profile_obj.CO2_layer_names]
    
    storage_df = pd.DataFrame(index = df.index)

    for i, var in enumerate(profile_obj.CO2_layer_names):

        molar_density = (df['ps'] * 10**3 / 
                         (8.3143 * (273.15 + df[profile_obj.Tair_layer_names[i]])))     
        molar_density_CO2 = molar_density * df[var] * 10**-6                                           
        layer_moles_CO2 = molar_density_CO2 * profile_obj.CO2_layers[i]
        delta_layer_moles_CO2 = ((layer_moles_CO2 - layer_moles_CO2.shift()) 
                                 * 10**6)                                              
        delta_layer_moles_CO2_per_sec = delta_layer_moles_CO2 / 1800                                      
        storage_df[storage_names_list[i]] = delta_layer_moles_CO2_per_sec
        
    storage_df['Sc_total'] = storage_df[storage_names_list].sum(axis = 1)
        
    return storage_df
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def check_ps(df, site_alt):
    if not 'ps' in df.columns:
        if site_alt is None:
            print ('Warning: there are no pressure data available in the raw '
                   'data file and a site altitude has not been specified; '
                   'standard sea level pressure will be used for subsequent '
                   'calculations but may result in substantial storage '
                   'underestimation for high altitude sites (by a factor of '
                   '1-p/p0!')
            df['ps'] = 101.3
        else:
            p0 = 101325
            L = 0.0065
            R = 8.3143
            T0 = 288.15
            g = 9.80665
            M = 0.0289644
            
            A = (g * M) / (R * L)
            B = L / T0
            
            p = (p0 * (1 - B * site_alt) ** A) / 1000
                          
            df['ps'] = p
            
            print p
            
    return
            
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def dir_select_dialog():
    
    """ Open a file select dialog to get path for file retrieval"""
    
    root = Tkinter.Tk(); root.withdraw()
    dir_in = tkFileDialog.askdirectory()
    root.destroy()   
    return dir_in
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------    
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
    current_freq = pd.infer_freq(local_df.index)
    
    if not smooth_window == 0:
        upsample_df = local_df.resample('1T').interpolate()
        smooth_df = upsample_df.rolling(window = smooth_window, 
                                        center = True).mean()
        downsample_df = smooth_df.resample(current_freq).pad()
    else:
        downsample_df = local_df
    
    return downsample_df.resample(output_freq_string).pad()
    
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def file_select_dialog():
    
    """ Open a file select dialog to get path for file retrieval"""
    
    root = Tkinter.Tk(); root.withdraw()
    file_in = tkFileDialog.askopenfilename(initialdir='')
    root.destroy()   
    return file_in
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def get_formatted_data():
    
    file_path = file_select_dialog()
    df = pd.read_csv(file_path)
    df.index = pd.to_datetime(df.Datetime)
    return df
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def get_layer_means(df, level_names, levels, layer_names):
    
    mean_df = pd.DataFrame(index = df.index)
    for i in range(len(levels)):   
        if i == 0:
            level_name = level_names[i]
            layer_name = layer_names[i]
            mean_df[layer_name] = df[level_name]
        else:
            upper_level_name = level_names[i]
            lower_level_name = level_names[i - 1]
            layer_name = layer_names[i]
            mean_df[layer_name] = (df[upper_level_name] + 
                                     df[lower_level_name]) / 2
    
    return mean_df
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def get_raw_data(site, output_dir = False):
    
    df = spdp.get_site_data(site)
    if output_dir:
        file_name = site + '_OzFlux_format_profile.csv'
        write_data_to_file(df, output_dir, file_name)
    return df

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def make_layers_df(df, profile_obj):
    
    layers_df = get_layer_means(df,  
                                profile_obj.CO2_level_names, 
                                profile_obj.CO2_levels,
                                profile_obj.CO2_layer_names)

    if profile_obj.use_Tair == None:
        layers_df = layers_df.join(get_layer_means(df,
                                                   profile_obj.Tair_level_names, 
                                                   profile_obj.Tair_levels,
                                                   profile_obj.Tair_layer_names))
    else:
        layers_df[profile_obj.use_Tair] = df[profile_obj.use_Tair]
        
    layers_df['ps'] = df['ps']
    
    return layers_df
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def plot_diurnal(df):

    # No idea why this is needed
    for var in df.columns:
        df[var] = df[var].astype(np.float64)

    diurnal_df = df.groupby([lambda x: x.hour, lambda y: y.minute]).mean()
    diurnal_df.index = np.arange(48) / 2.0
 
    vars_list = list(df.columns)
    vars_list.remove('Sc_total')
    strip_vars_list = [var.split('_')[1] for var in vars_list]
    fig, ax = plt.subplots(1, 1, figsize = (12, 8))
    fig.patch.set_facecolor('white')
    colour_idx = np.linspace(0, 1, len(vars_list))
    ax.set_xlim([0, 24])
    ax.set_xticks([0,4,8,12,16,20,24])
    ax.tick_params(axis = 'x', labelsize = 14)
    ax.tick_params(axis = 'y', labelsize = 14)
    ax.set_xlabel('$Time$', fontsize = 18)
    ax.set_ylabel('$S_c\/(\mu mol\/CO_2\/m^{-2}\/s^{-1})$', fontsize = 18)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    for i, var in enumerate(vars_list):
        color = plt.cm.cool(colour_idx[i])
        plt.plot(diurnal_df.index, diurnal_df[var], 
                 label = strip_vars_list[i], color = color)
    plt.plot(diurnal_df.index, diurnal_df.Sc_total, 
             label = 'Total', color = 'grey')
    plt.legend(loc=[0.65, 0.18], frameon = False, ncol = 2)    
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def plot_time_series(df):
    
    vars_list = list(df.columns)
    vars_list.remove('Sc_total')
    strip_vars_list = [var.split('_')[1] for var in vars_list]
    fig, ax = plt.subplots(1, 1, figsize = (12, 8))
    fig.patch.set_facecolor('white')
    colour_idx = np.linspace(0, 1, len(vars_list))
    ax.tick_params(axis = 'x', labelsize = 14)
    ax.tick_params(axis = 'y', labelsize = 14)
    ax.set_xlabel('$Date$', fontsize = 18)
    ax.set_ylabel('$S_c\/(\mu mol\/CO_2\/m^{-2}\/s^{-1})$', fontsize = 18)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    for i, var in enumerate(vars_list):
        color = plt.cm.cool(colour_idx[i])
        plt.plot(df.index, df[var], label = strip_vars_list[i], color = color)
    plt.plot(df.index, df.Sc_total, label = 'Total', color = 'grey')
    plt.legend(loc='lower left', frameon = False, ncol = 2)    
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def write_data_to_file(df, path, file_name):
    
    if not os.path.isdir(path):
        home_dir = os.path.expanduser('~')
        output_dir = os.path.join(home_dir, 'profile_data')
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        print ('Warning - the specified directory for data output could not be'
               ' found! \nWriting to the following directory: {0}'
               .format(output_dir))
    else:
        output_dir = path

    target = os.path.join(output_dir, file_name)    
     
    df.to_csv(target, index_label = 'Datetime')
    
    return
#------------------------------------------------------------------------------

###############################################################################
# Start main program                                                          #
###############################################################################   
def main(site = None, site_alt = None, use_Tair = None, 
         output_freq = 30, plot_ts = False, plot_diurnal_avg = False, 
         output_dir = None):

    """
    Processes profile CO2 data to storage estimates;
    
    Kwargs:
        - site (str; default None): if a site name is passed, the script 
          searches for the site raw data processing script in the site-specific
          profile module (profile_data_processing.py); if the site string is 
          valid, the user must then select the valid file or directory in the 
          file select dialog; if the site string is invalid, an IOError will be 
          raised 
        - site_alt (int; default None): script searches for a pressure variable 
          named 'ps' in the input file; if present, it uses the data from that 
          variable; if absent, it calculates the mean expected pressure for the 
          given site altitude; if site_alt is None AND there is no pressure 
          variable in the dataset, mean sea level pressure is used as input 
          to the equation of state; note that this may lead to large errors for
          high-altitude sites (a user warning to this effect is printed to 
          screen)
        - use_Tair (str: default None): if a temperature variable name is 
          passed, this variable will be used as the temperature input to the 
          equation of state (which converts mole fraction to molar density); 
          otherwise, the script searches for either a set of temperature 
          variables with heights matching those for the CO2 variables, or a 
          single temperature variable; if there are multiple temperature 
          variables in the file but they do not match the heights for the CO2
          variables, or there are more than 1 but less than the number of CO2 
          variables, an exception is raised
        - output_freq (int; default 30 [units - minutes]): the desired output 
          frequency of the storage calculation
        - plot_ts (bool; default False): output a plot of the entire time 
          series for storage
        - plot_diurnal_avg (bool; default False): output a plot of the 
          diurnally averaged data
        - output_dir (str; default None): if not None, writes data to a csv 
          file in the nominated directory (if the directory is invalid, a
          warning message is printed to screen and the data is written to the 
          location documented in the message)
    """
    
    # Either prompt for already formatted file, or prompt for directory 
    # containing unformatted file/s
    if site is None:
        data_df = get_formatted_data()
    else:
        data_df = get_raw_data(site)

    profile_obj = storage(data_df.columns, use_Tair = use_Tair)
    
    downsample_df = downsample_data(data_df, output_freq = output_freq)

    check_ps(downsample_df, site_alt = site_alt)
    
    layers_df = make_layers_df(downsample_df, profile_obj)
        
    storage_df = calculate_CO2_storage(layers_df, profile_obj)
    
    if plot_ts:
        plot_time_series(storage_df)
    if plot_diurnal_avg:
        plot_diurnal(storage_df)
    
    if not output_dir is None:
        if site is None: site = 'generic'
        file_name = site + '_storage.csv'
        write_data_to_file(storage_df, output_dir, file_name)
    
    return storage_df
###############################################################################