# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 16:20:41 2015

@author: imchugh
"""

import numpy as np
import os
import pdb
import netCDF4
from scipy import stats
import datetime as dt

import DataIO as io
import random_error as rand_err
import model_error as mod_err

# Open configuration and build dictionaries of config file contents
def get_random_error_configs():
    
    configs_dict = {'measurement_interval': 30,
                    'neg_averaging_bins': 10,
                    'pos_averaging_bins': 10,
                    'radiation_difference_threshold': 35,
                    'temperature_difference_threshold': 3,
                    'windspeed_difference_threshold': 1,
                    'noct_threshold': 5,
                    'ustar_threshold': 0,
                    'num_trials': 10**4,
                    'nan_value': -9999,
                    'QC_accept_code': 0,
                    'mean_flux_series': 'Fc_model',
                    'propagation_series': 'Fc_model'}
                    
    return configs_dict

def get_data(configs_dict):

    # Specify file location and name
    data_input_path = '/home/imchugh/Ozflux/Sites/Whroo/Data/Processed/all'
    data_input_file = 'Whroo_2011_to_2014_L5.nc'
    data_input_target = os.path.join(data_input_path, data_input_file)

    # Initialise dictionaries
    vars_dict = {'carbon_flux':'Fc',
                 'modelled_carbon_flux': 'Fc_SOLO',
                 'solar_radiation':'Fsd',
                 'temperature':'Ta',
                 'wind_speed': 'Ws_CSAT',
                 'friction_velocity': 'ustar'}
    newNames_dict = {'carbon_flux':'Fc',
                     'modelled_carbon_flux': 'Fc_model',
                     'solar_radiation':'Fsd',
                     'temperature':'Ta',
                     'wind_speed': 'ws',
                     'friction_velocity': 'ustar'}

    data_dict, attr = io.OzFluxQCnc_to_data_structure(data_input_target,
                                                      var_list = vars_dict.values(),
                                                      QC_accept_codes = [0])

    new_dict = {}
    
    for key in vars_dict.keys():
        new_dict[newNames_dict[key]] = data_dict.pop(vars_dict[key])  
    return new_dict

def main():    

    # User options
    QC_accept_code = 0
    propagation_series = 'Fc_model'
    

    # Specify file locations and names
    data_input_path = '/home/imchugh/Ozflux/Sites/Whroo/Data/Processed/all'
    data_input_file = 'Whroo_2011_to_2014_L5.nc'
    data_input_target = os.path.join(data_input_path, data_input_file)

    configs_input_path = '/home/imchugh/Code/Python/Config_files/'
    random_error_configs_input_file = 'random_error_config_new.txt'
    configs_input_target = os.path.join(configs_input_path, 
                                        random_error_configs_input_file)
    
    results_out_path = '/home/imchugh/Documents'

    # Update
    reload(rand_err)
    reload(mod_err)

    # Get configurations and data
    configs_dict = io.config_to_dict(configs_input_target)
    data_dict = get_data(configs_dict)
    return data_dict
    # Get list of years to iterate over
    years_array = np.array([date_.year for date_ in data_dict['date_time']])
    years_set_list = list(set(years_array))
    
    # Subset data and calculate daytime model error for each year
    model_error_dict = {'day': {}, 'night': {}}
    for year in years_set_list:
        year_index = years_array == year
        year_data_dict = {var: data_dict[var][year_index] 
                          for var in ['Fc', 'Fc_model', 'Fsd']}
        model_error_dict['day'][str(year)] = (mod_err.daytime_model_error
                                              (year_data_dict, configs_dict))
    
    # Calculate the linear regression parameters of sigma_delta as a function 
    # of flux magnitude
    fig, stats_dict = rand_err.regress_sigma_delta(data_dict, configs_dict)    
    fig.savefig(os.path.join(results_out_path, 'Random_error_plots.jpg'))
    
    # Calculate estimated sigma_delta for each data point, then remove records 
    # where no observational estimate is available (only has an effect if the 
    # propagation series is a model - which is recommended!!!)
    sig_del_array = (rand_err.estimate_sigma_delta
                        (data_dict[configs_dict['propagation_series']], 
                         stats_dict))
    sig_del_array = sig_del_array[~np.isnan(data_dict['Fc'])]
    obs_years_array = years_array[~np.isnan(data_dict['Fc'])]
    
    # Calculate random error estimates for each record where obs are available, 
    # then sum and scale appropriately to estimate annual error in gC m-2
    random_error_dict = {}
    for year in years_set_list:
        year_sig_del_array = sig_del_array[obs_years_array == year]
        print len(year_sig_del_array)
        random_error_dict[str(year)] = (rand_err.propagate_random_error
                                            (year_sig_del_array,
                                             configs_dict))

    return random_error_dict
