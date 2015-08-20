# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 16:20:41 2015

@author: imchugh
"""

# Standard modules
import numpy as np
import os
import pdb

# My modules
import DataIO as io
import random_error as rand_err
import model_error as mod_err

def get_data(configs_dict):

    # Initialise name change dictionary with new names via common keys
    vars_dict = configs_dict['variables']
    newNames_dict = {'carbon_flux':'Fc',
                     'modelled_carbon_flux': 'Fc_model',
                     'solar_radiation':'Fsd',
                     'temperature':'Ta',
                     'wind_speed': 'ws',
                     'friction_velocity': 'ustar'}

    # Get data (screen only the Fc data to obs only)
    data_input_target = os.path.join(configs_dict['files']['data_path'],
                                     configs_dict['files']['data_file'])
    Fc_dict = io.OzFluxQCnc_to_data_structure(data_input_target,
                                              var_list = [vars_dict['carbon_flux']],
                                              QC_accept_codes = [0])

    date_time = Fc_dict.pop('date_time')
    ancillary_vars = [vars_dict[var] for var in vars_dict.keys() 
                      if not var == 'carbon_flux']
    ancillary_dict, global_attr = io.OzFluxQCnc_to_data_structure(
                                      data_input_target,
                                      var_list = ancillary_vars,
                                      return_global_attr = True)
    data_dict = dict(Fc_dict, **ancillary_dict)                                                     

    # Reconstruct data dict with standard names used by algorithms
    new_dict = {}
    for key in vars_dict.keys():
        new_dict[newNames_dict[key]] = data_dict.pop(vars_dict[key])  
    new_dict['date_time'] = date_time
    return new_dict, global_attr

def main():    

    # Update
    reload(rand_err)
    reload(mod_err)
    reload(io)

    # Unpack configs and open data file
    configs_master_dict = io.config_to_dict(io.file_select_dialog())
    data_dict, attr = get_data(configs_master_dict)
    global_configs_dict = configs_master_dict['global_options']
    global_configs_dict['measurement_interval'] = int(attr['time_step'])
    random_configs_dict = dict(global_configs_dict,
                               **configs_master_dict['random_error_options'])
    model_configs_dict = dict(global_configs_dict,
                              **configs_master_dict['model_error_options'])

    # Get list of years to iterate over
    years_array = np.array([date_.year for date_ in data_dict['date_time']])
    years_set_list = list(set(years_array))
    
    test_0, test_1 = mod_err.nocturnal_model_error(data_dict, model_configs_dict)
    
    return test_1 
    
#    # Subset data and calculate daytime model error for each year
#    model_error_dict = {'day': {}, 'night': {}}
#    for year in years_set_list:
#        year_index = years_array == year
#        year_data_dict = {var: data_dict[var][year_index] 
#                          for var in ['Fc', 'Fc_model', 'Fsd']}
#        model_error_dict['day'][str(year)] = (mod_err.daytime_model_error
#                                              (year_data_dict, 
#                                               model_configs_dict))
#    
#        
#    
#    # Calculate the linear regression parameters of sigma_delta as a function 
#    # of flux magnitude
#    fig, stats_dict = rand_err.regress_sigma_delta(data_dict, random_configs_dict)
#    fig.savefig(os.path.join(global_configs_dict['output_path'], 
#                             'Random_error_plots.jpg'))
#    
#    # Calculate estimated sigma_delta for each data point, then remove records 
#    # where no observational estimate is available (only has an effect if the 
#    # propagation series is a model - which is recommended!!!)
#    sig_del_array = (rand_err.estimate_sigma_delta
#                        (data_dict[random_configs_dict['propagation_series']], 
#                         stats_dict))
#    sig_del_array = sig_del_array[~np.isnan(data_dict['Fc'])]
#    obs_years_array = years_array[~np.isnan(data_dict['Fc'])]
#    
#    # Calculate random error estimates for each record where obs are available, 
#    # then sum and scale appropriately to estimate annual error in gC m-2
#    random_error_dict = {}
#    for year in years_set_list:
#        year_sig_del_array = sig_del_array[obs_years_array == year]
#        random_error_dict[str(year)] = (rand_err.propagate_random_error
#                                            (year_sig_del_array,
#                                             random_configs_dict))
#
#    # Sum all error
#    results_dict = {}
#    results_dict['random_error'] = random_error_dict
#    results_dict['daytime_model_error'] = model_error_dict['day']
#    results_dict['combined_error'] = {}
#    for year in years_set_list:
#        results_dict['combined_error'][str(year)] = np.sqrt(
#             (results_dict['random_error'][str(year)])**2 +
#             (results_dict['daytime_model_error'][str(year)])**2)

    return results_dict    