# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 10:41:38 2015

@author: imchugh
"""

# Python standard modules
import os
import numpy as np

# My modules
import DataIO as io
import data_filtering as filt
import respiration as re
import random_error as ra
import datetime_functions as dtf

# This gets the data
def get_data(configs_dict):

    # Initialise name change dictionary with new names via common keys
    vars_dict = configs_dict['variables']
    newNames_dict = {'carbon_flux':'Fc',
                     'modelled_carbon_flux': 'Fc_model',
                     'temperature': 'TempC',
                     'solar_radiation': 'Fsd',
                     'vapour_pressure_deficit': 'VPD',
                     'friction_velocity': 'ustar',
                     'wind_speed': 'ws'}

    # Get file extension and target
    paths_dict = configs_dict['files']
    ext = os.path.splitext(paths_dict['input_file'])[1]
    data_input_target = os.path.join(paths_dict['input_path'],
                                     paths_dict['input_file'])

    # get data (screen only the Fc data to obs only)
    if ext == '.nc':
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
    elif ext == '.df':
        data_dict, global_attr = io.DINGO_df_to_data_structure(
                                     data_input_target,
                                     var_list = vars_dict.values(),
                                     return_global_attr = True)
        date_time = data_dict.pop('date_time')
    
    # Create continuous model series from ER and Fc series
    vars_dict['modelled_carbon_flux'] = 'Fc_model'
    temp_dict = {'Fc_model': np.concatenate([data_dict['ER_SOLO_all']
                                             [data_dict['Fsd'] < 10],
                                             data_dict['Fc_SOLO']
                                             [data_dict['Fsd'] >= 10]]),
                 'date_time': np.concatenate([data_dict['date_time']
                                             [data_dict['Fsd'] < 10],
                                             data_dict['date_time']
                                             [data_dict['Fsd'] >= 10]])}
    temp_dict = filt.sort_dict_on_index_variable(temp_dict, 'date_time')
    data_dict['Fc_model'] = temp_dict['Fc_model']
                                  
    # Reconstruct data dict with standard names used by algorithms
    new_dict = {}
    for key in vars_dict.keys():
        if key in newNames_dict.keys():
            new_dict[newNames_dict[key]] = data_dict.pop(vars_dict[key])  
    new_dict['date_time'] = date_time

    return new_dict, global_attr
    
#------------------------------------------------------------------------------    

# Get configurations
configs_dict = io.config_to_dict(io.file_select_dialog())

# Get data
data_dict, attr = get_data(configs_dict)

# Assign observational Fc to the 'Fc_series' var
data_dict['Fc_series'] = data_dict['Fc']

# Get the datetime variable so can construct a new partitioned dataset later
datetime_array = data_dict['date_time']

# Make local var names for config items
re_configs_dict = configs_dict['respiration_configs']
re_configs_dict['measurement_interval'] = configs_dict['globals']['measurement_interval']
step = re_configs_dict['step_size_days']
window = re_configs_dict['window_size_days']
num_trials = configs_dict['partitioning_uncertainty']['num_trials']

# Calculate Re by sending data to main respiration function
re_dict, params_dict = re.main(data_dict, configs_dict)
data_dict['Re'] = re_dict['Re']
                                                       
#------------------------------------------------------------------------------
                                                       
# Get the indices of the start and end rows of each unique date in the source 
# data array
dates_input_index_dict = dtf.get_day_indices(datetime_array)

# Get the indices of the start and end rows of each window in the source 
# data array
step_dates_input_index_dict = dtf.get_moving_window_indices(datetime_array, 
                                                            window, step)
                                                            
# Get the indices of the start and end rows of each year in the source 
# data array                                                            
years_input_index_dict = dtf.get_year_indices(datetime_array)

#------------------------------------------------------------------------------

# Get random error estimate using model
ra_configs_dict = configs_dict['random_error_configs']
ra_configs_dict['measurement_interval'] = configs_dict['globals']['measurement_interval']
ra_fig, ra_stats_dict = ra.regress_sigma_delta(data_dict, ra_configs_dict)
sigma_delta = ra.estimate_sigma_delta(data_dict['Re'], ra_stats_dict)

#------------------------------------------------------------------------------

# Initalise parameter dicts with prior estimates
all_noct_dict = re.filtering(data_dict)
params_in_dict = {'Eo_prior': 100,
                  'rb_prior': all_noct_dict['Fc_series'].mean()}

# Make results arrays
annual_re_sums_dict = {year: np.zeros([num_trials]) for year in 
                       years_input_index_dict.keys()}

#------------------------------------------------------------------------------

# Loop through number of trials
for i in xrange(num_trials):
    
    # Generate a results dictionary for the parameter values (1 for each day)
    params_out_dict = re.generate_results_array(datetime_array)

    # Generate noise estimate
    noise = ra.estimate_random_error(sigma_delta)

    # Sum noise with model estimate of Re
    data_dict['Fc_series'] = data_dict['Re'] + noise

    # Partition data into year and step
    years_data_dict = re.segment_data(data_dict, years_input_index_dict)
    step_data_dict = re.segment_data(data_dict, step_dates_input_index_dict)

    # Get Eo for all years
    re.calculate_Eo(years_data_dict, 
                    re_configs_dict,
                    params_in_dict,
                    params_out_dict)
    
    # Get rb for all steps
    re.calculate_rb(step_data_dict,
                    re_configs_dict,
                    params_in_dict,
                    params_out_dict)
    
    # Estimate Re for all data
    this_dict = re.estimate_Re(data_dict,
                               params_out_dict,
                               dates_input_index_dict)

    # Calculate annual sum for each year                     
    for j, year in enumerate(years_input_index_dict.keys()):
        indices = years_input_index_dict[year]
        annual_re_sums_dict[year][i] = (this_dict['Re'][indices[0]: 
                                        indices[1] + 1] * 12 * 0.0018).sum()
                                        
                                        