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
import data_handling as data
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

    # Get file extension and target
    ext = os.path.splitext(configs_dict['files']['data_file'])[1]
    data_input_target = os.path.join(configs_dict['files']['data_path'],
                                     configs_dict['files']['data_file'])

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
    reload(data)

    # Unpack configs and open data file
    configs_master_dict = io.config_to_dict(io.file_select_dialog())
    data_dict, attr = get_data(configs_master_dict)
    global_configs_dict = configs_master_dict['global_configs']
    global_configs_dict['measurement_interval'] = int(attr['time_step'])
    random_configs_dict = dict(global_configs_dict,
                               **configs_master_dict['random_error_configs'])
    model_configs_dict = dict(global_configs_dict,
                              **configs_master_dict['model_error_configs'])

    # Calculate the linear regression parameters of sigma_delta as a function 
    # of flux magnitude
    fig, stats_dict = rand_err.regress_sigma_delta(data_dict, random_configs_dict)
    fig.savefig(os.path.join(global_configs_dict['output_path'], 
                             'Random_error_plots.jpg')) 
                             
    # Calculate estimated sigma_delta for each data point, and remove records 
    # where no observational estimate is available (only has an effect if the 
    # propagation series is a model - which is recommended!!!);                             
    sig_del_array = (rand_err.estimate_sigma_delta
                        (data_dict[random_configs_dict['propagation_series']], 
                         stats_dict))
    sig_del_array[np.isnan(data_dict['Fc'])] = np.nan
    data_dict['sigma_delta'] = sig_del_array 

    # Create copy dataset separated into years, screen out low ustar and 
    # separate into day and night, ready for uncertainty calculations
    years_data_dict = data.subset_datayear_from_arraydict(data_dict.copy(), 
                                                          'date_time')                                                          
    for year in years_data_dict.keys():
        ustar_threshold = global_configs_dict['ustar_threshold'][str(year)]
        noct_threshold = global_configs_dict['noct_threshold']
        years_data_dict[year]['Fc'][(years_data_dict[year]['Fsd'] < 
                                     noct_threshold) & 
                                    (years_data_dict[year]['ustar'] < 
                                     ustar_threshold)] = np.nan
        subset_dict = (data.subset_onthreshold_from_arraydict(
                           years_data_dict[year],
                           'Fsd',
                           global_configs_dict['noct_threshold']))
        subset_dict = {'day': subset_dict['>'], 'night': subset_dict['<']}
        years_data_dict[year] = subset_dict

    # Calculate errors
    error_dict = {}        
    for year in years_data_dict.keys():
        error_dict[year] = {}
        total_obs = 0
        for cond in years_data_dict[year].keys():
            
            # Grab data            
            temp_dict = years_data_dict[year][cond]
            
            # Do availability stats
            stats_dict = data.count_nans_in_array(temp_dict['Fc'])
            total_obs = total_obs + stats_dict['Total_obs']
            error_dict[year]['Avail_obs_' + cond] = stats_dict['Avail_obs']
            error_dict[year]['Pct_avail_obs_' + cond] = stats_dict['Pct_avail_obs']
            
            # Do random error            
            error_dict[year]['random_error_' + cond] = (
                rand_err.propagate_random_error
                    (temp_dict['sigma_delta'][~np.isnan(temp_dict['sigma_delta'])],
                     random_configs_dict))

            # Do model error
            error_dict[year]['model_error_' + cond] = (
                mod_err.model_error(temp_dict, model_configs_dict))  

        # Do combined stats
        error_dict[year]['Year'] = year
        error_dict[year]['Avail_obs_total'] = (
            error_dict[year]['Avail_obs_day'] +
            error_dict[year]['Avail_obs_night'])
        error_dict[year]['Pct_avail_obs_total'] = (
            round(error_dict[year]['Avail_obs_total'] / float(total_obs) * 100, 1))
        error_dict[year]['combined_error'] = np.round(np.sqrt(
            error_dict[year]['model_error_day']**2 +
            error_dict[year]['model_error_night']**2 +
            error_dict[year]['random_error_day']**2 +
            error_dict[year]['random_error_night']**2), 2)

    # Ouput
    io.dict_to_csv(error_dict, 
                   ['Year', 'Avail_obs_day', 'Avail_obs_night', 
                    'Avail_obs_total', 'Pct_avail_obs_day', 
                    'Pct_avail_obs_night', 'Pct_avail_obs_total', 
                    'model_error_day', 'model_error_night', 'random_error_day', 
                    'random_error_night', 'combined_error'],
                   os.path.join(global_configs_dict['output_path'], 
                                'Uncertainty_output.csv'))
    
    return error_dict