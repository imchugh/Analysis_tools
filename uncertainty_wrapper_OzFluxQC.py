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
#    vars_dict['date_and_time_stamp'] = 'date_time'
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
    global_configs_dict = configs_master_dict['global_options']
    global_configs_dict['measurement_interval'] = int(attr['time_step'])
    random_configs_dict = dict(global_configs_dict,
                               **configs_master_dict['random_error_options'])
    model_configs_dict = dict(global_configs_dict,
                              **configs_master_dict['model_error_options'])

    # Create copy dataset separated into years
    years_data_dict = data.subset_datayear_from_arraydict(data_dict.copy(), 
                                                          'date_time')
    
    # Create a second copy dataset that is further subdivided into night and day
    night_day_years_data_dict = {}
    for year in years_data_dict:
        subset_dict = (data.subset_onthreshold_from_arraydict
                           (years_data_dict[year],
                           'Fsd',
                           model_configs_dict['noct_threshold']))
        subset_dict = {'day': subset_dict['>'], 'night': subset_dict['<']}
        subset_dict['night'] = (data.set_arraydict_to_nan_conditional
                                    (subset_dict['night'], 
                                    'ustar', 
                                     model_configs_dict['ustar_threshold'][year], 
                                     set_vars = ['Fc']))
        night_day_years_data_dict[year] = subset_dict

    ################
    # RANDOM ERROR #
    ################

    # Calculate the linear regression parameters of sigma_delta as a function 
    # of flux magnitude
    fig, stats_dict = rand_err.regress_sigma_delta(data_dict, random_configs_dict)
    fig.savefig(os.path.join(global_configs_dict['output_path'], 
                             'Random_error_plots.jpg'))
    
    # Calculate estimated sigma_delta for each data point, and remove records 
    # where no observational estimate is available (only has an effect if the 
    # propagation series is a model - which is recommended!!!); then compound
    # random error estimates for all records where obs are available, 
    # then sum and scale appropriately to estimate annual error in gC m-2

    random_error_dict = {}
    prop_series = random_configs_dict['propagation_series']
    for year in years_data_dict.keys():
        sig_del_array = (rand_err.estimate_sigma_delta
                            (years_data_dict[year][prop_series], 
                             stats_dict))
        sig_del_array = sig_del_array[~np.isnan(years_data_dict[year]['Fc'])]
        random_error_dict[year] = (rand_err.propagate_random_error
                                      (sig_del_array,
                                       random_configs_dict))
                                    

    ###############
    # MODEL ERROR #
    ###############

    # Subset data and calculate model error for each year
    model_error_dict = {'day': {}, 'night': {}}
    for year in years_data_dict.keys():
        subset_dict = (data.subset_onthreshold_from_arraydict(
                           years_data_dict[year],
                           'Fsd',
                           model_configs_dict['noct_threshold']))
        subset_dict = {'day': subset_dict['>'], 'night': subset_dict['<']}
        subset_dict['night'] = (data.set_arraydict_to_nan_conditional(
                                    subset_dict['night'], 'ustar', 0.42, 
                                    set_vars = ['Fc']))
        for cond in subset_dict:
            print 'For ' + str(year) + ' ' + cond
            model_error_dict[cond][year] = (mod_err.model_error
                                                (subset_dict[cond], 
                                                 model_configs_dict))    



    ######################
    # PARTITIONING ERROR #
    ######################

#    test_0, test_1 = mod_err.nocturnal_model_error(data_dict, model_configs_dict)
    
#    return test_1 

    # Sum all error
    results_dict = {}
    results_dict['random_error'] = random_error_dict
    results_dict['daytime_model_error'] = model_error_dict['day']
    results_dict['nocturnal_model_error'] = model_error_dict['night']
    results_dict['combined_error'] = {}
    for year in years_data_dict.keys():
        results_dict['combined_error'][year] = np.sqrt(
             (results_dict['random_error'][year])**2 +
             (results_dict['daytime_model_error'][year])**2 + 
             (results_dict['nocturnal_model_error'][year])**2)

    return results_dict    