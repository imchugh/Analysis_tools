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
import data_filtering as data_filter
import random_error as rand_err
import model_error as mod_err
import data_formatting as dt_fm

reload(dt_fm)
#------------------------------------------------------------------------------
# Fetch data from configurations
def get_data(configs_dict):

    # Get data (screen Fc data to obs only - keep gap-filled drivers etc)
    data_input_target = os.path.join(configs_dict['files']['input_path'],
                                     configs_dict['files']['input_file'])
    Fc_dict = io.OzFluxQCnc_to_data_structure(data_input_target,
                                              var_list = [configs_dict['variables']
                                                                      ['carbon_flux']],
                                              QC_accept_codes = [0])
    Fc_dict.pop('date_time')
    ancillary_vars = [configs_dict['variables'][var] for var in 
                      configs_dict['variables'] if not var == 'carbon_flux']
    ancillary_dict, attr = io.OzFluxQCnc_to_data_structure(
                               data_input_target,
                               var_list = ancillary_vars,
                               return_global_attr = True)
    data_dict = dict(Fc_dict, **ancillary_dict)
    
    # Rename to generic names used by scripts
    old_names_dict = configs_dict['variables']
    std_names_dict = dt_fm.standard_names_dictionary()
    map_names_dict = {old_names_dict[key]: std_names_dict[key] 
                      for key in old_names_dict}
    data_dict = dt_fm.rename_data_dict_vars(data_dict, map_names_dict)
    
    return data_dict, attr
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Screen out low ustar and separate into day and night
def filter_data(data_dict, configs_dict):
    
    ustar_threshold = configs_dict['ustar_threshold']
    noct_threshold = configs_dict['noct_threshold']
    data_dict['Fc'][(data_dict['Fsd'] < noct_threshold) & 
                    (data_dict['ustar'] < ustar_threshold)] = np.nan
    subset_dict = {}
    subset_dict['day'] = data_filter.subset_arraydict_on_threshold(
                             data_dict, 'Fsd', noct_threshold, '>', drop = True)
    subset_dict['night'] = data_filter.subset_arraydict_on_threshold(
                               data_dict, 'Fsd', noct_threshold, '<', drop = True)                                 
    
    return subset_dict
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def main():    

    # Update
    reload(rand_err)
    reload(mod_err)
    reload(io)
    reload(data_filter)

    # Unpack master config file and rebuild single config file
    configs_master_dict = io.config_to_dict(io.file_select_dialog())
    configs_dict = {'variables': dict(configs_master_dict['random_error_configs']['variables'],
                                      **configs_master_dict['model_error_configs']['variables']),
                    'model_options': configs_master_dict['model_error_configs']['options'],
                    'random_options': configs_master_dict['random_error_configs']['options'],
                    'global_options': configs_master_dict['global_configs']['options'],
                    'files': configs_master_dict['global_configs']['files'],
                    'uncertainty_options': configs_master_dict['NEE_uncertainty_configs']} 

    # Get the data
    data_dict, attr = get_data(configs_dict)

    # Sum Fc and Sc if storage is to be included, otherwise if requested, 
    # remove all Fc where Sc is missing
    if configs_dict['global_options']['use_storage']:
        data_dict['NEE_series'] = (data_dict['NEE_series'] + 
                                   data_dict['Fc_storage'])
    elif configs_dict['global_options']['unify_flux_storage_cases']:
        data_dict['NEE_series'][np.isnan(data_dict['Fc_storage'])] = np.nan

    # Save the time step information into the configuration file
    configs_dict['random_options']['measurement_interval'] = int(attr['time_step'])
    configs_dict['model_options']['measurement_interval'] = int(attr['time_step'])

    # Calculate the linear regression parameters of sigma_delta as a function 
    # of flux magnitude
    fig, stats_dict = rand_err.regress_sigma_delta(data_dict, configs_dict['random_options'])
    fig.savefig(os.path.join(configs_dict['files']['output_path'], 
                             'Random_error_plots.jpg')) 

    # Calculate estimated sigma_delta for each data point, and remove records 
    # where no observational estimate is available (only has an effect if the 
    # propagation series is a model - which is recommended!!!);                             
    sig_del_array = (rand_err.estimate_sigma_delta
                        (data_dict[configs_dict['random_options']['propagation_series']], 
                         stats_dict))
    sig_del_array[np.isnan(data_dict['NEE_series'])] = np.nan
    data_dict['sigma_delta'] = sig_del_array 

    # Create copy dataset separated into years
    years_data_dict = data_filter.subset_datayear_from_arraydict(data_dict, 
                                                                'date_time')   

    return

    # do the uncertainty analysis for each year
    num_trials = configs_dict['uncertainty_options']['num_trials']
    for this_year in years_data_dict.keys():

        # Generate an array of ustar values based on mu and sigma from change 
        # point detection analysis
        if configs_dict['uncertainty_options']['do_ustar_uncertainty']:
            ustar_threshold = configs_dict['global_options']['ustar_threshold'][str(this_year)]
            ustar_uncertainty = configs_dict['global_options']['ustar_uncertainty'][str(this_year)]
            ustar_array = np.random.normal(loc = ustar_threshold,
                                           scale = ustar_uncertainty,
                                           size = num_trials)

        # Set a filter flag to prevent repeat filtering for same ustar threshold 
        # when ustar_uncertainty is set to false
        filter_flag = False

        for this_trial in xrange(configs_dict['uncertainty_options']['num_trials']):

            # If including ustar uncertainty, set ustar threshold, then filter,
            # then generate model estimates
            if configs_dict['uncertainty_options']['do_ustar_uncertainty']:
                this_config_dict = {'ustar_threshold': ustar_array[this_trial],
                                    'noct_threshold': (configs_dict['Global_options']
                                                           ['noct_threshold'])}
                this_dict = filter_data(years_data_dict[this_year].copy(), 
                                        this_config_dict)
            # ... otherwise just do the filtering with the best estimate ustar
            # but only do it once (flip the filter flag on the first pass)!
            else:
                if not filter_flag:
                    this_config_dict = {'ustar_threshold': (configs_dict
                                                            ['global_options']
                                                            ['ustar_threshold']
                                                            [this_year]),
                                        'noct_threshold': (configs_dict
                                                           ['Global_options']
                                                           ['noct_threshold'])}
                    this_dict = filter_data(years_data_dict[this_year].copy(), 
                                            configs_dict)
                    filter_flag = True

    return ustar_array

    # Calculate errors
    error_dict = {}        
    for year in years_data_dict.keys():
        error_dict[year] = {}
        total_obs = 0
        for cond in years_data_dict[year].keys():
            
            # Grab data            
            temp_dict = years_data_dict[year][cond]
            
            # Do availability stats
            stats_dict = data_filter.count_nans_in_array(temp_dict['Fc'])
            total_obs = total_obs + stats_dict['Total_obs']
            error_dict[year]['Avail_obs_' + cond] = stats_dict['Avail_obs']
            error_dict[year]['Pct_avail_obs_' + cond] = stats_dict['Pct_avail_obs']
            
            # Do random error            
            error_dict[year]['random_error_' + cond] = (
                rand_err.propagate_random_error
                    (temp_dict['sigma_delta'][~np.isnan(temp_dict['sigma_delta'])],
                     configs_dict['random_options']))

            # Do model error
            error_dict[year]['model_error_' + cond] = (
                mod_err.model_error(temp_dict, configs_dict['model_options']))  

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
                   os.path.join(configs_dict['files']['output_path'], 
                                'Uncertainty_output.csv'))
    
    return error_dict
#------------------------------------------------------------------------------