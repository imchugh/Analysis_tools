# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 16:20:41 2015

@author: imchugh
"""

# Standard modules
import numpy as np
import os
import copy as cp
import calendar
import sys
import pdb

# My modules
import DataIO as io
import data_filtering as filt
import random_error as rand_err
import model_error as mod_err
import data_formatting as dt_fm
import respiration as re
import photosynthesis as ps

reload(dt_fm)
reload(filt)
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
# Rebuild the master configuration file for passing to respiration and light 
# response (if requested) functions
def build_config_file(configs_master_dict):

    # Build a specific configuration file
    configs_dict = {'files': configs_master_dict['global_configs']['files'],
                    'global_options': (configs_master_dict['global_configs']
                                                          ['options']),
                    'uncertainty_options': configs_master_dict['NEE_uncertainty_configs']}                                                          
    configs_dict['variables'] = {}
    dict_list = [configs_master_dict['random_error_configs']['variables'],
                 configs_master_dict['model_error_configs']['variables'],
                 configs_master_dict['respiration_configs']['variables'],
                 configs_master_dict['photosynthesis_configs']['variables']]
    for d in dict_list:
        configs_dict['variables'].update(d)
    return configs_dict                                                         
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Set all sigma delta values to nan where there are no observational data
def filter_sigma_delta(data_dict):
    data_dict['sigma_delta'][np.isnan(data_dict['NEE_series'])] = np.nan    
#------------------------------------------------------------------------------    

#------------------------------------------------------------------------------
# Screen out low ustar
def filter_ustar(data_dict, ustar_threshold, noct_threshold):
    data_dict['NEE_series'][(data_dict['Fsd'] < noct_threshold) & 
                            (data_dict['ustar'] < ustar_threshold)] = np.nan
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Split into day and night
def separate_night_day(data_dict, noct_threshold):
    subset_dict = {}
    subset_dict['day'] = filt.subset_arraydict_on_threshold(
                             data_dict, 'Fsd', noct_threshold, '>', drop = True)
    subset_dict['night'] = filt.subset_arraydict_on_threshold(
                               data_dict, 'Fsd', noct_threshold, '<', drop = True)
    return subset_dict
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Make a results dictionary
def init_interm_rslt_dict(num_trials):
    var_list = ['random_error_day', 'model_error_day', 'obs_avail_day',  
                'random_error_night', 'model_error_night', 'obs_avail_night',
                'ustar_error']
    nan_array = np.zeros(num_trials)
    nan_array[:] = np.nan
    return {var: cp.copy(nan_array) for var in var_list} 
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def init_final_rslt_dict(years_data_dict, configs_dict):

    years = years_data_dict.keys()
    final_dict = {this_year: {} for this_year in years}
    for this_year in years:
        days = 366 if calendar.isleap(this_year) else 365
        records = days * 1440 / configs_dict['measurement_interval']
        final_dict[this_year]['total_records'] = records
    return final_dict
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def run_model(data_dict, re_configs_dict, ps_configs_dict):
    
    re_rslt_dict, re_params_dict = re.main(data_dict, re_configs_dict)[0: 2]
    ps_rslt_dict = ps.main(data_dict, ps_configs_dict, re_params_dict)[0]
    data_dict['NEE_model'] = ps_rslt_dict['GPP'] + ps_rslt_dict['Re']
    data_dict['NEE_filled'] = np.where(np.isnan(data_dict['NEE_series']),
                                       data_dict['NEE_model'],
                                       data_dict['NEE_series'])
    return    
#------------------------------------------------------------------------------

def main():    

    # Update
    reload(rand_err)
    reload(mod_err)
    reload(io)
    reload(filt)

    #---------------------------
    # Preparation and formatting
    #---------------------------

    # Get master config file
    configs_master_dict = io.config_to_dict(io.file_select_dialog())

    # Build custom configuration file for this script
    configs_dict = build_config_file(configs_master_dict)

    # Get the data
    data_dict, attr = get_data(configs_dict)

    # Build required configuration files for imported scripts (random error,
    # model error, respiration, light response)
    rand_err_configs_dict = configs_master_dict['random_error_configs']['options']
    mod_err_configs_dict = configs_master_dict['model_error_configs']['options']
    re_configs_dict = configs_master_dict['respiration_configs']['options']
    ps_configs_dict = configs_master_dict['photosynthesis_configs']['options']

    # Save the time step information into the individual configuration files
    for d in [configs_dict, rand_err_configs_dict, mod_err_configs_dict, 
              re_configs_dict, ps_configs_dict]: 
        d['measurement_interval'] = int(attr['time_step'])
    
    # For respiration and light response, turn off the output option for 
    # window fits of the functions, even if requested in the configuration file
    # - they are WAY too computationally expensive!!!
    if re_configs_dict['output_fit_plots']:
        re_configs_dict['output_fit_plots'] = False 
    if ps_configs_dict['output_fit_plots']:
        ps_configs_dict['output_fit_plots'] = False 
    
    # Sum Fc and Sc if storage is to be included, otherwise if requested, 
    # remove all Fc where Sc is missing
    if configs_dict['global_options']['use_storage']:
        data_dict['NEE_series'] = (data_dict['NEE_series'] + 
                                   data_dict['Sc'])
    elif configs_dict['global_options']['unify_flux_storage_cases']:
        data_dict['NEE_series'][np.isnan(data_dict['Sc'])] = np.nan

    # Convert insolation to PPFD for light response calculations
    data_dict['PAR'] = data_dict['Fsd'] * 0.46 * 4.6       

    #---------------------------------------------
    # Initial model estimation for random error
    # (note: low u* data is left in intentionally)
    #---------------------------------------------

    # Generate initial model series for Re and GPP then combine
    run_model(data_dict, re_configs_dict, ps_configs_dict)

    #----------------------------------------
    # Random error calculation and statistics
    #----------------------------------------

    # Calculate the linear regression parameters of sigma_delta as a function 
    # of flux magnitude
    fig, stats_dict, rslt_dict = rand_err.regress_sigma_delta(
                                     data_dict, rand_err_configs_dict)
    fig.savefig(os.path.join(configs_dict['files']['output_path'], 
                             'Random_error_plots.jpg'))

    # Calculate estimated sigma_delta for each data point, and remove records 
    # where no observational estimate is available (only has an effect if the 
    # propagation series is a model - which is recommended!!!);                             
    sig_del_array = (rand_err.estimate_sigma_delta
                        (data_dict[rand_err_configs_dict['propagation_series']], 
                         stats_dict))
    data_dict['sigma_delta'] = sig_del_array

    #---------------------
    # Uncertainty analysis
    #---------------------

    print '\nRunning uncertainty analysis'
    print '----------------------------\n'

    # Write universal config items to local variables
    num_trials = configs_dict['uncertainty_options']['num_trials']
    noct_threshold = configs_dict['global_options']['noct_threshold']
    ustar_threshold = configs_dict['global_options']['ustar_threshold']
    do_ustar_uncertainty = configs_dict['uncertainty_options']['do_ustar_uncertainty']
    do_random_uncertainty = configs_dict['uncertainty_options']['do_random_uncertainty']
    do_model_uncertainty = configs_dict['uncertainty_options']['do_model_uncertainty']

#    # Generate a final results dictionary
#    final_rslt_dict = init_final_rslt_dict(years_data_dict, configs_dict)

#    # Generate a standard data dictionary and calculate annual NEE for all years
#    # using best estimate of ustar
#    this_dict = cp.deepcopy(data_dict)
#    filt.screen_low_ustar(this_dict, {'noct_threshold': noct_threshold,
#                                      'ustar_threshold': ustar_threshold})
#    print ('This much after return to main function: ' + 
#            str(len(data_dict['NEE_series'][~np.isnan(data_dict['NEE_series'])])))
#    filter_sigma_delta(this_dict)
#    run_model(this_dict, re_configs_dict, ps_configs_dict)
#    years_data_dict = filt.subset_datayear_from_arraydict(this_dict, 
#                                                          'date_time')
#    annual_NEE_sum_dict = {}
#    for this_year in years_data_dict.keys():
#        annual_NEE_sum_dict[this_year] = (years_data_dict[this_year]['NEE_filled'] * 
#                                          configs_dict['measurement_interval'] *
#                                          60 * 12 * 10**-6).sum()
#
#    # 
#    for this_trial in xrange(num_trials):
#
#        # If including ustar uncertainty, generate a u* estimate for each year
#
#        continue

    # Create dataset separated into years
    years_data_dict = filt.subset_datayear_from_arraydict(data_dict, 
                                                          'date_time')   

    final_rslt_dict = {}
        
    # Do the uncertainty analysis for each year        
    for this_year in years_data_dict.keys():

        print 'Running analysis for ' + str(this_year) + ':'

        # Make an intermediate results dictionary
        interm_rslt_dict = init_interm_rslt_dict(num_trials)
        
        # Write ustar thresholds for years to local variables
        if isinstance(configs_dict['global_options']['ustar_threshold'],
                      dict):
            ustar_threshold = (configs_dict['global_options']
                                           ['ustar_threshold']
                                           [str(this_year)])
        else:
            ustar_threshold = configs_dict['global_options']['ustar_threshold']
        
        # Write ustar uncertainty for years to local variables
        if isinstance(configs_dict['global_options']['ustar_uncertainty'],
                      dict):
            ustar_uncertainty = (configs_dict['global_options']
                                             ['ustar_uncertainty']
                                             [str(this_year)])
        else:
            ustar_uncertainty = (configs_dict['global_options']
                                             ['ustar_uncertainty'])

        # Generate a standard data dictionary; this will be overwritten if
        # ustar uncertainty is set to True, but the reference value for NEE
        # will be retained
        this_dict = cp.deepcopy(years_data_dict[this_year])
        filt.screen_low_ustar(this_dict, {'noct_threshold': noct_threshold,
                                          'ustar_threshold': ustar_threshold})
        try:
            run_model(this_dict, re_configs_dict, ps_configs_dict)
        except:
            print ('    - Excluding the year ' + str(this_year) + 
                   ' - insufficient data!')
            continue # Do the next year
        NEE_sum = (this_dict['NEE_filled'] * 
                   configs_dict['measurement_interval'] *
                   60 * 12 * 10**-6).sum()

        # If including ustar uncertainty:
        #   1) generate an array of ustar values based 
        #      on mu and sigma from change point detection analysis
        #   2) add the resulting array to the intermediate results dict
        #   3) add an empty array to keep NEE error due to ustar 
        if do_ustar_uncertainty:
            ustar_array = np.random.normal(loc = ustar_threshold,
                                           scale = ustar_uncertainty,
                                           size = num_trials)
            interm_rslt_dict['u_star'] = ustar_array
        else:
            interm_rslt_dict['u_star'] = np.tile(ustar_threshold, num_trials)
           
        # If not including ustar uncertainty, just count the available n 
        # for day and night and write to the intermediate results
        if not do_ustar_uncertainty:
            for cond in this_dict.keys():
                interm_rslt_dict['obs_avail_' + cond][:] = (
                    this_dict[cond]['NEE_series']
                        [~np.isnan(this_dict[cond]['NEE_series'])])
            
        # Do trials
        for this_trial in xrange(num_trials):

            # Print progress
            if this_trial == 0:
                print '    - Trial: ' + str(this_trial + 1),
            elif this_trial == num_trials - 1:
                print str(this_trial + 1) + ' ... Done!'
            else:
                print this_trial + 1,

            # If including ustar uncertainty:
            #   1) make a deep copy of the original data so it doesn't get overwritten
            #   2) set ustar threshold 
            #   3) filter for low ustar
            #   4) gap fill the filtered dataset
            #   5) sum, calculate difference relative to best u* and output to dict
            #   6) separate into night and day, ready for random and model error
            if do_ustar_uncertainty:
                this_dict = cp.deepcopy(years_data_dict[this_year])
                ustar_threshold = ustar_array[this_trial]
                filt.screen_low_ustar(this_dict, {'noct_threshold': noct_threshold,
                                                  'ustar_threshold': ustar_threshold})
                try:
                    run_model(this_dict, re_configs_dict, ps_configs_dict)
                except:
                    next
                this_sum = (this_dict['NEE_filled'] * 
                            configs_dict['measurement_interval'] * 60 *
                            12 * 10**-6).sum()
                interm_rslt_dict['ustar_error'][this_trial] = NEE_sum - this_sum
            
            # Before random and model error estimation:
            #   1) screen all sigma_delta values where observations are missing
            #   2) split into day and night
            filter_sigma_delta(this_dict)
            this_dict = separate_night_day(this_dict, noct_threshold)            

            # For each of day and night
            for cond in this_dict.keys():

                # If including ustar uncertainty, write the available n to the 
                # intermediate results dictionary for day and night
                if do_ustar_uncertainty:
                    interm_rslt_dict['obs_avail_' + cond][this_trial] = len(
                        this_dict[cond]['NEE_series']
                            [~np.isnan(this_dict[cond]['NEE_series'])])

                # Do the random error and write to correct position in 
                # intermediate results dict
                if do_random_uncertainty:
                    sig_del_array = (this_dict[cond]['sigma_delta']
                                     [~np.isnan(this_dict[cond]['sigma_delta'])])
                    error_array = rand_err.estimate_random_error(sig_del_array)                
                    interm_rslt_dict['random_error_' + cond][this_trial] = (
                        error_array.sum() * configs_dict['measurement_interval'] 
                                          * 60 * 12 * 10 ** -6)

                # Do the model error and write to correct position in 
                # intermediate results dict
                if do_model_uncertainty:
                    sub_dict = cp.deepcopy(this_dict[cond])
                    interm_rslt_dict['model_error_' + cond][this_trial] = (
                        mod_err.estimate_model_error(sub_dict, 
                                                     mod_err_configs_dict))

        final_rslt_dict[this_year] = interm_rslt_dict

    return final_rslt_dict

#        # Do statistics
#        lst = [var for var in interm_rslt_dict.keys() if 'error' in var]
#        this_dict = {var: np.round(interm_rslt_dict[var].std() * 2, 1) for 
#                     var in lst}
#        total_error_float = np.sqrt((np.array(this_dict.values())**2).sum())
#        this_dict['total_error'] = np.round(total_error_float, 2)
#        final_rslt_dict[this_year].update(this_dict)
#        
##        test_array = (interm_rslt_dict['random_error_day'] + interm_rslt_dict['random_error_night'] + 
##                      interm_rslt_dict['model_error_day'])
##        print 'This is calc' + str(round(test_array.std() * 2, 2))
#        
#    return final_rslt_dict, interm_rslt_dict