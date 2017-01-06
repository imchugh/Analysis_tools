# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 16:20:41 2015

@author: imchugh
"""

# Standard modules
import numpy as np
import os
import copy as cp
from scipy import stats
import matplotlib.pyplot as plt
from scipy.stats import norm
import matplotlib.gridspec as gridspec
import warnings
import logging
import pdb
import datetime as dt

# My modules
import DataIO as io
import data_filtering as filt
import datetime_functions as dtf
import random_error as rand_err
import model_error as mod_err
import data_formatting as dt_fm
import respiration as re
import photosynthesis as ps
import gap_filling as gf

#------------------------------------------------------------------------------
# Fetch data from configurations
def get_data(configs_dict):

    # Get data (screen Fc data to obs only - keep gap-filled drivers etc)
    data_input_target = os.path.join(configs_dict['files']['input_path'],
                                     configs_dict['files']['input_file'])

    ext = os.path.splitext(data_input_target)[1]
    if ext == '.nc':
        Fc_dict = io.OzFluxQCnc_to_data_structure(data_input_target,
                                                  var_list = [configs_dict['variables']
                                                                          ['carbon_flux']],
                                                  QC_accept_code = 0)
        Fc_dict.pop('date_time')
        ancillary_vars = [configs_dict['variables'][var] for var in 
                          configs_dict['variables'] if not var == 'carbon_flux']
        ancillary_dict, attr = io.OzFluxQCnc_to_data_structure(
                                   data_input_target,
                                   var_list = ancillary_vars,
                                   return_global_attr = True)
        data_dict = dict(Fc_dict, **ancillary_dict)

    elif ext == '.df':
        data_dict, attr = io.DINGO_df_to_data_structure(data_input_target,
                              var_list = configs_dict['variables'].values(),
                              return_global_attr = True)

    # Rename to generic names used by scripts
    old_names_dict = configs_dict['variables']
    std_names_dict = dt_fm.get_standard_names()
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
# Check whether ustar and Fsd (both used for filtering of the dataset) contain
# missing data when NEE data is available - if so, exclude these cases from 
# analysis
def check_data_consistency(data_dict):
    
    warnings.simplefilter('always')
    for var in ['Fsd', 'ustar']:
        
        flag_index = np.where(~np.isnan(data_dict['NEE_series']) & 
                              np.isnan(data_dict[var]))
        count_str = str(len(flag_index[0]))
        if not len(flag_index[0]) == 0:
            warnings.warn('There are %s' %count_str + ' instances where NEE ' 
                          'element contains valid data and %s' %var + ' element ' 
                          'is not a number - NEE values for these ' 
                          'instances will be excluded from analysis!')
        data_dict['NEE_series'][flag_index] = np.nan

#------------------------------------------------------------------------------
# Check whether all model drivers are complete - if not, warn user (may expand 
# this to force exception, since results will be nan)
def check_driver_consistency(data_dict):
    
    warnings.simplefilter('always')
    arr = np.array([])
    for var in ['Fsd', 'TempC', 'VPD']:
        
        flag_index = np.where(np.isnan(data_dict['NEE_series']) & 
                              np.isnan(data_dict[var]))
        arr = np.concatenate([arr, flag_index[0]])
        count_str = str(len(flag_index[0]))                              
        if not len(flag_index[0]) == 0:
            warnings.warn('There are %s' %count_str + ' instances where neither ' \
                          'the NEE nor model driver %s' %var + ' element ' \
                          'contains valid data - model estimates cannot be ' \
                          'calculated for these instances!')
        data_dict['NEE_series'][flag_index] = np.nan        
    arr = np.unique(arr)
    if not len(arr) == 0:
        print 'Total number of instances with missing driver data is ' + str(len(arr))
#------------------------------------------------------------------------------
        
#------------------------------------------------------------------------------
# Set all sigma delta values to nan where there are no observational data
def filter_sigma_delta(data_dict):
    data_dict['sigma_delta'][np.isnan(data_dict['NEE_series'])] = np.nan
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
def init_interm_rslt_dict(num_trials, do_ustar, do_random, do_model):

    var_list = ['obs_avail_day', 'obs_avail_night']
    if do_ustar:
        var_list = var_list + ['u_star', 'ustar_error']
    if do_random:
        var_list = var_list + ['random_error_day', 'random_error_night']
    if do_model:
        var_list = var_list + ['model_error_day', 'model_error_night']

    nan_array = np.zeros(num_trials)
    nan_array[:] = np.nan
    return {var: cp.copy(nan_array) for var in var_list} 
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def run_model(data_dict, NEE_model, re_configs_dict, ps_configs_dict):
    
    if NEE_model == 'LT':    
        
        try:
            re_rslt_dict, re_params_dict = re.main(data_dict, 
                                                   re_configs_dict)[0: 2]
        except Exception:
            raise
        try:
            ps_rslt_dict = ps.main(data_dict, ps_configs_dict, re_params_dict)[0]    
        except Exception:
            raise
        data_dict['NEE_model'] = ps_rslt_dict['GPP'] + ps_rslt_dict['Re']
        data_dict['NEE_filled'] = np.where(np.isnan(data_dict['NEE_series']),
                                           data_dict['NEE_model'],
                                           data_dict['NEE_series'])
                                           
    elif NEE_model == 'ANN':
        
        len_int = len(data_dict['NEE_series'])
        input_array = np.empty([len_int, 4])
        for i, var in enumerate(['TempC', 'Sws', 'Fsd', 'VPD']):
            input_array[:, i] = data_dict[var]
        target_array = np.empty([len_int, 1])
        target_array[:, 0] = data_dict['NEE_series']
        
        data_dict['NEE_model'] = gf.train_ANN(input_array, target_array, 
                                              100, 
                                              [4, 24, 16, 1])[:, 0]
        data_dict['NEE_filled'] = np.where(np.isnan(data_dict['NEE_series']),
                                           data_dict['NEE_model'],
                                           data_dict['NEE_series'])                                                     

    else:
        
        raise Exception('\'' + NEE_model + '\' is not a valid model type! ' \
                        'Valid choices are \'ANN\' or \'LT\'')
                                           
    return    
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def plot_data(data_d):
    """
    Pass a dictionary containing the following key / value pairs:
    keys: year(int or str) / values: dictionary containing following 
    key / value pairs: 
    keys: 'ustar', 'ustar_error', 'model_error_day', 'model_error_night', 
    'random_error_day' and 'random_error_night' as keys / values: equal length 
    numpy arrays (may contain np.nan - will be filtered)
    """
#    
#    key_id_list = [isinstance(key, int) for key in data_d.keys()]    
#    if not all(key_id_list):
#        raise Exception('Expected integer years as outer dictionary key!')    
    
    error_list = list(set([key.split('_')[0] for 
                           key in data_d.keys() if 'error' in key]))

    results_d = {}
    
    if 'ustar' in error_list:
        results_d['ustar'] = data_d['ustar_error']
    
    if 'random' in error_list:
        results_d['random'] = (data_d['random_error_day'] + 
                               data_d['random_error_night'])

    if 'model' in error_list:
        results_d['model'] = (data_d['model_error_day'] + 
                              data_d['model_error_night'])

    if len(error_list) > 1: 
        results_d['total'] = np.zeros(len(results_d['ustar']))
        for var in error_list:
            results_d['total'] = results_d['total'] + results_d[var]        
        bool_array = ~np.isnan(results_d['total'])
        error_list.append('total')
        for var in error_list:
            results_d[var] = results_d[var][bool_array]


    results_d['random'] = results_d['random'] + results_d['total'].mean()
    results_d['model'] = results_d['model'] + results_d['total'].mean()

    colors_d = {'total': 'grey',
                'ustar': 'blue',
                'random': 'cyan',
                'model': 'magenta'}

    pos_d = {'total': 0.9,
             'ustar': 0,
             'random': 0.3,
             'model': 0.6}

    # Do the stats
    mu_dict = {}
    sig_dict = {}
    for var in error_list:
        mu_dict[var] = results_d[var].mean()
        sig_dict[var] = results_d[var].std()

    # Create the plot
    fig = plt.figure(figsize = (12, 10))
    fig.patch.set_facecolor('white')
    gs = gridspec.GridSpec(2, 1, height_ratios=[4,1.5])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])

    # Set up the first subplot
    ax1.set_xlabel('$Uncertainty\/(g\/C\/m^{-2}a^{-1})$',
                  fontsize=18)
    ax1.set_ylabel('$Frequency$', fontsize=18)
    ax1.tick_params(axis = 'y', labelsize = 14)
    ax1.tick_params(axis = 'x', labelsize = 14)
    ax1.axvline(mu_dict['total'], color = 'black', 
                linewidth = 2, linestyle = '--')
    ax1.xaxis.set_ticks_position('bottom')
    ax1.yaxis.set_ticks_position('left')
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    
    # Plot the histogram
    for var in error_list:
        print var
        if var == 'total':
            ec = 'none'
            fc = colors_d[var]
            htype = 'bar'
            al = 0.5
        else:
            ec = colors_d[var]
            fc = 'none'
            htype = 'step'
            al = 1
        ax1.hist(results_d[var], 50, facecolor = fc, edgecolor = ec,
                 orientation = 'vertical', label = var, 
                 histtype = htype, normed = True, alpha = al)
    ax1.legend(loc='upper right', frameon = False)

    # Plot the normal distribution
    xmin, xmax = ax1.get_xlim()
    x = np.linspace(xmin, xmax, 100)
    p = stats.norm.pdf(x, mu_dict['total'], sig_dict['total'])
    ax1.plot(x, p, color = 'black')

    # Set up the second plot
    ax2.axes.get_yaxis().set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.xaxis.set_ticks_position('bottom')
    ax2.set_xticklabels([])                    
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_ylim([0, 1])
    the_buffer = 0.12
    
    # Plot the confidence intervals
    for var in error_list:
        if var == 'ustar': continue
        ax2.plot((mu_dict[var] - sig_dict[var] * 2, 
                  mu_dict[var] + sig_dict[var] * 2), 
                 (pos_d[var], pos_d[var]), color = colors_d[var], linewidth = 2)
        ax2.plot(mu_dict[var] - sig_dict[var] * 2, pos_d[var], 
                 marker = '|', color = colors_d[var], markersize = 10, mew = 2)
        ax2.plot(mu_dict[var] + sig_dict[var] * 2, pos_d[var], 
                 marker = '|', color = colors_d[var], markersize = 10,
                 mew = 2)
        ax2.plot(mu_dict[var], pos_d[var], 
                 marker = 'o', color = colors_d[var], markersize = 10,
                 mec = 'none')
        ax2.text(mu_dict[var] - sig_dict[var] * 2, pos_d[var] - the_buffer, 
                 str(round(mu_dict[var] - sig_dict[var] * 2, 1)),
                 verticalalignment = 'center',
                 horizontalalignment = 'center',
                 fontsize = 14)
        ax2.text(mu_dict[var] + sig_dict[var] * 2, pos_d[var] - the_buffer, 
                 str(round(mu_dict[var] + sig_dict[var] * 2, 1)),
                 verticalalignment = 'center',
                 horizontalalignment = 'center',
                 fontsize = 14)
        if var == 'total':
            ax2.text(mu_dict[var], pos_d[var] - the_buffer, 
                     str(round(mu_dict[var], 1)),
                     verticalalignment = 'center',
                     horizontalalignment = 'center',
                     fontsize = 14)
            

    plt.show()
    
    return
            
def main(output_trial_results = True, output_plot = True):    

    # Update
    reload(rand_err)
    reload(mod_err)
    reload(io)
    reload(filt)
    reload(re)
    reload(gf)
    reload(dt_fm)

    #-----------------------------------
    # General preparation and formatting
    #-----------------------------------

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
    
    # Write universal config items to local variables
    noct_threshold = configs_dict['global_options']['noct_threshold']
    ustar_threshold = configs_dict['global_options']['ustar_threshold']
    num_trials = configs_dict['uncertainty_options']['num_trials']
    do_ustar_uncertainty = configs_dict['uncertainty_options']['do_ustar_uncertainty']
    do_random_uncertainty = configs_dict['uncertainty_options']['do_random_uncertainty']
    do_model_uncertainty = configs_dict['uncertainty_options']['do_model_uncertainty']
    NEE_model = configs_dict['uncertainty_options']['NEE_model']
    measurement_interval = configs_dict['measurement_interval']
    if do_ustar_uncertainty: ustar_uncertainty = (configs_dict['global_options']
                                                              ['ustar_uncertainty'])

    # Print stuff
    print '---------------------------------'
    print 'Running uncertainty analysis for:'
    error_list = ['ustar', 'random', 'model']
    mode_count = 0
    for i, var in enumerate([do_ustar_uncertainty, do_random_uncertainty, 
                             do_model_uncertainty]):
         if var:
            mode_count = mode_count + 1
            print '- ' + error_list[i] + ' error'
    if mode_count == 0:
        raise Exception('Processing flags for all uncertainty sources ' \
                        'currently set to False: set at least one ' \
                        'uncertainty source to True in configuration file ' \
                        'before proceeding!')
    print '---------------------------------'

    # Open log file
    logf = open('/home/imchugh/Documents/log.txt', 'w')

    
    
    #-----------------
    # Data preparation
    #-----------------

    # Sum Fc and Sc if storage is to be included, otherwise if requested, 
    # remove all Fc where Sc is missing
    if configs_dict['global_options']['use_storage']:
        data_dict['NEE_series'] = (data_dict['NEE_series'] + 
                                   data_dict['Sc'])
    elif configs_dict['global_options']['unify_flux_storage_cases']:
        data_dict['NEE_series'][np.isnan(data_dict['Sc'])] = np.nan

    # Convert insolation to PPFD for light response calculations
    data_dict['PAR'] = data_dict['Fsd'] * 0.46 * 4.6       

    # Check no NEE values with missing ustar values
    check_data_consistency(data_dict)

    # Check no drivers missing where NEE is missing
    check_driver_consistency(data_dict)

    #----------------------------------------
    # Random error calculation and statistics
    #----------------------------------------

    if do_random_uncertainty:

        NEE_temp = cp.deepcopy(data_dict['NEE_series'])
        
        # Generate initial model series for Re and GPP then combine
        # (note: low u* data is left in intentionally)
        run_model(data_dict, NEE_model, re_configs_dict, ps_configs_dict)
    
        # Calculate the linear regression parameters of sigma_delta as a function 
        # of flux magnitude
        fig, stats_dict, rslt_dict = rand_err.regress_sigma_delta(
                                         data_dict, rand_err_configs_dict)
    
        # Calculate estimated sigma_delta for each data point, and remove records 
        # where no observational estimate is available (only has an effect if the 
        # propagation series is a model - which is recommended!!!);                             
        sig_del_array = (rand_err.estimate_sigma_delta
                            (data_dict[rand_err_configs_dict['propagation_series']], 
                             stats_dict))
        data_dict['sigma_delta'] = sig_del_array
        
        data_dict['NEE_series'] = NEE_temp
    
    filt.screen_low_ustar(data_dict, ustar_threshold, noct_threshold, configs_dict['global_options']['ustar_filter_day'])
    run_model(data_dict, NEE_model, re_configs_dict, ps_configs_dict)

    #---------------------
    # Uncertainty analysis
    #---------------------

    # Extract a list of years from the dataset
    years_array = np.array([date_.year for date_ in data_dict['date_time']])
    years_list = list(set(years_array))

    # Create a results dictionary
    final_rslt_dict = {this_year: init_interm_rslt_dict(num_trials,
                                                        do_ustar_uncertainty,
                                                        do_random_uncertainty,
                                                        do_model_uncertainty)
                       for this_year in years_list}

    # Get the t-statistic for the 95% CI
    t_stat = stats.t.ppf(0.975, num_trials)

    print '----------------------------'
    print 'Starting uncertainty trials:',

    # Do trials
    for this_trial in xrange(num_trials):

        # Set first_pass flag to prevent repetitive assignment
        first_pass = True if this_trial == 0 else False

        # Print progress
        if not this_trial == num_trials - 1:
            print str(this_trial + 1),
        else:
            print str(this_trial + 1) + ' ... Done!'

        # If doing ustar uncertainty, make a ustar dictionary by randomly 
        # sampling from the normal distribution and scaling according to
        # ustar threshold uncertainty
        if do_ustar_uncertainty:
            z_score = np.random.normal(loc = 0, scale = 1, size = 1)[0]
            if isinstance(ustar_threshold, dict):
                if not isinstance(ustar_uncertainty, dict):
                    raise Exception('ustar_threshold and ustar_uncertainty ' \
                                    'must be specified as same object type ' \
                                    'in configuration file! Exiting...')
                this_ustar = {}
                for key in ustar_threshold.keys():
                    this_ustar[key] = (ustar_threshold[key] +
                                       ustar_uncertainty[key] / 2 * 
                                       z_score)
            elif isinstance(ustar_threshold, int):
                if not isinstance(ustar_uncertainty, int):
                    raise Exception('ustar_threshold and ustar_uncertainty ' \
                                    'must be specified as same object type ' \
                                    'in configuration file! Exiting...')
                this_ustar = ustar_threshold + ustar_uncertainty * z_score
            else:
                raise Exception('ustar_threshold variable must be specified ' \
                                'as type either dict or int in configuration ' \
                                'file... exiting')
        
            # Make a copy of the data dictionary
            this_dict = cp.deepcopy(data_dict)
            
            # Screen low ustar then model and gap fill
            filt.screen_low_ustar(this_dict, this_ustar, noct_threshold, True)
            try:
                run_model(this_dict, NEE_model, re_configs_dict, ps_configs_dict)
            except Exception, e:
                logf.write('Model optimisation for trial {0} failed with the '
                           ' following message: ' + e[0] + '\n')
                continue # Do the next trial

            # If doing random uncertainty, screen out any estimates of 
            # sigma_delta where there are no obs
            if do_random_uncertainty: filter_sigma_delta(this_dict)
            
        # If not doing ustar uncertainty, just assign this_dict to data_dict
        elif first_pass:
                this_dict = data_dict
                filt.screen_low_ustar(this_dict, ustar_threshold, noct_threshold)
                if do_random_uncertainty: filter_sigma_delta(this_dict)
                
        # Create dataset separated into years
        years_data_dict = dtf.subset_datayear_from_arraydict(this_dict, 
                                                             'date_time')
        
        # Do calculations for each year
        for this_year in years_list:
            
            if do_ustar_uncertainty:
                final_rslt_dict[this_year]['u_star'][this_trial] = (
                    this_ustar[str(this_year)])
                NEE_annual_sum = (years_data_dict[this_year]['NEE_filled'] * 
                                  measurement_interval * 60 * 12 * 10**-6).sum()
                final_rslt_dict[this_year]['ustar_error'][this_trial] = NEE_annual_sum

            # Split the dictionary into day and night and calculate the uncertainty for each
            split_dict = separate_night_day(years_data_dict[this_year], noct_threshold)

            # For each of day and night
            for cond in split_dict.keys():

                # If including ustar uncertainty, write the available n to the 
                # intermediate results dictionary for day and night; otherwise,
                # just write it once (since available n won't vary if ustar
                # doesn't)
                if do_ustar_uncertainty:
                    final_rslt_dict[this_year]['obs_avail_' + cond][this_trial] = (
                        len(split_dict[cond]['NEE_series']
                                [~np.isnan(split_dict[cond]['NEE_series'])]))
                else:
                    if first_pass:
                        final_rslt_dict[this_year]['obs_avail_' + cond][:] = (
                            len(split_dict[cond]['NEE_series']
                                    [~np.isnan(split_dict[cond]['NEE_series'])]))
                    else:
                        first_pass = False
                    
                # Do the random error and write to correct position in 
                # intermediate results dict
                if do_random_uncertainty:
                    sig_del_array = (split_dict[cond]['sigma_delta']
                                     [~np.isnan(split_dict[cond]['sigma_delta'])])
                    error_array = rand_err.estimate_random_error(sig_del_array)
                    random_annual_sum = (error_array.sum() * 
                                         measurement_interval
                                         * 60 * 12 * 10 ** -6)
                    final_rslt_dict[this_year]['random_error_' + cond][this_trial] = (
                        random_annual_sum)

                # Do the model error and write to correct position in 
                # intermediate results dict
                if do_model_uncertainty:
                    sub_dict = cp.deepcopy(split_dict[cond])
                    model_annual_sum = mod_err.estimate_model_error(
                                           sub_dict, mod_err_configs_dict)
                    final_rslt_dict[this_year]['model_error_' + cond][this_trial] = (
                        model_annual_sum)

    # Output plots
    if output_plot:
        for year in final_rslt_dict.keys():
            try:
                plot_data(final_rslt_dict[year])
            except:
                continue

    return final_rslt_dict