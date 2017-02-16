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
#    warnings.showwarning = dt_fm.custom_warning
    for var in ['Fsd', 'ustar']:
        
        flag_index = np.where(~np.isnan(data_dict['NEE_series']) & 
                              np.isnan(data_dict[var]))
        count_str = str(len(flag_index[0]))
        if not len(flag_index[0]) == 0:
            warnings.warn('There are {0} instances where NEE element contains '
                          'valid data and {1} element is not a number - NEE '
                          'values for these instances will be excluded from '
                          'analysis!'.format(count_str, var))
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
            warnings.warn('There are {0} instances where neither the NEE nor '
                          'model driver {1} element contains valid data - '
                          'model estimates cannot be calculated for these '
                          'instances!'.format(count_str, var))
        data_dict['NEE_series'][flag_index] = np.nan        
    arr = np.unique(arr)
    if not len(arr) == 0:
        print 'Total number of instances with missing driver data is ' + str(len(arr))
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def check_path(path, generate = False):
    if not os.path.exists(path):
        if generate:
            try: 
                os.makedirs(path)
            except:
                warnings.warn('The specified path is incompatible with your OS!')
                new_path = os.path.join(os.path.expanduser('~'), 'Uncertainty')
                if not os.path.exists(new_path):
                    os.makedirs(new_path)
        else:
            raise IOError('Path {0} is invalid'.format(path))

    return
    
#------------------------------------------------------------------------------    
        
#------------------------------------------------------------------------------
# Set all sigma delta values to nan where there are no observational data
def filter_sigma_delta(data_dict):
    data_dict['sigma_delta'][np.isnan(data_dict['NEE_series'])] = np.nan
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------    
# Generate a random estimate of ustar threshold constrained by the mean and 
# confidsence interval for u* threshold derived from CPD
def generate_random_ustar_threshold(ustar_threshold, ustar_uncertainty):
    
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
    return this_ustar
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Fetch data from configurations
def get_data(configs_dict):
    
    data_file = os.path.join(configs_dict['files']['input_file'])
    var_list = configs_dict['variables'].values()

    file_ext = os.path.splitext(data_file)

    if file_ext[1] == '.nc':    
        data_dict, attr = io.OzFluxQCnc_to_data_structure(data_file, 
                                                          var_list = var_list, 
                                                          QC_var_list = ['Fc'], 
                                                          return_global_attr = True)
    elif file_ext[1] == '.df':
        data_dict, attr = io.DINGO_df_to_data_structure(data_file,
                                                        var_list = var_list,
                                                        return_global_attr = True)

    configs_dict['global_options']['measurement_interval'] = int(attr['time_step'])
        
    names_dict = dt_fm.get_standard_names(convert_dict = configs_dict['variables'])
    data_dict = dt_fm.rename_data_dict_vars(data_dict, names_dict)

    if configs_dict['global_options']['use_storage']:
        data_dict['NEE_series'] = data_dict['NEE_series'] + data_dict['Sc']
    elif configs_dict['global_options']['unify_flux_storage_cases']:
        data_dict['NEE_series'][np.isnan(data_dict['Sc'])] = np.nan 

    data_dict['PAR'] = data_dict['Fsd'] * 0.46 * 4.6       
        
    return data_dict    
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def get_ustar_data(path, years_list):
    
    try:
        check_path(path)
    except:
        raise IOError('Invalid path to ustar_threshold file')
    struct_array = np.genfromtxt(path, delimiter = ',', names = True)
    year_array = struct_array['Year']
    ustar_mean_array = struct_array['ustar_mean']
    ustar_sig_array = struct_array['ustar_sig']
    
    ustar_dict = {'ustar_threshold': {},
                  'ustar_uncertainty': {}}
                  
    crossref_years_list = []              
    for i in xrange(len(year_array)):
        this_year = int(year_array[i])
        crossref_years_list.append(this_year)
        if this_year in years_list:
            ustar_dict['ustar_threshold'][str(this_year)] = ustar_mean_array[i]
            ustar_dict['ustar_uncertainty'][str(this_year)] = ustar_sig_array[i]
            
    missing_years_list = [item for item in years_list if not item in crossref_years_list]
    
    if missing_years_list == years_list:
        raise Exception('None of the years in the data file were found in the '
                        'ustar_threshold file - defaulting to configuration file '
                        'values!')
    elif len(missing_years_list) != 0:
        missing_years_string = ', '.join([str(item) for 
                                          item in missing_years_list])
        print ('Year(s) {0} present in data file but not found in the '
               'ustar_threshold file - calculating mean of all available years '
               'to fill missing year(s)!'.format(missing_years_string))
        mean_ustar = ustar_mean_array.mean()
        mean_sig = ustar_sig_array.mean()
        for this_year in missing_years_list:
            ustar_dict['ustar_threshold'][str(this_year)] = mean_ustar
            ustar_dict['ustar_uncertainty'][str(this_year)] = mean_sig
            
    return ustar_dict
    
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Make a results dictionary
def init_interm_rslt_dict(num_trials, do_ustar, do_random, do_model):

    var_list = ['obs_avail_day', 'obs_avail_night']
    if do_ustar:
        var_list = var_list + ['ustar', 'ustar_error']
    if do_random:
        var_list = var_list + ['random_error_day', 'random_error_night']
    if do_model:
        var_list = var_list + ['model_error_day', 'model_error_night']

    nan_array = np.zeros(num_trials)
    nan_array[:] = np.nan
    return {var: cp.copy(nan_array) for var in var_list} 
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
    
    if len(error_list) == 0:
        raise KeyError('Passed dictionary contained no error data')    
    elif len(error_list) > 1: 
        results_d['total'] = np.zeros(len(results_d[error_list[0]]))
        for var in error_list:
            results_d['total'] = results_d['total'] + results_d[var]        
        bool_array = ~np.isnan(results_d['total'])
        error_list.append('total')
        for var in error_list:
            results_d[var] = results_d[var][bool_array]
        for var in error_list:
            if not var in ['total', 'ustar']:
                results_d[var] = results_d[var] + results_d['total'].mean()                        
                
    colors_d = {'total': 'grey',
                'ustar': 'blue',
                'random': 'cyan',
                'model': 'magenta'}

    # Set positions
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
    if 'total' in error_list:
        ax1.axvline(mu_dict['total'], color = 'black', 
                    linewidth = 2, linestyle = '--')
    else:
        ax1.axvline(mu_dict[error_list[0]], color = 'black', 
                    linewidth = 2, linestyle = '--')
    ax1.xaxis.set_ticks_position('bottom')
    ax1.yaxis.set_ticks_position('left')
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    
    # Plot the histogram
    for var in error_list:
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
    if 'total' in error_list:
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
            
    return fig
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
def main(input_config_file = False,
         input_data_file = False, 
         input_ustar_file = False, 
         num_trials = False, 
         output_directory = False, 
         output_plot = True):  
    
    # Update
    reload(logging)
    reload(rand_err)
    reload(mod_err)
    reload(io)
    reload(filt)
    reload(re)
    reload(gf)
    reload(dt_fm)
    
    # Set plotting to off in case running from iPython with interactive plotting
    if plt.isinteractive():
        is_on = True
        plt.ioff()
    else:
        is_on = False

                
    #-----------------------------------
    # General preparation and formatting
    #-----------------------------------

    # Get master config file (if not supplied)
    if input_config_file:
        configs_master_dict = io.config_to_dict(input_config_file)
    else:
        configs_master_dict = io.config_to_dict(io.file_select_dialog())

    # Build custom configuration file for this script
    configs_dict = build_config_file(configs_master_dict)

    # Build required configuration files for imported scripts (random error,
    # model error, respiration, light response)
    rand_err_configs_dict = configs_master_dict['random_error_configs']['options']
    mod_err_configs_dict = configs_master_dict['model_error_configs']['options']
    re_configs_dict = configs_master_dict['respiration_configs']['options']
    ps_configs_dict = configs_master_dict['photosynthesis_configs']['options']

          
    #--------------------------------------------------------------------
    # Write user overrides to configuration file OR to local namespace as 
    # appropriate, and force internal override to switch off plotting in
    # external functions
    #--------------------------------------------------------------------
    
    # User overrides
    if input_data_file:
        configs_dict['files']['input_file'] = input_data_file
    if not output_directory:
        output_directory = configs_dict['files']['output_path']  
    check_path(output_directory, generate = True)
    if not num_trials:
        num_trials = configs_dict['uncertainty_options']['num_trials']

    # Internal override (turn off respiration and light response plot functions)
    if re_configs_dict['output_fit_plots']:
        re_configs_dict['output_fit_plots'] = False 
    if ps_configs_dict['output_fit_plots']:
        ps_configs_dict['output_fit_plots'] = False 

        
    #---------------------------------
    # Logging setup and initialisation
    #---------------------------------
    
    log_file = 'log.txt'
    full_fname = os.path.join(output_directory, log_file)
    logging.basicConfig(filename = full_fname, level = logging.DEBUG)
    time_str = dt.datetime.strftime(dt.datetime.now(), '%Y-%m-%d %H:%M:%S')
    logging.info('\nRunning uncertainty analysis: {}\n'.format(time_str))
    warnings.showwarning = dt_fm.send_warnings_to_log        

        
    #-------------------------------    
    # Data retrieval and preparation
    #-------------------------------

    # Get data    
    data_dict = get_data(configs_dict)

    # Save the time step information into the individual configuration files
    measurement_interval = (configs_dict['global_options']
                                        ['measurement_interval'])
    for d in [rand_err_configs_dict, mod_err_configs_dict, 
              re_configs_dict, ps_configs_dict]: 
        d['measurement_interval'] = measurement_interval    
    
    # Check no NEE values with missing ustar or Fsd values
    check_data_consistency(data_dict)

    # Check no drivers missing where NEE is missing
    check_driver_consistency(data_dict)

    # Get ustar data
    if input_ustar_file:
        years_list = list(set([this_date.year for 
                               this_date in data_dict['date_time']]))
        try:
            ustar_dict = get_ustar_data(input_ustar_file, years_list)
            ustar_threshold = ustar_dict['ustar_threshold']
            ustar_uncertainty = ustar_dict['ustar_uncertainty']
        except IOError, e:
            print('Retrieval of ustar_threshold file failed with the '
                  'following message: {0}. Using defaults provided '
                  'in configuration file...'.format(e))
    else:
        ustar_threshold = configs_dict['global_options']['ustar_threshold']    
        if do_ustar_uncertainty: ustar_uncertainty = (
            configs_dict['global_options']['ustar_uncertainty'])


    #----------------------------------------------
    # Assign configuration items to local namespace
    #----------------------------------------------
    
    noct_threshold = configs_dict['global_options']['noct_threshold']
    ustar_filter_day = configs_dict['global_options']['ustar_filter_day']
    do_ustar_uncertainty = (configs_dict['uncertainty_options']
                                        ['do_ustar_uncertainty'])
    do_random_uncertainty = (configs_dict['uncertainty_options']
                                         ['do_random_uncertainty'])
    do_model_uncertainty = (configs_dict['uncertainty_options']
                                        ['do_model_uncertainty'])
    NEE_model = configs_dict['uncertainty_options']['NEE_model']
            
    
    #---------------------------------------------
    # Print info and check for configuration error
    #---------------------------------------------
    
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
        raise Exception('Processing flags for all uncertainty sources '
                        'currently set to False: set at least one '
                        'uncertainty source to True in configuration file '
                        'before proceeding!')
    print '---------------------------------'
    

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

    # Print stuff
    print '----------------------------'
    print 'Starting uncertainty trials:',
    
    # If there are > 100 trials, create arrays for percentiles
    if num_trials > 100:
        percentile_array = np.linspace(0, 95, 20)
        num_array = (np.round(percentile_array / 100 
                              * num_trials)).astype(int)
    
    # Do trials
    for this_trial in xrange(num_trials):

        logging.info('\nRunning trial #: {0}\n'.format(str(this_trial)))
        
        # Set first_pass flag to prevent repetitive assignment
        first_pass = True if this_trial == 0 else False

        # Print progress
        if not this_trial == num_trials - 1:
            if num_trials <= 100:
                print str(this_trial + 1) + ',',
            else:
                if this_trial in num_array:
                    this_pctl = int(percentile_array[np.where(num_array == 
                                                              this_trial)])
                    print str(this_pctl) + '%,',
        else:
            print str(this_trial + 1) + ' ... Done!'

        # If doing ustar uncertainty, make a ustar dictionary by randomly 
        # sampling from the normal distribution and scaling according to
        # ustar threshold uncertainty
        if do_ustar_uncertainty:
            this_ustar = generate_random_ustar_threshold(ustar_threshold,
                                                         ustar_uncertainty)
              
            # Make a copy of the data dictionary
            this_dict = cp.deepcopy(data_dict)
            
            # Screen low ustar then model and gap fill
            filt.screen_low_ustar(this_dict, this_ustar, noct_threshold, 
                                  ustar_filter_day)
            try:
                run_model(this_dict, NEE_model, re_configs_dict, ps_configs_dict)
                fail_flag = False
            except Exception, e:
                warnings.warn('Model optimisation for trial {0} failed with '
                              'the following message: {1}\n'.format(str(this_trial), 
                                                                    e[0]))
                fail_flag = True
    
            # If doing random uncertainty, screen out any estimates of 
            # sigma_delta where there are no obs
            if do_random_uncertainty: filter_sigma_delta(this_dict)
            
        # If not doing ustar uncertainty, just assign this_dict to data_dict
        elif first_pass:
            fail_flag = False    
            this_dict = data_dict
            filt.screen_low_ustar(this_dict, ustar_threshold, noct_threshold)
            run_model(this_dict, NEE_model, re_configs_dict, ps_configs_dict)
            if do_random_uncertainty: filter_sigma_delta(this_dict)
                
        # Create dataset separated into years
        years_data_dict = dtf.subset_datayear_from_arraydict(this_dict, 
                                                             'date_time')
        
        # Do calculations for each year
        for this_year in years_list:
            
            if do_ustar_uncertainty:
                final_rslt_dict[this_year]['ustar'][this_trial] = (
                    this_ustar[str(this_year)])
                if not fail_flag:
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
                    
                if not fail_flag:
                            
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
        plot_dict = {}
        for year in final_rslt_dict.keys():
            try:
                plot_dict[year] = plot_data(final_rslt_dict[year])
            except:
                continue

    if is_on:
        plt.ion()                    
    
    if output_plot:
        return final_rslt_dict, plot_dict
    else:
        return final_rslt_dict