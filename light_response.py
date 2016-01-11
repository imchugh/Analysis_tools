# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 11:23:04 2015

@author: imchugh
"""
# Python modules
import os
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import pdb

# My modules
import datetime_functions as dtf
import data_filtering as filt
import gap_filling as gf
import light_and_T_response_functions as light

def calculate_light_response(data_dict,
                             configs_dict,
                             params_in_dict,
                             params_out_dict):

    # Create local vars from configs
    meas_int = configs_dict['measurement_interval']    
    min_pct = configs_dict['minimum_pct_window']    

    # Calculate rb for windows    
    for date in data_dict.keys():

        # Make a temporary dict with nans and daytime dropped
        bool_filter = data_dict[date]['all_bool']
        temp_dict = {var: data_dict[date][var][bool_filter] 
                     for var in ['TempC', 'VPD', 'PAR', 'NEE_series']}

        # Get the index for the current date
        date_index = np.where(params_out_dict['date'] == date)

        # Check data availability 
        data_pct = int(len(temp_dict['NEE_series']) / 
                       float((1440 / meas_int) / 2) * 100)

        # Do (or do not - there is no try... actually there is!) the fit
        if not data_pct < min_pct:

            # Get values of respiration parameters
            params_in_dict['Eo_default'] = params_out_dict['Eo'][date_index]
            
            # If using nocturnal rb...
            if configs_dict['use_nocturnal_rb']:
                
                params_in_dict['rb_default'] = params_out_dict['rb'][date_index]
                fit_dict = light.optimise_fixed_rb(temp_dict, params_in_dict)
            
            # If using daytime rb...                                                             
            else:
                fit_dict = light.optimise_free_rb(temp_dict, params_in_dict)
            
            # Write data to results arrays
            for key in fit_dict.keys():
                params_out_dict[key][date_index] = fit_dict[key]

            # Write alpha default values to default dictionary if valid
            if fit_dict < 2:
                params_in_dict['alpha_default'] = fit_dict['alpha']
            else:
                params_in_dict['alpha_default'] = 0

        else:

            # Error code for not enough data
            params_out_dict['error_code'][date_index] = 10

    # Interpolate
    for key in fit_dict.keys():
        if not key == 'error_code':
            params_out_dict[key] = gf.generic_2d_linear(params_out_dict[key])

    # Rename the error code variable
    params_out_dict['light_response_error_code'] = params_out_dict.pop('error_code')
    
    return

def estimate_GPP_Re(data_dict, 
                    all_params_dict, 
                    datetime_input_index_dict):
    
    # Create output arrays for estimates
    results_dict = {}
    for var in ['Re', 'GPP']:
        results_dict[var] = np.empty(len(data_dict['TempC']))
        results_dict[var][:] = np.nan
    
    # Estimate time series GPP and Re and write to results dictionary
    for i, date in enumerate(all_params_dict['date']):
        
        indices = datetime_input_index_dict[date]
        this_Eo = all_params_dict['Eo'][i]
        this_rb = all_params_dict['rb'][i]
        this_alpha = all_params_dict['alpha'][i]
        this_beta = all_params_dict['beta'][i]
        this_k = all_params_dict['k'][i]
        
        this_dict = {var: data_dict[var][indices[0]: indices[1] + 1]
                     for var in ['NEE_series', 'VPD', 'PAR', 'TempC']}
        GPP, Re = light.LRF_part(this_dict, 
                                 this_Eo,
                                 this_rb, 
                                 this_alpha,
                                 this_beta,
                                 this_k) 
        results_dict['GPP'][indices[0]: indices[1] + 1] = GPP
        results_dict['Re'][indices[0]: indices[1] + 1] = Re

    return results_dict

#def append_results_array(dates_input_index_dict):
#    
#    date_array = np.array(dates_input_index_dict.keys())
#    date_array.sort()
#    generic_array = np.empty(len(dates_input_index_dict))
#    generic_array[:] = np.nan
#    rb_error_code_array = np.ones(len(generic_array)) * 20
#    return {'date': date_array,
#            'Eo': cp.copy(generic_array),
#            'Eo_error_code': cp.copy(generic_array),
#            'rb': cp.copy(generic_array),
#            'rb_error_code': rb_error_code_array}  

# Nocturnal fits for each window
def plot_windows(step_data_dict, configs_dict, params_dict):

    # Don't send plots to screen
    if plt.isinteractive():
        is_on = True
        plt.ioff()
    else:
        is_on = False

    # Set parameters from dicts
    window = configs_dict['window_size_days']
    
    x_lab = '$PPFD\/(\mu mol\/photons\/m^{-2} s^{-1})$'
    
    for date in step_data_dict.keys():

        Eo = params_dict['Eo'][params_dict['date'] == date]
        rb = params_dict['rb'][params_dict['date'] == date]
        alpha = params_dict['alpha'][params_dict['date'] == date]
        beta = params_dict['beta'][params_dict['date'] == date]
        k = params_dict['k'][params_dict['date'] == date]

        bool_filter = step_data_dict[date]['day_bool']
        x_var = step_data_dict[date]['PAR'][bool_filter]
        y_var1 = step_data_dict[date]['NEE_series'][bool_filter]
        index = x_var.argsort()
        x_var = x_var[index]
        y_var1 = y_var1[index]
        GPP, Re = light.LRF_part(step_data_dict[date], 
                                 Eo, 
                                 rb,
                                 alpha,
                                 beta,
                                 k)
        y_var2 = (GPP + Re)[bool_filter]
        y_var2 = y_var2[index]

        # Plot
        date_str = dt.datetime.strftime(date,'%Y-%m-%d')
        fig = plt.figure(figsize = (12,8))
        fig.patch.set_facecolor('white')
        ax = plt.gca()
        ax.plot(x_var, y_var1, 'o' , markerfacecolor = 'none',
                 markeredgecolor = 'black', label = 'Observed', color = 'black')
        ax.plot(x_var, y_var2, '^', color = 'black', label = 'Estimated')
        ax.set_title('Fit for ' + str(window) + ' day window centred on ' + 
                      date_str + '\n', fontsize = 22)
        ax.set_xlabel(x_lab, fontsize = 18)
        ax.set_ylabel('$NEE\/(\mu mol C\/m^{-2} s^{-1}$)', fontsize = 18)
        ax.axhline(y = 0, color = 'black')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.legend(loc = 'upper right', fontsize = 16, frameon = False, 
                  numpoints = 1)
        plot_out_name = 'day' + '_' + date_str + '.jpg'
        plt.tight_layout()
        fig.savefig(os.path.join(configs_dict['output_path'],
                                 plot_out_name))
        plt.close(fig)

    if is_on:
        plt.ion()
        
    return
                
#------------------------------------------------------------------------------
#############
# Main code #
#############        

def main(data_dict, configs_dict, params_out_dict):

    """
    Calculates Re using Lloyd and Taylor function where Eo is fitted to annual
    data and rb is fitted to specified combination of step and window;
    Pass: 1) a data dictionary containing the following key / value pairs (it 
             is assumed that there are no gaps in the time series, but this is 
             not yet enforced!; also, all values in dict must be numpy arrays, 
             and all must be of same length):
                 - 'date_time': numpy array of Python datetimes
                 - 'NEE_series': numpy array of the NEE time series to be used 
                                as the optimisation target (float); note:
                                - missing values must be np.nan
                                - no QC is done on values - this must be done 
                                  prior to passing the data
                 - 'TempC': numpy array of temperatures to be used as 
                            optimisation input (float)
          2) a configs dict containing the following key / value pairs:
                 - 'step_size_days': step size in days between fitting windows
                                     (int; range 0 < x < n days)
                 - 'window_size_days': width of fitting window in days 
                                     (int; range 0 < x < n days)
                 - 'minimum_pct_annual': minimum acceptable percentage of 
                                         available annual data for fitting of 
                                         Eo (int; range 0 <= x <= 100)
                 - 'minimum_pct_window': minimum acceptable percentage of 
                                         available window data for fitting of 
                                         rb (int; range 0 <= x <= 100)
                 - 'measurement_interval': measurement interval (minutes) of 
                                           the input data (int)                                          
    Returns: 1) a results dictionary containing 2 key / value pairs:
                    - 'date_time': numpy array of Python datetimes for each
                      datum in the original time series
                    - 'Re': numpy array of half-hourly estimates of Re
             2) a results dictionary containing 5 key / value pairs:
                    - 'date_time': numpy array of Python dates for each day in 
                                   the original time series 
                    - 'Eo': numpy array of Eo estimates for each day (note that 
                            each year has a constant value)
                    - 'Eo_error_code': numpy array of Eo diagnostic errors for 
                                       each day (note that each year has a 
                                       constant value)
                    - 'rb': numpy array of rb estimates for each day (note that 
                            linear interpolation is used to gap fill between 
                            steps)
                    - 'rb_error_code': numpy array of rb diagnostic errors for 
                                       each day (including whether the value is 
                                       interpolated or calculated)                    
    """
    
    #------------------------------------------------------------------------------

    # Create boolean indices for masking daytime and nan values
    day_mask = data_dict['Fsd'] > 5
    nan_mask = filt.subset_arraydict_on_nan(data_dict,
                                            var_list = ['NEE_series', 'TempC',
                                                        'VPD', 'PAR'],
                                            subset = False)
    all_mask = [all(rec) for rec in zip(day_mask, nan_mask)]
    data_dict['day_bool'] = np.array(day_mask)
    data_dict['all_bool'] = np.array(all_mask)

    # Partition the data into year and step pieces
    step_data_dict = dtf.get_moving_window(data_dict, 
                                           'date_time', 
                                           configs_dict['window_size_days'],
                                           configs_dict['step_size_days'])

    # Get the indices of the start and end rows of each unique date in the source 
    # data array - no data dict is built from this, since these indices are used to
    # assign Re estimates to the estimated time series output array only
    dates_input_index_dict = dtf.get_day_indices(data_dict['date_time'])
    
    # Generate a results dictionary for the parameter values (1 for each day)
    if configs_dict['use_nocturnal_rb']:
        rslt_var_list = ['alpha', 'beta', 'k', 'error_code']
    else:
        rslt_var_list = ['rb_day', 'alpha', 'beta', 'k', 'error_code']
    for var in rslt_var_list:
        params_out_dict[var] = np.empty(len(params_out_dict['rb']))
        if not var == 'error_code':
            params_out_dict[var][:] = np.nan
        else:
            params_out_dict[var][:] = 20

    # Initalise parameter dicts with prior estimates    
    params_in_dict = {'k_prior': 0,
                      'alpha_prior': -0.01,
                      'rb_prior': params_out_dict['rb'].mean(),
                      'beta_prior': (np.percentile(data_dict['NEE_series']
                                                   [data_dict['all_bool']], 5) - 
                                     np.percentile(data_dict['NEE_series']
                                                   [data_dict['all_bool']], 95)),
                      'alpha_default': 0,
                      'beta_default': 0,
                      'k_default': 0 }

    # Write the light response parameters to the existing parameters dict
    calculate_light_response(step_data_dict,
                             configs_dict,
                             params_in_dict,
                             params_out_dict)
    
    # Estimate Re and GPP for all data
    rslt_dict = estimate_GPP_Re(data_dict,
                                params_out_dict,
                                dates_input_index_dict)

    # Do plotting if specified
    if configs_dict['output_fit_plots']:
        plot_windows(step_data_dict, configs_dict, params_out_dict)
#
#    return rslt_dict, params_out_dict