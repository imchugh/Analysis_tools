# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 11:23:04 2015

@author: imchugh
"""
# Python modules
import os
import numpy as np
import copy as cp
import calendar
import matplotlib.pyplot as plt
import datetime as dt
from scipy.optimize import curve_fit
import pdb

# My modules
import datetime_functions as dtf
import data_filtering as filt
import gap_filling as gf


#------------------------------------------------------------------------------
# Data optimisation algorithm
def TRF(data_dict, rb, Eo, theta1, theta2):
    
    T_func = rb  * np.exp(data_dict['TempC'] * Eo)
    sws_func = 1 / (1 + np.exp(theta1 - theta2 * data_dict['Sws']))
    
    return sws_func * T_func
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Write error messages to dictionary with codes as keys
def error_codes():
    
    d = {0:'Optimisation successful',
         1:'Value of Eo failed range check - rejecting all parameters',
         2:'Value of rb has wrong sign - rejecting all parameters',
         3:'Optimisation reached maximum number of iterations ' \
           'without convergence',
         10:'Data did not pass minimum percentage threshold - ' \
            'skipping optimisation',
         20:'Data are linearly interpolated from nearest non-interpolated ' \
            'neighbour'}
    
    return d
#------------------------------------------------------------------------------    

#------------------------------------------------------------------------------
# run optimisation and raise error code for combined rb and Eo
def optimise_all(data_dict, params_dict):

    # Initialise error state variable
    error_state = 0              
    drivers_dict = {driver: data_dict[driver] for driver in ['TempC', 'Sws']}
    response_array = data_dict['NEE_series']

    try:
        params = curve_fit(lambda x, a, b, c, d:
                           TRF(x, a, b, c, d),
                           drivers_dict, 
                           response_array, 
                           p0 = [params_dict['rb_prior'], 
                                 params_dict['Eo_prior'],
                                 params_dict['theta1_prior'],
                                 params_dict['theta2_prior']])[0]
    except RuntimeError:
        params = [np.nan, np.nan, np.nan, np.nan]
        error_state = 3

    # If negative rb returned, set to nan
    if params[0] < 0:
        error_state = 2
        params = [np.nan, np.nan, np.nan, np.nan]
    elif params[1] < -50 or params[1] > 400: 
        error_state = 1
        params = [np.nan, np.nan, np.nan, np.nan]


    return {'rb': params[0], 'Eo': params[1], 'theta1': params[2], 
            'theta2': params[3], 'error_code': error_state}
#------------------------------------------------------------------------------
    
#------------------------------------------------------------------------------    
# run optimisation and raise error code for rb with fixed Eo
def optimise_rb(data_dict, params_dict):

    # Initialise error state variable
    error_state = 0              
    
    # Get drivers and response
    drivers_dict = {driver: data_dict[driver] for driver in ['TempC', 'Sws']}
    response_array = data_dict['NEE_series']        
    
    try:
        params = curve_fit(lambda x, b:
                           TRF(x, b, params_dict['Eo_default'],
                               params_dict['theta1_default'],
                               params_dict['theta2_default']),
                           drivers_dict, 
                           response_array, 
                           p0 = [params_dict['rb_prior']])[0]                           
    except RuntimeError:
        params = [np.nan]

    # If negative rb returned, set to nan
    if params[0] < 0:
        error_state = 2
        params = [np.nan]
       
    return {'rb': params[0], 'error_code': error_state}
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Iterate through dates and write optimisation results for rb parameter to
# parameter results dictionary
def calculate_rb(data_dict,
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
                     for var in ['TempC', 'NEE_series', 'Sws']}
        
        # Get the index for the current date
        date_index = np.where(params_out_dict['date'] == date)
        
        # Check data availability 
        data_pct = int(len(temp_dict['NEE_series']) / 
                       float((1440 / meas_int) / 2) * 100)
        
        # If enough data, go ahead
        if not data_pct < min_pct:

            # Do fit
            for var in ['Eo', 'theta1', 'theta2']:
                params_in_dict[var + '_default'] = params_in_dict[var + '_calc'][date.year]
#            params_in_dict['Eo_default'] = params_in_dict['Eo_calc'][date.year]
            fit_dict = optimise_rb(temp_dict, params_in_dict)
            fit_dict['rb_error_code'] = fit_dict.pop('error_code')

            # Write data to results arrays
            for key in fit_dict.keys():
                params_out_dict[key][date_index] = fit_dict[key]
            
        else:

            # Error code for not enough data
            params_out_dict['rb_error_code'][date_index] = 10

    # Interpolate rb
    params_out_dict['rb'] = gf.generic_2d_linear(params_out_dict['rb'])

    return
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Iterate through years and write optimisation results for Eo parameter to
# parameter results dictionary
def calculate_Eo(data_dict,
                 configs_dict,
                 params_in_dict,
                 params_out_dict):

    # Create local vars from configs
    meas_int = configs_dict['measurement_interval']
    min_pct = configs_dict['minimum_pct_annual']

    # Do annual fits for Eo
    years_list = data_dict.keys()
    params_temp_dict = {}
    fail_years_list = []
    vars_list = ['rb', 'Eo', 'theta1', 'theta2']
    for this_year in years_list:

        # Dict for this year only (overwritten on subsequent pass)
        fit_dict = {}
        
        # Make a temporary dict with nans and daytime dropped
        bool_filter = data_dict[this_year]['all_bool']
        temp_dict = {var: data_dict[this_year][var][bool_filter] 
                     for var in ['TempC', 'NEE_series', 'Sws']}

        # Calculate number of nocturnal recs for year
        days = 366 if calendar.isleap(this_year) else 365
        recs = days * (1440 / meas_int) / 2

        # Input indices
        data_pct = int(len(temp_dict['NEE_series']) / float(recs) * 100)
        if not data_pct < min_pct:
            fit_dict = optimise_all(temp_dict, params_in_dict)
        else:
            fit_dict = {'Eo': np.nan, 'theta1': np.nan, 
                        'theta2': np.nan, 'error_code' : 10}
        if not fit_dict['error_code'] == 0: 
            fail_years_list.append(this_year)
            
        params_temp_dict[this_year] = fit_dict
            
    # Fill gaps 
    pass_years_list = list(set(years_list) - set(fail_years_list))
    if len(pass_years_list) == 0:
        raise Exception('Could not find any values of Eo for selected year(s)!')
    if len(pass_years_list) < len(data_dict.keys()):
        for param in vars_list:
            param_sum = 0
            for this_year in pass_years_list:
                param_sum = param_sum + params_temp_dict[this_year][param]
            param_mean = param_sum / len(pass_years_list)
            for this_year in fail_years_list:
                params_temp_dict[this_year][param] = param_mean

    # Attach the yearly Eo to the parameter dictionary
    for var in vars_list:
        params_in_dict[var + '_calc'] = {this_year: params_temp_dict[this_year][var] 
                                         for this_year in years_list}

    # Project Eo to the appropriate indices of the results array
    vars_list.remove('rb')
    years_array = np.array([date.year for date in params_out_dict['date']])
    for this_year in years_list:
        index = np.where(years_array == this_year)[0]
        for var in vars_list:
            params_out_dict[var][index] = (params_temp_dict[this_year][var])
        params_out_dict['annual_error_code'][index] = (
            params_temp_dict[this_year]['error_code'])

    return
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Use observed meteorology and optimised respiration parameters to estimate Re
def estimate_Re(data_dict, 
                all_params_dict, 
                datetime_input_index_dict):
    
    # Create output dicts for estimates
    results_dict = {}
    results_dict['Re'] = np.empty(len(data_dict['TempC']))
    results_dict['Re'][:] = np.nan
    results_dict['date_time'] = data_dict['date_time']
    
    # Estimate time series Re
    for i, date in enumerate(all_params_dict['date']):
        
        indices = datetime_input_index_dict[date]
        this_Eo = all_params_dict['Eo'][i]
        this_rb = all_params_dict['rb'][i]
        this_theta1 = all_params_dict['theta1'][i]
        this_theta2 = all_params_dict['theta2'][i]
        this_dict = {'TempC': data_dict['TempC'][indices[0]: indices[1] + 1],
                     'Sws': data_dict['Sws'][indices[0]: indices[1] + 1]}
        Re = TRF(this_dict, this_rb, this_Eo, this_theta1, this_theta2)                                                         
        results_dict['Re'][indices[0]: indices[1] + 1] = Re

    return results_dict
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Build the dictionary for 
def generate_results_dict(datetime_array):
    
    dates_input_index_dict = dtf.get_day_indices(datetime_array)
    date_array = np.array(dates_input_index_dict.keys())
    date_array.sort()
    generic_array = np.empty(len(dates_input_index_dict))
    generic_array[:] = np.nan
    rb_error_code_array = np.ones(len(generic_array)) * 20
    return {'date': date_array,
            'Eo': cp.copy(generic_array),
            'annual_error_code': cp.copy(generic_array),
            'rb': cp.copy(generic_array),
            'rb_error_code': rb_error_code_array,
            'theta1': cp.copy(generic_array),
            'theta2': cp.copy(generic_array)}  
#------------------------------------------------------------------------------
            
#------------------------------------------------------------------------------            
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
    
    x_lab = '$Temperature\/(^{o}C$)'
    
    for date in step_data_dict.keys():

        Eo = params_dict['Eo'][params_dict['date'] == date]
        rb = params_dict['rb'][params_dict['date'] == date]

        bool_filter = step_data_dict[date]['night_bool']
        x_var = step_data_dict[date]['TempC'][bool_filter]
        y_var1 = step_data_dict[date]['NEE_series'][bool_filter]
        index = x_var.argsort()
        x_var = x_var[index]
        y_var1 = y_var1[index]
        y_var2 = TRF({'TempC': x_var}, Eo, rb)
          
        # Plot
        date_str = dt.datetime.strftime(date,'%Y-%m-%d')
        fig = plt.figure(figsize = (12,8))
        fig.patch.set_facecolor('white')
        ax = plt.gca()
        ax.plot(x_var, y_var1, 'o' , markerfacecolor = 'none',
                 markeredgecolor = 'black', label = 'NEE_obs', color = 'black')
        ax.plot(x_var, y_var2, linestyle = ':', color = 'black', 
                 label = 'NEE_est')
        ax.set_title('Fit for ' + str(window) + ' day window centred on ' + 
                      date_str + '\n', fontsize = 22)
        ax.set_xlabel(x_lab, fontsize = 18)
        ax.set_ylabel('$NEE\/(\mu mol C\/m^{-2} s^{-1}$)', fontsize = 18)
        ax.axhline(y = 0, color = 'black')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plot_out_name = 'noct' + '_' + date_str + '.jpg'
        plt.tight_layout()
        fig.savefig(os.path.join(configs_dict['output_path'],
                                 plot_out_name))
        plt.close(fig)
    
    if is_on:
        plt.ion()
        
    return
#------------------------------------------------------------------------------
                
#------------------------------------------------------------------------------
#############
# Main code #
#############        
def main(data_dict, configs_dict):

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
                 - 'output_fit_plots': whether to output plots showing the fit
                                       of the parameters to the data (boolean)
                 - 'output_path': path for output of plots (str)
                 
    Returns: 1) a results dictionary containing 2 key / value pairs:
                    - 'date_time': numpy array of Python datetimes for each
                      datum in the original time series
                    - 'Re': numpy array of half-hourly estimates of Re
                    
             2) a results dictionary containing 5 key / value pairs:
                    - 'date_time': numpy array of Python dates for each day in 
                                   the original time series 
                    - 'Eo': numpy array of Eo estimates for each day (note that 
                            each year has a constant value)
                    - 'annual_error_code': numpy array of Eo diagnostic errors for 
                                           each day (note that each year has a 
                                           constant value)
                    - 'rb': numpy array of rb estimates for each day (note that 
                            linear interpolation is used to gap fill between 
                            steps)
                    - 'rb_error_code': numpy array of rb diagnostic errors for 
                                       each day (including whether the value is 
                                       interpolated or calculated)                    
    """

    # Create boolean indices for masking daytime and nan values
    night_mask = data_dict['Fsd'] < 5
    nan_mask = filt.subset_arraydict_on_nan(data_dict,
                                            var_list = ['NEE_series', 'TempC'],
                                            subset = False)
    all_mask = [all(rec) for rec in zip(night_mask, nan_mask)]
    data_dict['night_bool'] = np.array(night_mask)
    data_dict['all_bool'] = np.array(all_mask)

    # Partition the data into year and step pieces
    years_data_dict = dtf.get_year_window(data_dict,
                                          'date_time')

    step_data_dict = dtf.get_moving_window(data_dict, 
                                           'date_time', 
                                           configs_dict['window_size_days'],
                                           configs_dict['step_size_days'])

    # Get the indices of the start and end rows of each unique date in the source 
    # data array - no data dict is built from this, since these indices are used to
    # assign Re estimates to the estimated time series output array only
    dates_input_index_dict = dtf.get_day_indices(data_dict['date_time'])

    # Initalise parameter dicts with prior estimates
    params_in_dict = {'Eo_prior': 0.01,
                      'rb_prior': data_dict['NEE_series'][data_dict['all_bool']]
                      .mean(),
                      'theta1_prior': 1,
                      'theta2_prior': 10}
    
    # Generate a results dictionary for the parameter values (1 for each day)
    params_out_dict = generate_results_dict(data_dict['date_time'])

    # Get Eo for all years
    calculate_Eo(years_data_dict, 
                 configs_dict,
                 params_in_dict,
                 params_out_dict)
    
    # Get rb for all steps
    calculate_rb(step_data_dict,
                 configs_dict,
                 params_in_dict,
                 params_out_dict)

    # Estimate Re for all data
    rslt_dict = estimate_Re(data_dict,
                            params_out_dict,
                            dates_input_index_dict)

    # Get error codes
    error_dict = error_codes()


    # Do plotting if specified
    if configs_dict['output_fit_plots']:
        plot_windows(step_data_dict, configs_dict, params_out_dict)

    return rslt_dict, params_out_dict, error_dict
#------------------------------------------------------------------------------
