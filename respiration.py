# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 13:47:20 2015

@author: imchugh
"""

import sys
import numpy as np
import pdb
import calendar

sys.path.append('../Partitioning')
import DataIO as io
import gap_filling as gf
import dark_T_response_functions as dark
import datetime_functions as dtf
import data_filtering as filt

def get_LT_fit_params(data_dict, configs_dict):
    """
    DOCSTRING!!!!!!!!!!!!!!!!!!
    """
    def filtering(this_dict):
        noct_dict = filt.subset_arraydict_on_threshold(this_dict, 'Fsd', 5, '<', 
                                                       drop = True)    
        ustar_dict = filt.subset_arraydict_on_threshold(noct_dict, 'ustar', 0.42, 
                                                        '>', drop = True)
        sub_dict = filt.subset_arraydict_on_nan(ustar_dict)
        return sub_dict

    window = configs_dict['window_size_days']
    step = configs_dict['step_size_days']
    
    # Create data step indices for input arrays
    datetime_array = data_dict.pop('date_time')
    years_input_index_dict = dtf.get_year_indices(datetime_array)
    step_dates_input_index_dict = dtf.get_moving_window_indices(datetime_array, 
                                                                window, step)
    
    # Create date step indices for output arrays
    dates_output_index_dict = dtf.get_unique_dates(datetime_array)
    date_array = np.array(dates_output_index_dict.keys())
    date_array.sort()
    years_output_index_dict = dtf.get_year_indices(date_array, retro_stamp = False)

#------------------------------------------------------------------------------
    
    # Create output arrays for parameter time series
    Eo_array = np.empty([len(dates_output_index_dict)])
    Eo_array[:] = np.nan
    
    # Initalise dicts
    all_noct_dict = filtering(data_dict)
    params_dict = {'Eo_prior': 100,
                   'rb_prior': all_noct_dict['NEE'].mean()}
    
    # Do annual fits for Eo
    Eo_annual_data_dict = {}
    Eo_annual_error_dict = {}
    Eo_pass_keys = []
    Eo_fail_keys = []
    for year in years_input_index_dict.keys():
    
        # Calculate number of nocturnal recs for year
        days = 366 if calendar.isleap(year) else 365
        recs = days * (24 / configs_dict['measurement_interval']) / 2
    
        # Input and output indices
        in_indices = years_input_index_dict[year]
        
        this_dict = {var: data_dict[var][in_indices[0]: in_indices[1]] 
                     for var in data_dict.keys()}
        sub_dict = filtering(this_dict)
        data_pct = int(len(sub_dict['NEE']) / float(recs) * 100)
        if not data_pct < configs_dict['minimum_pct_annual']:
            params, error_code = dark.optimise_all(sub_dict, params_dict)
        else:
            params, error_code = [np.nan, np.nan], 10
        Eo_annual_data_dict[year] = params[0]
        Eo_annual_error_dict[year] = error_code
        if error_code == 0: 
            Eo_pass_keys.append(year)
        else:
            Eo_fail_keys.append(year)
            
    # Fill gaps and project Eo to the appropriate indices of the results array
    if np.all(np.isnan(Eo_annual_data_dict.values())):
        print 'Could not find any values of Eo for any years! Exiting...'
        sys.exit()
    if np.any(np.isnan(Eo_annual_data_dict.values())):
        Eo_mean = np.array([Eo_annual_data_dict[year] for year in Eo_pass_keys]).mean()    
        for year in Eo_fail_keys:
            Eo_annual_data_dict[year] = Eo_mean
            Eo_pass_keys.append(year)
    for year in Eo_pass_keys:
        out_indices = years_output_index_dict[year]
        Eo_array[out_indices[0]: out_indices[1] + 1] = Eo_annual_data_dict[year]

    rb_array = np.empty([len(dates_output_index_dict)])
    rb_array[:] = np.nan
        
    # Calculate rb for windows    
    for date in step_dates_input_index_dict.keys():
        
        # Specify Eo value for the relevant year
        params_dict['Eo_default'] = Eo_annual_data_dict[date.year]
        
        # Input and output indices
        in_indices = step_dates_input_index_dict[date]
        out_index = dates_output_index_dict[date]
    
        # Do optimisation
        this_dict = {var: data_dict[var][in_indices[0]: in_indices[1]] 
                     for var in data_dict.keys()}
        sub_dict = filtering(this_dict)
        data_pct = int(len(sub_dict['NEE']) / float((24 / configs_dict['measurement_interval']) / 2) * 100)
        if not data_pct < configs_dict['minimum_pct_noct_window']:
            params, error_code = dark.optimise_rb(sub_dict, params_dict)
        else:
            params, error_code = [np.nan], 10
        rb_array[out_index] = params
    
    # Interpolate rb
    rb_array = gf.generic_2d_linear(rb_array)
    
    # Return a dictionary
    return {'date': date_array, 'Eo': Eo_array, 'rb': rb_array}

def estimate_Re(data_dict, params_dict):

    # Create date step indices for input arrays
    datetime_array = data_dict.pop('date_time')
    datetime_input_index_dict = dtf.get_day_indices(datetime_array)
    
    # Create output arrays for Re estimates
    results_array = np.empty(len(data_dict['TempC']))
    results_array[:] = np.nan
    
    # Estimate time series Re
    for i, date in enumerate(params_dict['date']):
        
        indices = datetime_input_index_dict[date]
        this_Eo = params_dict['Eo'][i]
        this_rb = params_dict['rb'][i]
        this_dict = {'TempC': data_dict['TempC'][indices[0]: indices[1] + 1]}
        results_array[indices[0]: indices[1] + 1] = dark.TRF(this_dict, 
                                                             this_Eo,
                                                             this_rb)
    return {'Re': results_array}