# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 11:23:04 2015

@author: imchugh
"""
import sys
import numpy as np
import copy as cp
import calendar

#sys.path.append('../Analysis_tools')
#sys.path.append('../Partitioning')
import datetime_functions as dtf
import data_filtering as filt
import gap_filling as gf
import dark_T_response_functions as dark

def calculate_rb(data_dict,
                 configs_dict,
                 params_in_dict,
                 params_out_dict):

    # Create local vars from configs
    meas_int = configs_dict['measurement_interval']    
    min_pct = configs_dict['minimum_pct_noct_window']    
    
    # Calculate rb for windows    
    for date in data_dict.keys():
        
        # Specify Eo value for the relevant year
        params_in_dict['Eo_default'] = params_in_dict['Eo_years'][date.year]
        data_pct = int(len(data_dict[date]['Fc_series']) / 
                       float((1440 / meas_int) / 2) * 100)
        if not data_pct < min_pct:
            params, error_code = dark.optimise_rb(data_dict[date], 
                                                  params_in_dict)
        else:
            params, error_code = [np.nan], 10
        out_index = np.where(params_out_dict['date'] == date)
        params_out_dict['rb'][out_index] = params
        params_out_dict['rb_error_code'][out_index] = error_code

    # Interpolate rb
    params_out_dict['rb'] = gf.generic_2d_linear(params_out_dict['rb'])
    
    return

def calculate_Eo(data_dict,
                 configs_dict,
                 params_in_dict,
                 params_out_dict):

    # Create local vars from configs
    meas_int = configs_dict['measurement_interval']
    min_pct = configs_dict['minimum_pct_annual']

    # Do annual fits for Eo
    Eo_annual_data_dict = {}
    Eo_annual_error_dict = {}
    Eo_pass_keys = []
    Eo_fail_keys = []
    for year in data_dict.keys():
    
        # Calculate number of nocturnal recs for year
        days = 366 if calendar.isleap(year) else 365
        recs = days * (1440 / meas_int) / 2
    
        # Input indices
        data_pct = int(len(data_dict[year]['Fc_series']) / float(recs) * 100)
        if not data_pct < min_pct:
            params, error_code = dark.optimise_all(data_dict[year], 
                                                   params_in_dict)
        else:
            params, error_code = [np.nan, np.nan], 10
        Eo_annual_data_dict[year] = params[0]
        Eo_annual_error_dict[year] = error_code
        if error_code == 0: 
            Eo_pass_keys.append(year)
        else:
            Eo_fail_keys.append(year)
            
    # Fill gaps 
    if np.all(np.isnan(Eo_annual_data_dict.values())):
        print 'Could not find any values of Eo for any years! Exiting...'
        sys.exit()
    if np.any(np.isnan(Eo_annual_data_dict.values())):
        Eo_mean = np.array([Eo_annual_data_dict[year] for year in Eo_pass_keys]).mean()    
        for year in Eo_fail_keys:
            Eo_annual_data_dict[year] = Eo_mean
            Eo_pass_keys.append(year)

    # Attach the yearly Eo to the parameter dictionary
    params_in_dict['Eo_years'] = Eo_annual_data_dict

    # Project Eo to the appropriate indices of the results array
    years_array = np.array([date.year for date in params_out_dict['date']])
    for year in Eo_pass_keys:
        index = np.where(years_array == year)[0]
        out_indices = [index[0], index[-1]]
        params_out_dict['Eo'][out_indices[0]: out_indices[1] + 1] = (
            Eo_annual_data_dict[year])
        params_out_dict['Eo_error_code'][out_indices[0]: out_indices[1] + 1] = (
            Eo_annual_error_dict[year])

    return

def estimate_Re(data_dict, 
                all_params_dict, 
                datetime_input_index_dict):
    
    # Create output arrays for Re estimates
    results_array = np.empty(len(data_dict['TempC']))
    results_array[:] = np.nan
    
    # Estimate time series Re
    for i, date in enumerate(all_params_dict['date']):
        
        indices = datetime_input_index_dict[date]
        this_Eo = all_params_dict['Eo'][i]
        this_rb = all_params_dict['rb'][i]
        this_dict = {'TempC': data_dict['TempC'][indices[0]: indices[1] + 1]}
        results_array[indices[0]: indices[1] + 1] = dark.TRF(this_dict, 
                                                             this_Eo,
                                                             this_rb)
    return {'Re': results_array}

def filtering(this_dict):
    noct_dict = filt.subset_arraydict_on_threshold(this_dict, 'Fsd', 5, '<', 
                                                   drop = True)    
    ustar_dict = filt.subset_arraydict_on_threshold(noct_dict, 'ustar', 0.42, 
                                                    '>', drop = True)
    sub_dict = filt.subset_arraydict_on_nan(ustar_dict)
    return sub_dict

def generate_results_array(datetime_array):
    
    dates_input_index_dict = dtf.get_day_indices(datetime_array)
    date_array = np.array(dates_input_index_dict.keys())
    date_array.sort()
    generic_array = np.empty(len(dates_input_index_dict))
    generic_array[:] = np.nan
    rb_error_code_array = np.ones(len(generic_array)) * 20
    return {'date': date_array,
            'Eo': cp.copy(generic_array),
            'Eo_error_code': cp.copy(generic_array),
            'rb': cp.copy(generic_array),
            'rb_error_code': rb_error_code_array}  

def partition_by_date(data_dict, configs_dict, datetime_array):

    # Create local vars for step and window
    step = configs_dict['respiration_configs']['step_size_days']
    window = configs_dict['respiration_configs']['window_size_days']
    
    # Get the indices of the start and end rows of each year in the source data 
    # array, then build a dict containing a filtered (ustar-screened nocturnal data
    # with nans dropped) data dict for each year with year as key
    years_input_index_dict = dtf.get_year_indices(datetime_array)
    years_data_dict = segment_data(data_dict, years_input_index_dict)
    
    # Get the indices of the start and end rows of each window (depending on step
    # and window size) in the source data array, then build a dict containing a 
    # filtered (ustar-screened nocturnal data with nans dropped) data dict for each 
    # window with the central date for that window as key
    step_dates_input_index_dict = dtf.get_moving_window_indices(datetime_array, 
                                                                window, step)
    step_data_dict = segment_data(data_dict, step_dates_input_index_dict)    
    
    return years_data_dict, step_data_dict
        
def segment_data(data_dict, indices_dict):

    d = {}    
    for key in indices_dict.keys():
        start = indices_dict[key][0]
        end = indices_dict[key][1]
        this_dict = {var: data_dict[var][start: end + 1] 
                     for var in data_dict.keys()}
        d[key] = filtering(this_dict)
    return d        
        
#------------------------------------------------------------------------------
#############
# Main code #
#############        

def main(data_dict, configs_dict):

    # Make local var names for config items
    re_configs_dict = configs_dict['respiration_configs']
    re_configs_dict['measurement_interval'] = (configs_dict['globals']
                                               ['measurement_interval'])
    
    #------------------------------------------------------------------------------
    # Strip datetime from dict
    datetime_array = data_dict.pop('date_time')
    
    # Partition the data into year and step pieces
    years_data_dict, step_data_dict = partition_by_date(data_dict, 
                                                        configs_dict,
                                                        datetime_array)
    
    # Get the indices of the start and end rows of each unique date in the source 
    # data array - no data dict is built from this, since these indices are used to
    # assign Re estimates to the estimated time series output array only
    dates_input_index_dict = dtf.get_day_indices(datetime_array)
    
    # Generate a results dictionary for the parameter values (1 for each day)
    params_out_dict = generate_results_array(datetime_array)
    
    # Initalise parameter dicts with prior estimates
    all_noct_dict = filtering(data_dict)
    params_in_dict = {'Eo_prior': 100,
                      'rb_prior': all_noct_dict['Fc_series'].mean()}
    
    calculate_Eo(years_data_dict, 
                 re_configs_dict,
                 params_in_dict,
                 params_out_dict)
    
    # Get rb for all steps
    calculate_rb(step_data_dict,
                 re_configs_dict,
                 params_in_dict,
                 params_out_dict)
    
    # Estimate Re for all data
    Re_dict = estimate_Re(data_dict,
                          params_out_dict,
                          dates_input_index_dict)
    Re_dict['date_time'] = datetime_array

    return Re_dict, params_out_dict                         