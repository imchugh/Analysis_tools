# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 13:47:20 2015

@author: imchugh
"""

import sys
import numpy as np
import calendar

sys.path.append('../Partitioning')
import gap_filling as gf
import dark_T_response_functions as dark

#------------------------------------------------------------------------------

def calculate_Eo(data_dict,
                 params_in_dict,
                 params_out_dict,
                 years_output_index_dict,
                 meas_int,
                 min_pct):

    # Do annual fits for Eo
    Eo_annual_data_dict = {}
    Eo_annual_error_dict = {}
    Eo_pass_keys = []
    Eo_fail_keys = []
    for year in data_dict.keys():
    
        # Calculate number of nocturnal recs for year
        days = 366 if calendar.isleap(year) else 365
        recs = days * (24 / meas_int) / 2
    
        # Input indices
        data_pct = int(len(data_dict[year]['NEE']) / float(recs) * 100)
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
    for year in Eo_pass_keys:
        out_indices = years_output_index_dict[year]
        params_out_dict['Eo'][out_indices[0]: out_indices[1] + 1] = (
            Eo_annual_data_dict[year])
        params_out_dict['Eo_error_code'][out_indices[0]: out_indices[1] + 1] = (
            Eo_annual_error_dict[year])

    return

def calculate_rb(data_dict,
                 params_in_dict,
                 params_out_dict,
                 meas_int,
                 min_pct):
    
    # Calculate rb for windows    
    for date in data_dict.keys():
        
        # Specify Eo value for the relevant year
        params_in_dict['Eo_default'] = params_in_dict['Eo_years'][date.year]
        data_pct = int(len(data_dict[date]['NEE']) / 
                       float((24 / meas_int) / 2) * 100)
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