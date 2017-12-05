#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 17:07:47 2017

@author: ian
"""

import calendar
import copy as cp
import datetime as dt
from lmfit import Model
import numpy as np
import pandas as pd
import sys
import pdb

import DataIO as io

###############################################################################
# Functions                                                                   #
###############################################################################

#------------------------------------------------------------------------------
def get_slice(df, configs_dict, date_obj):
    
    interval_mins = int(filter(str.isdigit, configs_dict['interval']))
    if isinstance(date_obj, tuple):
        pct_required = configs_dict['min_data_pct_window']
        expected_length = 1440 / interval_mins * configs_dict['window_size']
        sub_df = cp.copy(df.loc[date_obj[0]:date_obj[1]])
    elif isinstance(date_obj, str):
        pct_required = configs_dict['min_data_pct_annual']
        n_days = 366 if calendar.isleap(int(date_obj)) else 365
        expected_length = 1440 / interval_mins * n_days
        sub_df = cp.copy(df.loc[date_obj])
    sub_df = sub_df.loc[sub_df.Fsd < 5, ['Fc', 'Ts', 'Sws']].dropna()
    pct_available = len(sub_df) / float(expected_length) * 100
    assert pct_available > pct_required
    return sub_df
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def generate_data(obs_df, params_df):
    
    df_list = []
    for this_date in params_df.index:
        date_str = dt.datetime.strftime(this_date, '%Y-%m-%d')
        sub_df = obs_df.loc[date_str]
        param_list = map(lambda x: params_df.loc[date_str, x], 
                         ['rb', 'Eo', 'theta_1', 'theta_2'])
        df_list.append(response_func(sub_df.Ts, sub_df.Sws, 
                                     param_list[0], param_list[1],
                                     param_list[2], param_list[3]))
    output_df = pd.DataFrame(pd.concat(df_list), columns = ['ER'])
    return output_df.reindex(obs_df.index).interpolate()
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def get_data(configs_dict):

    df = io.OzFluxQCnc_to_data_structure(configs_dict['file_path'], 
                                         output_structure='pandas') 
    try:
        freq = pd.infer_freq(df.index)
        assert freq == '30T' or '60T'
    except AssertionError:
        print 'File is not chronologically continuous, exiting...'
        sys.exit()
    configs_dict['interval'] = freq
    return df
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def make_date_iterator(df, configs_dict):
    
    size = configs_dict['window_size']
    step = configs_dict['window_step']
    interval_mins = int(filter(str.isdigit, configs_dict['interval']))
    start_date = df.index[0].to_pydatetime().date() + dt.timedelta(size / 2)
    end_date = df.index[-1].to_pydatetime().date()
    freq_str = '{}D'.format(str(int(step)))
    dates = pd.date_range(start_date, end_date, freq = freq_str)
    ref_dates = dates + dt.timedelta(0.5)
    date_bounds = zip(ref_dates - dt.timedelta(size / 2.0 - interval_mins / 1440.0),
                      ref_dates + dt.timedelta(size / 2.0))
    return dict(zip(dates, date_bounds))
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def make_params_df(df):
    return pd.DataFrame(index = pd.date_range(df.index[0].date(), 
                                              df.index[-1].date(), freq = 'D'),
                        columns = ['rb', 'Eo', 'theta_1', 'theta_2', 'QCFlag'],
                        dtype = 'float')
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def optimise(df, model, params):
    result = model.fit(df.Fc, t_series = df.Ts, vwc_series = df.Sws, 
                       params = params)
#    if params['Eo'] < 50 or params['Eo'] > 400:
        
    return result
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def process_data(df, configs_dict):

    # Generate a date iterator based on window size and step
    date_iterator_dict = make_date_iterator(df, configs_dict)
    dates_list = sorted(date_iterator_dict.keys())

    # Make a result dataframe to hold parameter values
    params_df = make_params_df(df)

    # Initialise the model
    model = Model(response_func, independent_vars = ['t_series', 
                                                     'vwc_series'])
    params = model.make_params(rb = 1, Eo = 200, theta_1 = 1, theta_2 = 10)
    
    # Iterate on year
    years_list = sorted(map(lambda x: str(x), list(set(df.index.year))))
    for this_year in years_list:
        try:
            sub_df = get_slice(df, configs_dict, this_year)
        except AssertionError:
            params_df.loc[this_year, 'QCFlag'] = 1
            continue
        result = optimise(sub_df, model, params)
        for param in result.best_values.keys():
            if not param == 'rb':
                params_df.loc[str(this_year), param] = result.best_values[param]
   
    # Fix parameters
    params['Eo'].vary = False
    params['theta_1'].vary = False
    params['theta_2'].vary = False
    
    # Iterate on date
    for this_date in dates_list:
        try:
            sub_df = get_slice(df, configs_dict, date_iterator_dict[this_date])
        except AssertionError:
            params_df.loc[this_date, 'QCFlag'] = 2
            continue
        result = optimise(sub_df, model, params)
        params_df.loc[this_date, 'rb'] = result.best_values['rb']
    
    return params_df.interpolate()
#------------------------------------------------------------------------------    

#------------------------------------------------------------------------------
def response_func(t_series, vwc_series, rb, Eo, theta_1, theta_2):
    
    return (rb  * np.exp(Eo * (1 / (10 + 46.02) - 1 / (t_series + 46.02))) *
            1 / (1 + np.exp(theta_1 - theta_2 * vwc_series)))
#------------------------------------------------------------------------------
    
###############################################################################
# Main program                                                                #
###############################################################################
    
# Set window size and step (both in units of days)
configs_dict = {'file_path': ('/home/ian/OzFlux/Sites/GatumPasture/Data/'
                              'Processed/All/GatumPasture_L4.nc'),
                'window_size': 7, 
                'window_step': 5,
                'min_data_pct_window': 20,
                'min_data_pct_annual': 10}

# Get the data and check timestamp integrity (exit if not continuous)
df = get_data(configs_dict)

# Get parameters
params_df = process_data(df, configs_dict)

# Calculate respiration for time series
result_df = generate_data(df, params_df)
