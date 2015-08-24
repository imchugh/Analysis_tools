# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 12:13:41 2015

@author: imchugh
"""

import numpy as np
import pdb

def subset_data_dict(data_dict, boolean_index, copy = True):
    
    work_dict = data_dict.copy() if copy else data_dict
    return {var: data_dict[var][boolean_index] for var in work_dict.keys()}

def subset_datayear_from_arraydict(data_dict, date_time_var, year = None):
    
    years_array = np.array([date_.year for date_ in data_dict[date_time_var]])
    if not year:
        year_list = set(list(years_array))    
    else:
        if not isinstance(year, list): year = [year]
    
    new_dict = {}
    for yr in year_list:
        year_index = years_array == yr            
        new_dict[yr] = subset_data_dict(data_dict, year_index)

    return new_dict
    
def subset_onthreshold_from_arraydict(data_dict, ind_var, threshold):
    
    hi_index = data_dict[ind_var] > threshold
    lo_index = data_dict[ind_var] < threshold
    return {'>': subset_data_dict(data_dict, hi_index),
            '<': subset_data_dict(data_dict, lo_index)}
            
def set_arraydict_to_nan_conditional(data_dict, ind_var, threshold, hilo = '>'):
    
    if hilo == '>':
        nan_index = data_dict[ind_var] < threshold
    elif hilo == '<':
        nan_index = data_dict[ind_var] > threshold
    else: return
    new_dict = data_dict.copy()
    for var in new_dict.keys():
        new_dict[var][nan_index] == np.nan
    return new_dict
    