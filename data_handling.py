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
    """
    Pass: 1) data_dict - a dictionary containing arrays, one of which must be a 
             python datetime;
          2) date_time_var - namestring of datetime variable
          3) year to be returned
    Returns: if year is specified, return the same dictionary structure with 
             only data for that year; if no year is specified, it returns a 
             dictionary with each data year contained as the value with the 
             year as the key
    """    
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
    """
    Pass a dictionary
    """
    hi_index = data_dict[ind_var] > threshold
    lo_index = data_dict[ind_var] < threshold
    return {'>': subset_data_dict(data_dict, hi_index),
            '<': subset_data_dict(data_dict, lo_index)}

def sort_dict_on_index_variable(data_dict, sort_var):
    """
    Sorts a dictionary of equal-length numpy arrays on the basis of a
    sorting variable
    Pass: 1) data_dict - a dictionary containing arrays
          2) sort_var - the variable whose sorted order will dictate the 
                        ordering of all others
    """
    index = sort_var.argsort()
    for key in data_dict.keys():
        data_dict[key] = data_dict[key][index]
    return data_dict
            
def set_arraydict_to_nan_conditional(data_dict, ind_var, value, 
                                     cond = 'less', set_vars = None):
    if cond == 'less':
        nan_index = data_dict[ind_var] < value
    elif cond == 'greater':
        nan_index = data_dict[ind_var] > value
    elif cond == 'equal':
        nan_index = data_dict[ind_var] == value
    else: return
    if set_vars:
        if not isinstance(set_vars, list): set_vars = [set_vars]
    else:
        set_vars = data_dict.keys()
    new_dict = data_dict.copy()
    for var in set_vars:
        new_dict[var][nan_index] = np.nan
    return new_dict
    