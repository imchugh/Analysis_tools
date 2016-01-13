# -*- coding: utf-8 -*-
"""
Created on Sat Jul 11 12:18:25 2015

@author: imchugh
"""

from scipy.optimize import curve_fit
import numpy as np
import pdb

#------------------------------------------------------------------------------
# Data optimisation algorithm

def TRF(data_dict, Eo, rb):
    return rb * np.exp(Eo * (1 / (10 + 46.02) - 1 / (data_dict['TempC'] + 46.02)))

#------------------------------------------------------------------------------    
# Write error messages to dictionary with codes as keys
def error_codes():
    
    d = {0:'Optimisation successful',
         1:'Value of Eo failed range check - rejecting all parameters',
         2:'Value of rb has wrong sign - rejecting all parameters',
         3:'Optimisation reached maximum number of iterations' \
           'without convergence',
         10:'Data did not pass minimum percentage threshold - ' \
            'skipping optimisation'}
    
    return d
#------------------------------------------------------------------------------    

# rb and Eo
def optimise_all(data_dict, params_dict):

    # Initialise error state variable
    error_state = 0              
    drivers_dict = {driver: data_dict[driver] for driver in ['TempC']}
    response_array = data_dict['NEE_series']

    try:
        params = curve_fit(lambda x, a, b:
                           TRF(x, a, b),
                           drivers_dict, 
                           response_array, 
                           p0 = [params_dict['Eo_prior'], 
                                 params_dict['rb_prior']])[0]
    except RuntimeError:
        params = [np.nan, np.nan]
        error_state = 3

    # If negative rb returned, set to nan
    if params[0] < 50 or params[0] > 400: 
        error_state = 1
        params = [np.nan, np.nan]
    elif params[1] < 0:
        error_state = 2
        params = [np.nan, np.nan]

    return {'Eo': params[0], 'rb': params[1], 'error_code': error_state}

# rb   
def optimise_rb(data_dict, params_dict):

    # Initialise error state variable
    error_state = 0              
    
    # Get drivers and response
    drivers_dict = {driver: data_dict[driver] for driver in ['TempC']}
    response_array = data_dict['NEE_series']        
    
    try:
        params = curve_fit(lambda x, b:
                           TRF(x, params_dict['Eo_default'], b),
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
