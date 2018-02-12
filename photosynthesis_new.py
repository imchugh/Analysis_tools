#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 13:01:47 2018

@author: ian
"""
from lmfit import Model
import numpy as np
import pandas as pd

import respiration_new as rn
import DataIO as io

#------------------------------------------------------------------------------
def response_func(fsd_series, vpd_series, alpha, beta, k):
    beta_VPD = beta * np.exp(-k * (vpd_series - 1))
    index = vpd_series <= 1
    beta_VPD[index] = beta
    GPP = (alpha * data_d['PAR']) / (1 - (data_d['PAR'] / 2000) + 
           (alpha * data_d['PAR'] / beta_VPD))
    index = data_d['PAR'] < 5
    GPP[index] = 0          
    pass

#------------------------------------------------------------------------------
    
###############################################################################
# Main program                                                                #
###############################################################################
def main(df, config_dict = None):
    
    if config_dict is None:      
        configs_dict = {'file_path': ('/home/ian/OzFlux/Sites/GatumPasture/'
                                      'Data/Processed/All/GatumPasture_L4.nc'),
                        'window_size': 10, 
                        'window_step': 5,
                        'min_data_pct_window': 20,
                        'min_data_pct_annual': 10}
    
    pass











    


path = '/home/ian/OzFlux/Sites/GatumPasture/Data/Processed/All/GatumPasture_L5.nc'

df = io.OzFluxQCnc_to_data_structure(path, output_structure = 'pandas')


#
#
## Get the data and check timestamp integrity (exit if not continuous)
#df = get_data(configs_dict)
#
## Get parameters
#params_df = process_data(df, configs_dict)
#
## Calculate respiration for time series
#result_df = generate_data(df, params_df)