# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 16:20:41 2015

@author: imchugh
"""

import numpy as np
import os
import pdb
import datetime as dt
import netCDF4
import xlrd

import random_error_3 as r_e

# Open configuration and build dictionaries of config file contents
def get_configs():
    
    configs_dict = {'measurement_interval': 30,
                    'neg_averaging_bins': 10,
                    'pos_averaging_bins': 10,
                    'radiation_difference_threshold': 35,
                    'temperature_difference_threshold': 3,
                    'windspeed_difference_threshold': 1,
                    'noct_threshold': 5,
                    'ustar_threshold': 0,
                    'num_trials': 10**4,
                    'nan_value': -9999,
                    'QC_accept_code': 0,
                    'mean_flux_series': 'Fc_model'}
                    
    return configs_dict

def get_data(configs_dict):

    # Specify file location and name
    data_input_path = '/home/imchugh/Ozflux/Sites/Whroo/Data/Processed/all'
    data_input_file = 'Whroo_2011_to_2014_L5.nc'
    data_input_target = os.path.join(data_input_path, data_input_file)

    # Initialise dictionaries
    vars_dict = {'carbon_flux':'Fc',
                 'modelled_carbon_flux': 'Fc_SOLO',
                 'solar_radiation':'Fsd',
                 'temperature':'Ta',
                 'wind_speed': 'Ws_CSAT',
                 'friction_velocity': 'ustar'}
    QC_dict = {'carbon_flux':'Fc_QCFlag',
               'solar_radiation':'Fsd_QCFlag',
               'temperature':'Ta_QCFlag',
               'wind_speed': 'Ws_CSAT_QCFlag',
               'friction_velocity': 'ustar_QCFlag'}
    newNames_dict = {'carbon_flux':'Fc',
                     'modelled_carbon_flux': 'Fc_model',
                     'solar_radiation':'Fsd',
                     'temperature':'Ta',
                     'wind_speed': 'ws',
                     'friction_velocity': 'ustar'}

    # Read .nc file
    d={}
    vars_list = vars_dict.values() + QC_dict.values()
    nc_obj = netCDF4.Dataset(data_input_target)
    date_time = np.array([dt.datetime(*xlrd.xldate_as_tuple(elem,0)) 
                          for elem in nc_obj.variables['xlDateTime']])
    for i in vars_list:
        ndims=len(nc_obj.variables[i].shape)
        if ndims==3:
            d[i]=nc_obj.variables[i][:,0,0]
        elif ndims==1:    
            d[i]=nc_obj.variables[i][:]
        d[i] = np.where(d[i] == configs_dict['nan_value'], np.nan, d[i])
    nc_obj.close()

    # Screen low ustar values
    index = np.where(d[vars_dict['friction_velocity']] < 
                     configs_dict['ustar_threshold'])
    d[vars_dict['carbon_flux']][index] = np.nan

    # Replace configured error values with NaNs and remove data with unacceptable 
    # QC codes, then drop QC flag variables and rename variable to new names
    for key in vars_dict.keys():
        if key in QC_dict.keys():
            d[vars_dict[key]] = np.where(d[QC_dict[key]] != 
                                             configs_dict['QC_accept_code'],
                                         np.nan, d[vars_dict[key]])
            d.pop(QC_dict[key])
        d[newNames_dict[key]] = d.pop(vars_dict[key])

    # Add date_time variable
    d['date_time'] = date_time
          
    return d

def main():    

    results_out_path = '/home/imchugh/Documents'

    # Update
    reload(r_e)

    # Get configurations and data
    configs_dict = get_configs()
    data_dict = get_data(configs_dict)

    # Calculate sigma_delta as a function of flux magnitude
    fig, stats_dict = r_e.calculate_sigma_delta(data_dict, configs_dict)    
    fig.savefig(os.path.join(results_out_path, 'Random_error_plots.jpg'))    
    sig_del = r_e.calculate_sigma_delta(data_dict['Fc'], stats_dict)
    
    return sig_del
    
    






