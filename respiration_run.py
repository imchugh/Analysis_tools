# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 16:08:41 2016

@author: imchugh
"""

# Python standard modules
import os
import copy as cp
import numpy as np

# My modules
import DataIO as io
import respiration as re
import data_formatting as dt_fm

reload(re)

# Get the data and format appropriately
def get_data(configs_dict):

    # Get file extension and target
    paths_dict = configs_dict['files']
    ext = os.path.splitext(paths_dict['input_file'])[1]
    data_input_target = os.path.join(paths_dict['input_path'],
                                     paths_dict['input_file'])

    # Initialise name change dictionary with new names via common keys
    oldNames_dict = configs_dict['variables']
    newNames_dict = {'carbon_flux':'Fc_series',
                     'temperature': 'TempC',
                     'solar_radiation': 'Fsd',
                     'vapour_pressure_deficit': 'VPD',
                     'friction_velocity': 'ustar',
                     'wind_speed': 'ws'}
    names_dict = {oldNames_dict[key]: newNames_dict[key] for key in oldNames_dict}                     

    # get data (screen only the Fc data to obs only)
    if ext == '.nc':
        Fc_dict = io.OzFluxQCnc_to_data_structure(data_input_target,
                                                  var_list = [oldNames_dict
                                                              ['carbon_flux']],
                                                  QC_accept_codes = [0])
        Fc_dict.pop('date_time')
        ancillary_vars = [oldNames_dict[var] for var in oldNames_dict.keys() 
                          if not var == 'carbon_flux']
        ancillary_dict, global_attr = io.OzFluxQCnc_to_data_structure(
                                          data_input_target,
                                          var_list = ancillary_vars,
                                          return_global_attr = True)
        data_dict = dict(Fc_dict, **ancillary_dict)
    elif ext == '.df':
        data_dict, global_attr = io.DINGO_df_to_data_structure(
                                     data_input_target,
                                     var_list = oldNames_dict.values(),
                                     return_global_attr = True)
    
    # Rename relevant variables    
    data_dict = dt_fm.rename_data_dict_vars(data_dict, names_dict)

    return data_dict, global_attr
    
#------------------------------------------------------------------------------    
# Do the respiration fit

# Get configurations
configs_dict = io.config_to_dict(io.file_select_dialog())

# Get data
data_dict, attr = get_data(configs_dict)

# Drop additional config items and add measurement interval to configs
configs_dict = configs_dict['respiration_configs']
configs_dict['measurement_interval'] = int(attr['time_step'])

# Remove low ustar values according to threshold
data_dict['Fc_series'][data_dict['ustar'] < 
                       configs_dict['ustar_threshold']] = np.nan

# Calculate Re by sending data to main respiration function
re_dict, params_dict = re.main(cp.copy(data_dict), configs_dict)
data_dict['Re'] = re_dict['Re']