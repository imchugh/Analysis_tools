# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 14:37:38 2015

@author: imchugh
"""
import pdb
import numpy as np
import sys
sys.path.append('../Partitioning')
import Partition_NEE_5 as pt

def daytime_model_error(data_dict, configs_dict):

    # Create data arrays with nocturnal data removed, then count all
    day_index = data_dict['Fsd'] > configs_dict['noct_threshold']
    obs_array = data_dict['Fc'][day_index]
    mod_array = data_dict['Fc_model'][day_index]
    total_records = len(obs_array)

    # Rescale to gC m-2
    for arr in [obs_array, mod_array]:
        arr[...] = arr * configs_dict['measurement_interval'] * 60 * 12 * 10 ** -6

    # Calculate annual sum for obs and model combined
    day_sum = np.where(np.isnan(obs_array), mod_array, obs_array).sum()

    # Subset arrays to remove nans and calculate proportion remaining
    nan_index = ~np.isnan(obs_array)
    obs_array = obs_array[nan_index]
    mod_array = mod_array[nan_index]
    avail_records = len(obs_array)
    
    # Get the amount of data that will be removed from the sample (based on the 
    # proportion missing from the complete dataset)
    sample_missing = 1000 - int(1000 * (avail_records / float(total_records)))
    
    # Draw a random sample of 1000 data from the timeseries, then calculate the
    # difference between the observed and model-spliced series (appropriately 
    # scaled to gC m-2)
    error_array = np.empty(configs_dict['num_trials'])
    for this_trial in xrange(configs_dict['num_trials']):
        random_index = np.random.randint(0, len(obs_array), 1000)
        subset_obs_array = obs_array[random_index]
        subset_mod_array = mod_array[random_index]
        subset_splice_array = np.concatenate([subset_mod_array[:sample_missing],
                                              subset_obs_array[sample_missing:]])
        obs_sum = subset_obs_array.sum()
        splice_sum = subset_splice_array.sum()
        error_array[this_trial] = (obs_sum - splice_sum) / obs_sum
        
    propn_error = error_array.std() * 2
    abs_error = abs(day_sum * propn_error)
                                   
    return abs_error
    
def nocturnal_model_error(data_dict, configs_dict):
    
    reload(pt)    
    
    # Return parameter and series dictionaries
    param_dict, series_dict = pt.main(data_dict, configs_dict)
    
    return series_dict