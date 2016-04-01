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
import light_response as li
import data_formatting as dt_fm
import data_filtering as data_filter

#------------------------------------------------------------------------------
# Fetch data from configurations
def get_data(configs_dict):

    # Get data (screen Fc data to obs only - keep gap-filled drivers etc)
    data_input_target = os.path.join(configs_dict['files']['input_path'],
                                     configs_dict['files']['input_file'])
    Fc_dict = io.OzFluxQCnc_to_data_structure(data_input_target,
                                              var_list = [configs_dict['variables']
                                                                      ['carbon_flux']],
                                              QC_accept_codes = [0])
    Fc_dict.pop('date_time')
    ancillary_vars = [configs_dict['variables'][var] for var in 
                      configs_dict['variables'] if not var == 'carbon_flux']
    ancillary_dict, attr = io.OzFluxQCnc_to_data_structure(
                               data_input_target,
                               var_list = ancillary_vars,
                               return_global_attr = True)
    data_dict = dict(Fc_dict, **ancillary_dict)
    
    # Rename to generic names used by scripts
    old_names_dict = configs_dict['variables']
    std_names_dict = dt_fm.standard_names_dictionary()
    map_names_dict = {old_names_dict[key]: std_names_dict[key] 
                      for key in old_names_dict}
    data_dict = dt_fm.rename_data_dict_vars(data_dict, map_names_dict)
    
    return data_dict, attr
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Rebuild the master configuration file for passing to respiration and light 
# response (if requested) functions
def build_config_file(configs_master_dict, do_light_response):
    
    # Build a specific configuration file
    configs_dict = {'files': configs_master_dict['global_configs']['files'],
                    'global_options': (configs_master_dict['global_configs']
                                                          ['options'])}                                                          
    if do_light_response:
        configs_dict['variables'] = dict(configs_master_dict
                                         ['respiration_configs']['variables'],
                                         ** configs_master_dict
                                            ['light_response_configs']
                                            ['variables'])
    else:
        configs_dict['variables'] = (configs_master_dict['respiration_configs']
                                                        ['variables'])
    return configs_dict                                                         
#------------------------------------------------------------------------------
    
#------------------------------------------------------------------------------
# Remove low ustar values
def screen_low_ustar(data_dict, configs_dict):
    
    ustar_threshold = configs_dict['global_options']['ustar_threshold']
    noct_threshold = configs_dict['global_options']['noct_threshold']
    if isinstance(ustar_threshold, dict):
        years_data_dict = data_filter.subset_datayear_from_arraydict(data_dict, 
                                                                     'date_time')
        threshold_keys = [int(key) for key in ustar_threshold.keys()]
        miss_list = [year for year in years_data_dict.keys() 
                     if not year in threshold_keys]
        if not len(miss_list) == 0:
            miss_string = ', '.join([str(this_year) for this_year in miss_list])
            raise Exception('Missing years: %s' %miss_string + '; please edit ' \
                            'your configuration file so that years specified ' \
                            'for ustar threshold match those available in ' \
                            'data file, or alternatively specify a single ' \
                            'value (float) in the configuration file under ' \
                            '[global_configs][options][ustar_threshold]. '\
                            'Exiting...')
        data_list = []
        for this_year in years_data_dict.keys():
            this_threshold = ustar_threshold[str(this_year)]
            this_NEE = years_data_dict[this_year]['NEE_series']
            this_NEE[(years_data_dict[this_year]['ustar'] < this_threshold) &
                     (years_data_dict[this_year]['Fsd'] < noct_threshold)] = np.nan
            data_list.append(this_NEE)
        data_dict['NEE_series'] = np.concatenate(data_list)
    else:
        data_dict['NEE_series'][(data_dict['ustar'] < ustar_threshold) &
                                (data_dict['Fsd'] < noct_threshold)] = np.nan
                                
    return
#------------------------------------------------------------------------------    

#------------------------------------------------------------------------------    
def main(use_storage = False, ustar_threshold = False, 
         config_file = False, do_light_response = False):
    """
    No positional arguments - prompts for a configuration file
    Kwargs: use_storage - if True then algorithm looks for a variable called 
                          Fc_storage and then sums that with Fc
            ustar_threshold - set to a particular value to override the ustar
                              threshold set in the global configs root item of
                              the configuration file (this is done so that can
                              be set in function call from another script)
    """
    
    # Do the respiration fit

    # Get master configuration file
    if not config_file:
        configs_master_dict = io.config_to_dict(io.file_select_dialog())
    else:
        configs_master_dict = io.config_to_dict(config_file)
    
    # Build custom configuration file for this script
    configs_dict = build_config_file(configs_master_dict, do_light_response)

    # Get data
    data_dict, attr = get_data(configs_dict)

    # Override default configuration file ustar_threshold if requested by user
    if not isinstance(ustar_threshold, bool):
        if isinstance(ustar_threshold, (int, float, dict)):
            configs_dict['global_configs']['ustar_threshold'] = ustar_threshold

    # Sum Fc and Sc if storage is to be included, otherwise if requested, 
    # remove all Fc where Sc is missing
    if configs_dict['global_options']['use_storage']:
        data_dict['NEE_series'] = (data_dict['NEE_series'] + 
                                   data_dict['Fc_storage'])
    elif configs_dict['global_options']['unify_flux_storage_cases']:
        data_dict['NEE_series'][np.isnan(data_dict['Fc_storage'])] = np.nan

    # Remove low ustar data
    screen_low_ustar(data_dict, configs_dict)     
    
    # Set up respiration configs and add measurement interval and output path
    re_configs_dict = configs_master_dict['respiration_configs']['options']
    re_configs_dict['measurement_interval'] = int(attr['time_step'])
    re_full_path = os.path.join(configs_dict['files']['output_path'],
                                re_configs_dict['output_folder'])
    if not os.path.isdir(re_full_path): os.makedirs(re_full_path)
    re_configs_dict['output_path'] = re_full_path
    
    # Calculate Re                            
    re_rslt_dict, re_params_dict, re_error_dict = re.main(cp.copy(data_dict), 
                                                          re_configs_dict)
    
    # Add original time series                                                              
    re_rslt_dict['NEE_series'] = data_dict['NEE_series']
    
    # Do the light response parameters (or not)
    if do_light_response:
    
        # Set up light response configs and add measurement interval and output path
        li_configs_dict = configs_master_dict['light_response_configs']['options']
        li_configs_dict['measurement_interval'] = int(attr['time_step'])
        li_full_path = os.path.join(configs_dict['files']['output_path'],
                                    li_configs_dict['output_folder'])
        if not os.path.isdir(li_full_path): os.makedirs(li_full_path)
        li_configs_dict['output_path'] = li_full_path
        
        # Convert insolation to PPFD
        data_dict['PAR'] = data_dict['Fsd'] * 0.46 * 4.6
        
        # Call light response function
        li_rslt_dict, li_params_dict, li_error_dict = li.main(data_dict, 
                                                              li_configs_dict, 
                                                              re_params_dict)

        # Add original time series                                                              
        li_rslt_dict['NEE_series'] = data_dict['NEE_series']
        
        # Write data to file
        var_order_list = ['date', 'Eo', 'Eo_error_code', 'rb', 'alpha', 'beta', 'k', 
                          'light_response_error_code']
        if li_configs_dict['use_nocturnal_rb']:
            var_order_list.insert(4, 'rb_error_code')
        io.array_dict_to_csv(li_params_dict, 
                             os.path.join(li_full_path, 'params.csv'), 
                             var_order_list)
        io.array_dict_to_csv(li_rslt_dict, os.path.join(li_full_path, 'Re_GPP.csv'), 
                             ['date_time', 'Re', 'GPP'])
        io.text_dict_to_text_file(li_error_dict, os.path.join(li_full_path, 'error_codes.txt'))

        return li_rslt_dict, li_params_dict
        
    else:
        
        # Write data to file
        io.array_dict_to_csv(re_params_dict, 
                             os.path.join(re_full_path, 'params.csv'), 
                             ['date', 'Eo', 'Eo_error_code', 'rb', 'rb_error_code'])
        io.array_dict_to_csv(re_rslt_dict, os.path.join(re_full_path, 'Re.csv'), 
                             ['date_time', 'Re'])
        io.text_dict_to_text_file(re_error_dict, os.path.join(re_full_path, 'error_codes.txt'))        
        
        return re_rslt_dict, re_params_dict

if __name__ == "__main__":

   main()