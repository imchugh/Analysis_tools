#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 10:07:35 2017

@author: ian
"""
import pandas as pd
import datetime as dt
import profile_data_processing as pdp
import os
import pdb

#------------------------------------------------------------------------------
# Common scripts
#------------------------------------------------------------------------------

def get_site_data(site_name):
    
    sites_dict = {'howard_springs': howard_springs,
                  'warra_avg': warra_average,
                  'warra_raw': warra_raw} 
    
    return sites_dict[site_name]()

def write_data_to_file(df, path, site):
    
    if not os.path.isdir(path):
        home_dir = os.path.expanduser('~')
        output_dir = os.path.join(home_dir, 'profile_data')
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        print ('The specified directory for data output could not be found!'
               '\nWriting to the following directory: {0}'
               .format(output_dir))
    else:
        output_dir = path

    fname = '{0}_standard_format_profile.csv'.format(site)
    target = os.path.join(output_dir, fname)    
     
    df.to_csv(target, index_label = 'Datetime')
    
    return

#------------------------------------------------------------------------------    
# Site scripts
#------------------------------------------------------------------------------

###############################################################################
# Warra
###############################################################################

def warra_raw(write_to_dir = None):

    ###########################################################################
    # User-setable options - set with CARE!!!
    # 1) Variable setting the number of seconds to drop to allow for manifold 
    #    flush
    drop_n_leading_seconds = 5
    # 2) List of profile valve numbers
    valve_list = [1, 2, 3, 4, 5, 6, 7, 8]
    # 3) List of profile heights to cross-match with valve numbers
    heights_list = [2, 4, 8, 16, 30, 42, 54, 70]
    ###########################################################################

    # Open files and concatenate (note that the file configuration
    # currently has files starting on day x 00:00:00.5 and ending on 
    # day x+1 00:00:00; this is a problem because the switch to the next valve 
    # occurs after 00:00:00; this means that the first set of measurements 
    # on valve 1 only has 14 instead of 15 instances, and there is a single 
    # observation on valve 1 at the end of the file; so we just move the time 
    # stamp by 0.5s)
    dir_str = pdp.dir_select_dialog()        
    dir_list = os.listdir(dir_str)
    df_list = []
    for f in dir_list:
        full_path = os.path.join(dir_str, f)
        df_list.append(pd.read_csv(full_path, skiprows = [0, 2, 3]))
    df = pd.concat(df_list)
    df.index = pd.to_datetime(df.TIMESTAMP)
    df.sort_index(inplace = True)
    df.index = df.index + dt.timedelta(seconds = 0.5)
    df = df.reindex(df.index - dt.timedelta(seconds = 0.5))
    df['modulo_15'] = df.index.second % 15

    # Generate the required date ranges for indexing and outputting data
    start = dt.datetime(df.index[0].year, df.index[0].month, df.index[0].day,
                        0, 0, 0, 500000)
    end = dt.datetime(df.index[-1].year, df.index[-1].month, df.index[-1].day)
    dt_range_1 = pd.date_range(start, end, freq = '2T')
    dt_range_2 = dt_range_1 + dt.timedelta(minutes = 1, seconds = 59, 
                                           microseconds = 500000)
    dt_range_out = dt_range_1 + dt.timedelta(minutes = 1, seconds = 59, 
                                             microseconds = 500000)
    
    # Make a reference dictionary cross-matching valve number with height
    str_heights_list = [str(height) for height in heights_list]
    
    # Make a CO2 names list and a reference dictionary for assigning the data 
    # to the output dataframe
    CO2_names_list = ['CO2_{0}m'.format(height) 
                      for height in str_heights_list]    
    CO2_names_dict = dict(zip(valve_list, CO2_names_list))
    
    # Make a T names list and a reference dictionary for assigning the data 
    # to the output dataframe
    T_names_list = ['Tair_{0}m'.format(height) 
                    for height in str_heights_list] 
    T_names_dict = dict(zip(valve_list, T_names_list))
    
    # Make an output dataframe
    rslt_df = pd.DataFrame(index = dt_range_out, 
                           columns = CO2_names_list + T_names_list)
    
    # Cycle through all time periods
    for i in xrange(len(dt_range_1)):
    
        # Subset the dataframe to get the required period (cut some of the data
        # to allow for flushing of manifold)
        this_df = df.loc[dt_range_1[i]: dt_range_2[i]]
        sub_df = (this_df[this_df.modulo_15 >= drop_n_leading_seconds]
                  .groupby('valve_number').mean())
    
        # Do the averaging for each valve and send CO2 and temp data to output 
        # df
        for valve in valve_list:
            try:
                rslt_df.loc[dt_range_out[i], 
                            CO2_names_dict[valve]] = sub_df.loc[valve, 'CO2']
            except:
                continue
            T_name = 'T_air({0})'.format(str(valve))
            rslt_df.loc[dt_range_out[i], 
                        T_names_dict[valve]] = sub_df.loc[valve, T_name]
    
    if not write_to_dir is None:
        write_data_to_file(rslt_df, write_to_dir, 'Warra')
        
    return rslt_df

#------------------------------------------------------------------------------
        
def warra_average(write_to_dir = None):
    
    # Create a dict to reference heights to valve numbers
    profile_n = [1, 2, 3, 4, 5, 6, 7, 8]
    profile_heights = [2, 4, 8, 16, 30, 42, 54, 70]
    heights_dict = dict(zip(profile_n, 
                                           [str(height) for height in 
                                            profile_heights]))
    
    # Create a dict to correct the lag due to valve cycling
    lag_dict = {1: 105,
                2: 90,
                3: 75,
                4: 60,
                5: 45,
                6: 30,
                7: 15,
                8: 0}

    # Prepare df
    file_in = pdp.file_select_dialog()
    df = pd.read_csv(file_in, skiprows = [0, 2, 3])
    df.index = pd.to_datetime(df.TIMESTAMP)
    idx = df[df.valve_number == 8].index
    rslt_df = pd.DataFrame(index = idx)
    
    # Cycle through time series and break out CO2 variable to individual 
    # heights on basis of valve number
    T_names = [var for var in df.columns if 'T_air' in var]
    for i in profile_n:
        sub_df = df[df.valve_number == i]
        sub_df.index = sub_df.index + dt.timedelta(seconds = lag_dict[i])
        sub_df = sub_df.reindex(idx)
        CO2_out_name = 'CO2_{0}m'.format(heights_dict[i])
        T_in_name = [var for var in T_names if str(i) in var]
        T_out_name = 'Tair_{0}m'.format(heights_dict[i])
        rslt_df[CO2_out_name] = sub_df['CO2_Avg']
        rslt_df[T_out_name] = sub_df[T_in_name]
        
    return rslt_df

###############################################################################



###############################################################################
# Howard Springs                                                              #
###############################################################################

def howard_springs(write_to_dir = None):

    file_str = '/home/ian/OzFlux/Sites/Howard_Springs/Data/Profile/Howard_profile_Slow_avg2.dat'
    #pdp.dir_select_dialog()        

    df = pd.read_csv(file_str, skiprows = [0, 2, 3])
    df.drop(df.index[1], inplace = True)
    df.index = pd.to_datetime(df.TIMESTAMP)
    df.drop_duplicates(inplace = True)

    old_names = [i for i in df.columns if 'Cc' in i]
    old_names.sort()

    new_names = ['CO2_{0}'.format(this_name.split('_')[2]) 
                 for this_name in old_names]
    
    old_names.append('T_air_Avg')
    new_names.append('Tair_2m')
    
    names_dict = dict(zip(old_names, new_names))
    new_df = pd.DataFrame(index = df.index)
    for name in names_dict.keys():
        new_df[names_dict[name]] = df[name]

    return new_df
