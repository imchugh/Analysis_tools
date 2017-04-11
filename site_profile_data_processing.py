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

def get_site_data(site_name):
    
    sites_dict = {'warra_avg': warra_average,
                  'warra_raw': warra_raw} 
    
    return sites_dict[site_name]()

# User configurations
path = '/home/ian/Downloads/TOA5_RawData392.dat'

def warra_raw():

    # Trim n leading seconds to allow manifold flush
    drop_n_leading_seconds = 5

    profile_n = [1, 2, 3, 4, 5, 6, 7, 8]
    profile_heights = [2, 4, 8, 16, 30, 42, 54, 70]
    heights_dict = heights_dict = dict(zip(profile_n, 
                                           [str(height) for height in 
                                            profile_heights]))

#    dir_in = pdp.dir_select_dialog()
#    pdb.set_trace()
#    return dir_in
        
    file_in = '/home/ian/OzFlux/Sites/Warra/Data/Profile/Raw/TOA5_RawData391.dat'
    dir_str = '/home/ian/OzFlux/Sites/Warra/Data/Profile/Raw'
    dir_list = os.listdir(dir_str)
    file_in = os.path.join('/home/ian/OzFlux/Sites/Warra/Data/Profile/Raw', dir_list[0])
    df_list = []
    for f in dir_list:
        full_path = os.path.join(dir_str, f)
        df_list.append(pd.read_csv(full_path, skiprows = [0, 2, 3]))
#        df_list.append(full_path)
#    a = df_list[1]
#    print a
#    df = pd.read_csv(a, skiprows = [0, 2, 3])
    df = df_list[0]
    df.index = pd.to_datetime(df.TIMESTAMP)
    df['modulo_15'] = df.index.second % 15
    
    # Generate the required date ranges for indexing and outputting data
    start = dt.datetime(df.index[0].year, df.index[0].month, df.index[0].day,
                        0, 0, 0)
    end = dt.datetime(df.index[0].year, df.index[0].month, df.index[0].day,
                      23, 58)
    dt_range_1 = pd.date_range(start, end, freq = '2T')
    dt_range_2 = dt_range_1 + dt.timedelta(minutes = 1, seconds = 59, 
                                           microseconds = 500000)
    dt_range_out = dt_range_1 + dt.timedelta(minutes = 2)
    
    # Make a CO2 names list and a reference dictionary for assigning the data to 
    # the output dataframe
    CO2_names_list = ['CO2_{0}m'.format(heights_dict[i]) for i in heights_dict.keys()]    
    CO2_names_dict = dict(zip(heights_dict.keys(), CO2_names_list))
    
    # Make a T names list and a reference dictionary for assigning the data to 
    # the output dataframe
    T_names_list = ['Tair_{0}m'.format(heights_dict[i]) for i in heights_dict.keys()] 
    T_names_dict = dict(zip(heights_dict.keys(), T_names_list))
    
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
    
        # Do the averaging for each valve and send CO2 and temp data to output df
        for j in profile_n:
            rslt_df.loc[dt_range_out[i], CO2_names_dict[j]] = sub_df.loc[j, 'CO2']
            T_name = 'T_air({0})'.format(str(j))
            rslt_df.loc[dt_range_out[i], T_names_dict[j]] = sub_df.loc[j, T_name]
    
    return rslt_df
    
def warra_average():
    
    # Create a dict to reference heights to valve numbers
    profile_n = [1, 2, 3, 4, 5, 6, 7, 8]
    profile_heights = [2, 4, 8, 16, 30, 42, 54, 70]
    heights_dict = heights_dict = dict(zip(profile_n, 
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
#
##------------------------------------------------------------------------------
## Start program
#
#profile_n = [1, 2, 3, 4, 5, 6, 7, 8]
#profile_heights = [2, 4, 8, 16, 30, 42, 54, 70]
#
#path_to_avg_file = '/home/ian/Downloads/TOA5_25569.SiteAvg.dat'
#path_to_raw_file = '/home/ian/Downloads/TOA5_RawData392.dat'
#
#heights_dict = dict(zip(profile_n, [str(height) for height in profile_heights]))
#        
##avg_df = pd.read_csv(path_to_avg_file, skiprows = [0, 2, 3])
##rslt_avg_df = process_avg_series(avg_df, heights_dict)
#
#raw_df = pd.read_csv(path_to_raw_file, skiprows = [0, 2, 3])
#rslt_raw_df = process_raw_series(raw_df, heights_dict)
#
#rslt_avg_df.to_csv('/home/ian/Documents/profile.csv', index_label = 'Datetime')