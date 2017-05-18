#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 10:07:35 2017

@author: ian
"""
import numpy as np
import pandas as pd
import datetime as dt
import profile_data_processing as pdp
import os
import pdb

#------------------------------------------------------------------------------
# Common scripts
#------------------------------------------------------------------------------

def get_site_data(site_name):
    
    sites_dict = {'cumberland plains': cumberland_plains, 
                  'howard springs': howard_springs,
                  'warra_avg': warra_average,
                  'warra_raw': warra_raw,
                  'whroo': whroo} 
    
    return sites_dict[site_name]()

#------------------------------------------------------------------------------    
# Site scripts
#------------------------------------------------------------------------------

###############################################################################
# Cumberland Plains                                                           #
###############################################################################

def cumberland_plains():
    
    def interpolate_T(T_series_lower, T_series_upper, T_height):
           
        T1_height = 7
        T2_height = 30
        dtdz = (T_series_upper - T_series_lower) / (T2_height - T1_height)
        
        return (T_height - T1_height) * dtdz + T_series_lower
    
    max_CO2 = 1000
    min_CO2 = 300    
    target_var = 'CO2_Avg'
    missing_data_float = -9999
    input_profile_file_path = '/home/ian/OzFlux/Sites/CumberlandPlains/Data/Profile/EFS_S00_CO2PROFILE_R_2017.csv'
    input_met_data_path = '/home/ian/OzFlux/Sites/CumberlandPlains/Data/Profile/EddyFlux_slow_met_2017.csv'
    
    # Get data and create index 
    df = pd.read_csv(input_profile_file_path)
    df.index = [dt.datetime.combine(dt.datetime.strptime(
                                        df.Date[i], '%d/%m/%Y').date(),
                                    dt.datetime.strptime(
                                        df.Time[i], '%H:%M:%S').time()) 
                for i in xrange(len(df))]
    
    # Remove missing_data_values and replace with NaN
    df[target_var] = np.where(df[target_var] == missing_data_float, 
                              np.nan, 
                              df[target_var])
    total_obs = len(df)
    missing_obs = total_obs - len(df[target_var].dropna())
    print ('{0}% of records contain no data'
           .format(str(round(missing_obs / float(total_obs),1))))
    
    # Remove data that is outside specified range and replace with NaN
    df[target_var] = np.where((df[target_var] < min_CO2)|
                              (df[target_var] > max_CO2),
                              np.nan, 
                              df[target_var])
    dropped_obs = len(df[target_var].dropna()) - missing_obs
    print ('{0}% of records contain data outside range limits'
           .format(str(round(dropped_obs / float(total_obs),1))))
    
    # Remove irregular time stamps then align all to half-hour (+59 seconds!!!)
    startdate_lst = []
    enddate_lst = []
    series_list = []
    valve_list = [1, 2, 3, 4, 5, 6, 7, 8]
    heights_list = [0.5, 1, 2, 3.5, 7, 12, 20, 29]
    for i, valve in enumerate(valve_list):
        new_name = 'CO2_{0}m'.format(str(heights_list[i]))
        expected_minute = valve - 1
        temp_s = df[target_var][df['ValveNo'] == valve].copy()
        index = [i for i in temp_s.index if i.minute %30 != expected_minute]
        if not len(index) == 0:
            print ('{0} instances of inconsistent time step found for valve {1} '
                   '({2}); removed!'
                   .format(str(len(index)), 
                           str(valve), 
                           ','.join([str(i) for i in index])))               
            temp_s.drop(index, axis = 0, inplace = True)
        new_index = [i - dt.timedelta(minutes = expected_minute,
                                      seconds = 59) for i in temp_s.index]    
        temp_s.index = new_index
        temp_s.name = new_name
        startdate_lst.append(temp_s.index[0])
        enddate_lst.append(temp_s.index[-1])
        series_list.append(temp_s)
    
    # Sort the start and end date lists and create a universal date range to be applied to all valve series'
    startdate_lst.sort()
    enddate_lst.sort()    
    new_index = pd.date_range(startdate_lst[0], enddate_lst[-1], freq='30T')
    
    # Reindex, write to dataframe and return
    for series in series_list:
        series = series.reindex(new_index)   
    output_df = pd.concat(series_list, axis = 1)
    
    # Import met data and align time index with profile data
    met_df = pd.read_csv(input_met_data_path)
    met_df.index = [dt.datetime.combine(dt.datetime.strptime(
                                            met_df.Date[i], '%d/%m/%Y').date(),
                                        dt.datetime.strptime(
                                            met_df.Time[i], '%H:%M:%S').time()) 
                    for i in xrange(len(met_df))]
    met_df = met_df[['Ta_HMP_01_Avg', 'Ta_HMP_155_Avg', 'ps_7500_Avg']]
    met_df.columns = ['Tair_lower', 'Tair_upper', 'ps']
    met_df = met_df.reindex(output_df.index)
    
    # Interpolate the temperature data to the same heights as the CO2 
    # measurements
    for this_height in heights_list:
        T_name = 'Tair_{0}m'.format(str(this_height))
        output_df[T_name] = interpolate_T(met_df.Tair_lower, met_df.Tair_upper, 
                                          this_height)
    output_df['ps'] = met_df.ps
    
    return output_df


###############################################################################
# Howard Springs                                                              #
###############################################################################

def howard_springs():

#    dir_str = pdp.dir_select_dialog()        
    dir_str = '/home/ian/OzFlux/Sites/Howard_Springs/Data/Profile/'
    dir_list = os.listdir(dir_str)
    dir_list = [f for f in dir_list if 'Howard_profile_Slow_avg' in f]

    df_list = []
    non_CO2_vars = ['TIMESTAMP', 'T_air_Avg']
    for f in dir_list:
        full_path = os.path.join(dir_str, f)
        this_df = pd.read_csv(full_path, skiprows = [0, 2, 3])
        these_cols = this_df.columns
        var_list = [var for var in these_cols if 'Cc' in var]
        CO2_df = this_df[var_list]
        CO2_df.columns = [var.split('_Avg')[0] for var in var_list]
        CO2_df = CO2_df.join(this_df.loc[:, non_CO2_vars])
        df_list.append(CO2_df)

    df = pd.concat(df_list)
    df.drop_duplicates('TIMESTAMP', inplace = True)
    df.index = pd.to_datetime(df.TIMESTAMP)   
    df.drop('TIMESTAMP', axis = 1, inplace = True)
    old_names = [i for i in df.columns if 'Cc' in i]
    new_names = ['CO2_{0}'.format(this_name.split('_')[2]) 
                 for this_name in old_names]
    new_names.append('Tair_2m')
    df.columns = new_names

    return df

###############################################################################
# Warra
###############################################################################

def warra_raw():

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
    dir_str = '/home/ian/OzFlux/Sites/Warra/Data/Profile/Raw' #pdp.dir_select_dialog()        
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
    
    return rslt_df

#------------------------------------------------------------------------------
        
def warra_average():
    
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
# Warra
###############################################################################

def whroo():
    
    """
    This script should be used to make any adjustments to CO2 data; it drops all
    other variables except the timestamp and the CO2 data because it assumes
    that these variables have been used to make corrections and are 
    redundant when the data is returned (as dictionary)
    """
 
    # Set locations   
    path = '/home/ian/OzFlux/Sites/Whroo/Data/Processed/Profile/'
    file_list = os.listdir(path)
    profile_file_list = [f for f in file_list if 'IRGA' in f]
    met_file_list = [f for f in file_list if not 'IRGA' in f]
    
    # Set var names
    CO2_vars_list = ['Cc_LI840_1m', 'Cc_LI840_2m', 'Cc_LI840_4m', 'Cc_LI840_8m', 
                     'Cc_LI840_16m', 'Cc_LI840_32m']
    
    # Set dates for correction
    last_1min_date = '2012-02-28 12:03:00'
    first_2min_date = '2012-02-28 12:10:00'
    baddata_dates = [['2013-08-24', '2013-10-29'],
                     ['2016-10-12 12:00', 
                      dt.datetime.strftime(dt.datetime.today().date(), 
                                           '%Y-%m-%d')]]
    badcoeff_dates = ['2012-06-28 11:00:00', '2012-10-17 12:50:00']
    instrument_dates = [['2011-12-02 12:00:00', '2012-06-28 10:58:00'],
                        ['2012-06-28 11:00:00', '2012-10-13 12:00:00'],
                        ['2012-10-13 12:02:00', '2013-08-23 23:58:00'],
                        ['2013-10-29 12:00:00', '2014-06-02 23:58:00']]
    
    # Set some other stuff    
    coeff_correct = 2.5
    true_heights = [0.5, 2, 4, 8, 16, 36]
    CO2_range = [300, 600]
    
    # Import datasets, adjust frequency for 1-minute data then concatenate,
    # drop duplicates and extraneous data and ensure no missing time stamps
    profile_df_list = []
    for f in profile_file_list:
        this_df = pd.read_csv(os.path.join(path, f), skiprows = [0,2,3],
                              na_values = 'NAN')
        profile_df_list.append(this_df)
    profile_df = pd.concat(profile_df_list)
    profile_df.drop_duplicates('TIMESTAMP', inplace = True)
    profile_df.index = pd.to_datetime(profile_df.TIMESTAMP)
    profile_df.sort_index(inplace = True)
    profile_df = pd.concat([profile_df.loc[:last_1min_date].resample('2T').mean(),
                            profile_df.loc[first_2min_date:]])
    profile_df = profile_df.reindex(pd.date_range(profile_df.index[0], 
                                                  profile_df.index[-1], 
                                                  freq = '2T'))
    profile_df = profile_df[CO2_vars_list]
    
    # Remove block bad data
    for date_span in baddata_dates:
        profile_df.loc[date_span[0]: date_span[1], CO2_vars_list] = np.nan
    
    # Correct bad program scalar
    profile_df.loc[badcoeff_dates[0]: 
                   badcoeff_dates[1], CO2_vars_list] *= coeff_correct
    
    # Impose range limits
    for this_var in CO2_vars_list:
        profile_df.loc[:, this_var][profile_df.loc[:, this_var] < CO2_range[0]] = np.nan
        profile_df.loc[:, this_var][profile_df.loc[:, this_var] > CO2_range[1]] = np.nan
    
    # Reverse names, which have always referred to wrong heights
    true_heights.reverse()
    new_CO2_vars_list = ['CO2_{0}m'.format(str(i)) for i in true_heights]
    reverse_dict = {CO2_vars_list[i]: new_CO2_vars_list[i] 
                    for i in range(len(CO2_vars_list))}
    profile_df = profile_df.rename(columns = reverse_dict)
    
    # Downsample CO2 to match temperature and pressure data
    trunc_profile_df = pdp.downsample_data(profile_df)
  
    # Open temperature and pressure series and align and join with CO2
    met_df_list = []
    for f in met_file_list:
        this_df = pd.read_csv(os.path.join(path, f), 
                              skiprows = [0,2,3], na_values = 'NAN')
        this_df.index = pd.to_datetime(this_df.TIMESTAMP)
        met_df_list.append(this_df)
    if len(met_df_list) > 1:
        met_df = pd.concat(met_df_list)
    else:
        met_df = met_df_list[0]
    met_df.sort_index(inplace = True)
    met_df.drop_duplicates(inplace = True)
    met_df = met_df.reindex(profile_df.index)
    trunc_profile_df['Tair_36m'] = met_df['Ta_HMP_Avg']
    trunc_profile_df['ps'] = met_df['ps_7500_Avg']
    
    return trunc_profile_df

###############################################################################