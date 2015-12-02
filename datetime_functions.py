# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 14:21:38 2015

@author: imchugh
"""

import datetime as dt
import numpy as np
import pdb

def get_timestep(datetime_array):
    """
    Checks timedelta for consistency and finds measurement interval
    Takes datetime array as argument
    Returns measurement interval in hours (float)
    """
    check_timedelta = datetime_array[1: ] - datetime_array[: -1]
    if not all(check_timedelta[0] == rest for rest in check_timedelta):
        print 'Time series is not continuous'
        return None
    else:
        return check_timedelta[0].seconds / 3600.0

def get_moving_window_indices(datetime_array, window, step, retro_stamp = True):
    """
    Finds available windows given window and step
    Takes datetime array, window (int) and step (int) as arguments
    Returns dictionary containing date (centre of window) as key and array
    location indices as value
    Note: retro_stamp indicates that the timestamp applies to the data 
          retrospectively i.e. the 00:00 timestamp indicates a data collection
          period of 23:30-00:00; so for example the 2012-01-01 data day is 
          correctly represented by 2012-01-01 00:30 to 2012-02-01 00:00; if 
          this is set to false, it will interpret the timestamps literally, so 
          the 2012-01-01 data day will be represented by 2012-01-01 00:00 to 
          2012-12-31 23:30. 
          This is only correct if timestamps are not retrospective!
    """

    # Check measurement interval    
    meas_int = get_timestep(datetime_array)
    
    if retro_stamp:
        shift_by = meas_int
    else:
        shift_by = 0

    # Find part days at beginning
    start_date = (dt.datetime.combine(datetime_array[0].date(), 
                                      dt.datetime.min.time()) +
                  dt.timedelta(hours = shift_by))
    if start_date > datetime_array[0]:
        start_index = np.where(datetime_array == start_date)[0].item()
    elif start_date < datetime_array[0]:
        start_date = start_date + dt.timedelta(1)
        start_index = np.where(datetime_array == start_date)[0].item()
    else:
        start_index = 0              

    # Find part days at end
    end_date = (dt.datetime.combine(datetime_array[-1].date(), 
                                    dt.datetime.min.time()) +
                                    dt.timedelta(1) - 
                                    dt.timedelta(hours = meas_int) +
                                    dt.timedelta(hours = shift_by))
    if end_date > datetime_array[-1]:
        end_date = end_date - dt.timedelta(1)
        end_index = np.where(datetime_array == end_date)[0].item()               
    else:
        end_index = len(datetime_array)

    # Slice a new working array
    work_array = datetime_array[start_index: end_index]

    # Generate dates representing begin, centre and end of window
    num_days = (work_array[-1].date() - 
                work_array[0].date()).days + 1 - window
    num_days_range = range(0, num_days + 1, step)
    centre_datetime_array = np.array([(work_array[0] + 
                                       dt.timedelta(i + window / 2.0))
                                      for i in num_days_range])
    begin_datetime_array = np.array(centre_datetime_array - 
                                    dt.timedelta(window / 2.0))
    end_datetime_array = np.array(centre_datetime_array + 
                                  dt.timedelta(window / 2.0) - 
                                  dt.timedelta(hours = meas_int))

    # Create dictionary with date as key and indices as values
    step_dates_index_dict = {}
    centre_date_array = np.array([date_time.date() 
                                  for date_time in centre_datetime_array])
    for i, date in enumerate(centre_date_array):
        begin_ind = np.where(datetime_array == begin_datetime_array[i])[0].item()
        end_ind = np.where(datetime_array == end_datetime_array[i])[0].item()
        step_dates_index_dict[date] = [begin_ind, end_ind]

    return step_dates_index_dict

def get_year_indices(datetime_array, retro_stamp = True):    
    """
    Finds the array location indices for the years;
    Takes datetime array as arg
    Returns dictionary containing year as key and indices as value
    Note: retro_stamp indicates that the timestamp applies to the data 
          retrospectively i.e. the 00:00 timestamp indicates a data collection
          period of 23:30-00:00; so for example the 2012 data year is correctly
          represented by 2012-01-01 00:30 to 2013-01-01 00:00; if this is set 
          to false, it will interpret the timestamps literally, so the 2012 
          data year will be represented by 2012-01-01 00:00 to 2012-12-31 23:30. 
          This is only correct if timestamps are not retrospective!
    """    
    # Check measurement interval    
    meas_int = get_timestep(datetime_array)
    
    if retro_stamp:
        shift_by = meas_int
    else:
        shift_by = 0

    datetime_array = datetime_array - dt.timedelta(hours = shift_by)
    years_index_dict = {}
    year_array = np.array([i.year for i in datetime_array])
    year_list = list(set(year_array))
    for yr in year_list:
        index = np.where(year_array == yr)[0]
        years_index_dict[yr] = [index[0], index[-1]]
        
    return(years_index_dict)

def get_time_indices(datetime_array, time_interval, retro_stamp = True):    
    """
    Finds the array location indices for the years;
    Takes datetime array as arg
    Returns dictionary containing year as key and indices as value
    Note: retro_stamp indicates that the timestamp applies to the data 
          retrospectively i.e. the 00:00 timestamp indicates a data collection
          period of 23:30-00:00; so for example the 2012 data year is correctly
          represented by 2012-01-01 00:30 to 2013-01-01 00:00; if this is set 
          to false, it will interpret the timestamps literally, so the 2012 
          data year will be represented by 2012-01-01 00:00 to 2012-12-31 23:30. 
          This is only correct if timestamps are not retrospective!
    """    
    # Check measurement interval    
    meas_int = get_timestep(datetime_array)
    
    if retro_stamp:
        shift_by = meas_int
    else:
        shift_by = 0



    datetime_array = datetime_array - dt.timedelta(hours = shift_by)
    interval_index_dict = {}
    interval_array = np.array([eval('i.' + time_interval) for i in datetime_array])
    interval_list = list(set(interval_array))
    for intvl in interval_list:
        index = np.where(interval_array == intvl)[0]
        interval_index_dict[intvl] = [index[0], index[-1]]
        
    return(interval_index_dict)

def get_day_indices(datetime_array, time_interval, retro_stamp = True):
    
    # Check measurement interval    
    meas_int = get_timestep(datetime_array)    

    if retro_stamp:
        shift_by = meas_int
    else:
        shift_by = 0

    datetime_array = datetime_array - dt.timedelta(hours = shift_by)
    
    year_list = np.arange(datetime_array[0].year, datetime_array[-1].year)
    month_list = np.arange(1, 13)
    
    date_list = []
    if time_interval == 'year':
        for year in year_list:
            date_list.append(dt.date(year, 1, 1))
    elif time_interval == 'month':
        for year in year_list:
            for month in month_list
    
    # Create a series of continuous whole day dates that will be used for output
    t_delta = (datetime_array[-1].date() - datetime_array[0].date())
    num_days = t_delta.days + 1
    num_months = t_delta.months + 1
    num_years = t_delta.years + 1
    all_dates_array = np.array([datetime_array[0].date() + dt.timedelta(day) 
                                for day in xrange(num_days)])
            
    return 