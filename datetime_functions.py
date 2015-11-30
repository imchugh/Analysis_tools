# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 14:21:38 2015

@author: imchugh
"""

import pandas as pd
import datetime as dt
import numpy as np

window = 10
step = 5
meas_int_hours = 0.5

test = np.array([dt.datetime(2012, 1, 1) + dt.timedelta(minutes = i * 30) 
                 for i in range(52608)])

# Check no time jumps and get measurement interval
check_timedelta = test[1: ] - test[: -1]
if not all(check_timedelta[0] == rest for rest in check_timedelta):
    print 'Time series is not continuous'
else:
    meas_int = check_timedelta[0].seconds / 3600.0

# Create a series of continuous whole day dates that will be used for output
# (parameter series will be interpolated between window centres)
num_days = (test[-1].date() - test[0].date()).days + 1
all_dates_array = np.array([test[0].date() + dt.timedelta(day) 
                            for day in xrange(num_days)])

# Create a shifted array
shift_mins = 60 * meas_int
shift_datetime_array = test - dt.timedelta(minutes = shift_mins)

# Trim to remove part days
start_date = shift_datetime_array[0].date()
num_first_day = len([i for i in shift_datetime_array if i.date() == start_date]
if num_first_day < 24: 
    start_date = start_date + dt.timedelta(1)
end_date = shift_datetime_array[-1].date()
num_last_day = len([i for i in shift_datetime_array if i.date() == end_date]
if num_first_day < 24: start_date = start_date + dt.timedelta(1)
    
    
else:
    start_date = shift_datetime_array[0].date()

if len(shift_datetime_array[shift_datetime_array.date() == 
                            shift_datetime_array[0].date()]) < 24:
    start_date = dt.datetime(shift_datetime_array[0].year,
                             shift_datetime_array[0].month,
                             shift_datetime_array[0].day) + dt
    print start_date                             
else:
    start_date = shift_datetime_array[0].date()
print start_date





#window_width = window * 24 * 1 / meas_int
#num_windows = int((len(all_dates_array) - window_width) / step)
#print num_windows
#start_minimum = shift_datetime_array[0] + dt.timedelta(window / 2.0)
#if not 
#if start_minimum.hour < 12 and start_minimum.hour > 0:
#    hour_split = 12
#else:
#    hour_split = 0

    


    

#for i in all_dates_array:
#    print len(np.where(test.date == i))

# Check that first and last days are complete and revise start and end dates if required
#temp_date = dt.datetime.combine((test[0] + dt.timedelta(1)).date(), 
#                                dt.datetime.min.time())
#num_obs = len(np.where(test < temp_date)[0])
#if num_obs < 24 * (1 / meas_int_hours):
#    start_date = start_date + dt.timedelta(1)
#temp_date = dt.datetime.combine(shift_datetime_array[-1].date(), 
#                                dt.datetime.min.time())
#num_obs = len(np.where(shift_datetime_array >= temp_date)[0])
#if num_obs < 24 * (1 / configs_dict['measurement_interval']):
#    end_date = end_date - dt.timedelta(1)