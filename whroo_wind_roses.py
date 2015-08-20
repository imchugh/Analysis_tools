# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 09:53:14 2015

@author: imchugh
"""

import pandas as pd
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import pdb

import DataIO as io

def get_data():
    
    reload(io)    
    
    file_in = '/home/imchugh/Ozflux/Sites/Whroo/Data/Processed/all/Whroo_2011_to_2014_L6.nc'
    
    return io.OzFluxQCnc_to_pandasDF(file_in)

def wind_roses():
    
    # Program error corrections
    sonic_az = 282
    program_az = 214.5
    correction = sonic_az - program_az
    
    # Get data
    df, attr = get_data()

    # Create variables lists
    wind_spd_list = ['Ws_CSAT',
                     'Ws_RMY_1m',
                     'Ws_RMY_2m', 
                     'Ws_RMY_4m', 
                     'Ws_RMY_8m',
                     'Ws_RMY_16m', 
                     'Ws_RMY_32m']
    wind_dir_list = ['Wd_CSAT',
                     'Wd_RMY_1m',
                     'Wd_RMY_2m', 
                     'Wd_RMY_4m', 
                     'Wd_RMY_8m',
                     'Wd_RMY_16m', 
                     'Wd_RMY_32m']

    # Remove bad data
    for var in wind_spd_list:
        df[var] = np.where(df[var + '_QCFlag'] == 0, df[var], np.nan)
        df.drop(var + '_QCFlag', axis = 1, inplace = True)
    for var in wind_dir_list:
        df[var] = np.where(df[var + '_QCFlag']==0, df[var], np.nan)    
        df.drop(var + '_QCFlag', axis = 1, inplace = True)

    # Subset df columns    
    df = df[wind_spd_list + wind_dir_list + ['ustar', 'Fsd']]

    # Remove 1) daytime, 2) zero wind speeds, 3) nans
    df = df[df.Fsd < 10]
    df.drop('Fsd', axis = 1, inplace = True)
    for i in range(len(wind_spd_list)):
        df[wind_spd_list[i]][df[wind_spd_list[i]] == 0] = np.nan
    df.dropna(axis = 0, inplace = True)

    # Correct CSAT wind direction for incorrect angle in program
    df['Wd_CSAT'] = df['Wd_CSAT'] + correction
    df['Wd_CSAT'] = np.where(df['Wd_CSAT'] > 360, df['Wd_CSAT'] - 360, 
                             df['Wd_CSAT'])

    # Correct 1m windspeed for error
    # Correct all others for magnetic

    # Create a dictionary of quantiles that will define ustar classes
    ustar_dict = {}
    ustar_dict['lower'] = [0, 0.6]
    ustar_dict['upper'] = [0.2, 1]

    # Create a directions dictionary
    directions_dict = {}
    directions_dict['lower'] = [sector * 22.5 - 11.25 for sector in range(0, 16)]
    directions_dict['upper'] = [i + 22.5 for i in directions_dict['lower']]
    directions_dict['lower'][0] = directions_dict['lower'][0] + 360

    # Create x polar coordinates against which to map wind direction frequencies 
    theta_1 = np.array([348.75])
    theta_2 = np.linspace(11.25, 326.25, 15)
    theta = np.radians(np.append(theta_1, theta_2))
    width = np.pi/ 8    

    # Create axes
    fig, ax = plt.subplots(5, 2, figsize = (8, 20), subplot_kw = dict(polar = True))
    fig.patch.set_facecolor('white')

    # Iterate through ustar quantile boundaries
    for i in range(2):

        # Retrieve quantile boundaries and subset the df
        lower_quantile = df['ustar'].quantile(ustar_dict['lower'][i])
        upper_quantile = df['ustar'].quantile(ustar_dict['upper'][i])
        test_df = df[(df.ustar >= lower_quantile) & (df.ustar < upper_quantile)]

        # Calculate frequencies for wind directions
        rslt_df = pd.DataFrame(columns = wind_dir_list[2:], index = range(0, 16))
        for var in rslt_df.columns:
            for sector in xrange(0, 16):
                lower_bound = directions_dict['lower'][sector]
                upper_bound = directions_dict['upper'][sector]
                if sector == 0:
                    rslt_df[var].iloc[sector] = len(test_df[var][(test_df[var] >= lower_bound) | 
                                                   (test_df[var] < upper_bound)])
                else:
                    rslt_df[var].iloc[sector] = len(test_df[var][(test_df[var] >= lower_bound) & 
                                                   (test_df[var] < upper_bound)])

        # Convert to percent
        rslt_df = rslt_df / rslt_df.sum() * 100
    
        # Plot
        for j, var in enumerate(rslt_df.columns):
            this_ax = ax[j, i]
            this_ax.set_theta_zero_location('N')
            this_ax.set_theta_direction(-1)
            this_ax.bar(theta, rslt_df[var], width = width, alpha=0.5,
                        facecolor='0.5', edgecolor = 'None')
            plt.setp(this_ax.get_xticklabels(), visible = False)
            [this_ax.axvline(line, color = '0.5') for line in 
             np.radians([0, 90, 180, 270])]
        
    plt.tight_layout()
    
    return df