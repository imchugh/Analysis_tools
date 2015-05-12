# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 02:57:37 2015

Routines for cleaning up profile data

@author: imchugh
"""

import pandas as pd
import os
import pdb
import numpy as np
from scipy import stats
import datetime as dt

def correctTa_data():
    
    """
    Requires temperature and pressure data as well as temperature and pressure 
    reference variables in the input data files; 
    these reference variables are used to fill data gaps
    """
    
    # Set IO
    path = '/home/imchugh/Processing/Whroo/Profile data'
    inname = 'profile_met_uncorrected.dat'
    outname = 'profile_met_corrected.dat'
    inpathname=os.path.join(path, inname)
    outpathname=os.path.join(path, outname)
    
    # Set var names
    # Ta
    Ta_names = ['Ta_HMP_1m_Avg', 'Ta_HMP_2m_Avg', 'Ta_HMP_4m_Avg', 'Ta_HMP_8m_Avg', 
                'Ta_HMP_16m_Avg', 'Ta_HMP_32m_Avg']
    Ta_reference = 'Ta_HMP_36m'
    # P            
    ps_names = ['ps']
    ps_reference = 'ps_EC'
    
    # Set bad data descriptors
    bad_data = {}
    bad_data['Ta_HMP_1m_Avg'] = ['2012-10-20', '2013-04-10']
    bad_data['Ta_HMP_8m_Avg'] = ['2013-07-18', '2014-05-08']
    bad_data['Ta_HMP_32m_Avg'] = [['2011-12-02','2012-06-28'],
 				          ['2013-02-01','2013-04-10']]
    bad_all_data = ['2014-06-11', None]
    
    # Set other stuff
    true_heights=[0.5,2,4,8,16,32]
    Ta_range_limits=[-20,50]
    p_range_limits=[90,110]

    # Import data
    df = pd.read_csv(inpathname, parse_dates=[0], index_col=0)

    # Remove known bad data
    for i in bad_data.keys():
        if not type(bad_data[i][0]) is list:
            df[i].loc[bad_data[i][0]: bad_data[i][1]] = np.nan
        else:
            for j in bad_data[i]:
                df[i].loc[j[0]: j[1]] = np.nan

    # remove data outside range limits or occurring during universal bad data period
    # for temps
    for i in Ta_names:
        df[i]=np.where(df[i] < Ta_range_limits[0], np.nan,
                       np.where(df[i] > Ta_range_limits[1] , np.nan, df[i]))
        df[i].loc[bad_all_data[0]: bad_all_data[1]] = np.nan
    # For pressure
    for i in ps_names:
        df[i] = df[i] / 10

        df[i]=np.where(df[i] < p_range_limits[0], np.nan, 
                       np.where(df[i] > p_range_limits[1], np.nan, df[i]))
        df[i].loc[bad_all_data[0]: bad_all_data[1]] = np.nan
  			
    # Rename the columns using the true heights
    new_Ta_names = ['Ta_'+str(i)+'m' for i in true_heights]
    new_labels = {Ta_names[i]: new_Ta_names[i] for i in range(6)}
    df = df.rename(columns = new_labels)

    # Fill temperature data for each level from level with best correlation
    for i in new_Ta_names:
        
        vars_list = list(new_Ta_names)
        vars_list.remove(i)
        results_df = pd.DataFrame(columns=['slope','int','r2','p_val','se'],index=vars_list)
        for j in vars_list:
            temp_df=df[[i,j]].dropna(how = 'any', axis = 0)
            results_df.ix[j] = stats.linregress(temp_df[j], temp_df[i])
        results_df.sort('r2', ascending=False, inplace=True)
        vars_list = results_df.index
        for j in vars_list:
            if len(df[j].dropna())==len(df):
                break
            temp_S=df[j] * results_df['slope'].loc[j] + results_df['int'].loc[j]
            df[i] = np.where(pd.isnull(df[i]), temp_S, df[i])
	
    # Fill and remaining missing temperature data using the reference temperature
    for i in new_Ta_names:
        
        temp_df = df[[i, Ta_reference]].dropna(how = 'any', axis = 0)
        result = stats.linregress(temp_df[Ta_reference], temp_df[i])
        temp_S = df[Ta_reference] * result[0] + result[1]
        df[i] = np.where(pd.isnull(df[i]), temp_S, df[i])

    # Fill missing pressure data using the reference pressure
    for i in ps_names:
        
        temp_df = df[[i, ps_reference]].dropna(how = 'any', axis = 0)
        result = stats.linregress(temp_df[ps_reference], temp_df[i])
        temp_S = df[ps_reference] * result[0] + result[1]
        df[i] = np.where(pd.isnull(df[i]), temp_S, df[i])
        
    # Fix the obviously wrong temperature near ground level (due to bug in CSI AM25T instruction?)
    df[new_Ta_names[0]] = df[new_Ta_names[0]] + (df[new_Ta_names[1]].mean() - df[new_Ta_names[0]].mean())
    
    # Output data
    df.to_csv(outpathname, index_label = 'TIMESTAMP')
    
    return df
    
def correctCO2_data():
    
    # Set IO
    path='/home/imchugh/Processing/Whroo/Profile data/'
    inname='profile_CO2_uncorrected.dat'
    outname='profile_CO2_corrected.dat'
    inpathname=os.path.join(path, inname)
    outpathname=os.path.join(path, outname)
    
    # Set var names
    CO2_names = ['Cc_LI840_1m', 'Cc_LI840_2m', 'Cc_LI840_4m', 'Cc_LI840_8m', 
                'Cc_LI840_16m', 'Cc_LI840_32m']
    
    # Set bad data descriptors
    baddata_dates = ['2013-08-24','2013-10-29']
    badcoeff_dates = ['2012-06-28 11:00:00','2012-10-17 12:50:00']
    instrument_dates = [['2011-12-02','2012-06-28'],
                        ['2012-06-28','2012-10-13'],
                        ['2012-10-13','2013-08-23'],
        		      ['2013-10-29','2014-06-02']]
    instrument_times = [['12:00:00','10:58:00'],
                        ['11:00:00','12:00:00'],
    			      ['12:02:00','23:58:00'],
    			      ['12:00:00','23:58:00']]
    
    # Set other stuff    
    coeff_correct=2.5
    true_heights=[0.5,2,4,8,16,32]
    CO2_range=[300,600]
    CO2_base=390
    		
    # Import data
    df = pd.read_csv(inpathname, parse_dates = [0], index_col = 0)
    
    # Get the CO2 columns
    CO2_cols_list = [i for i in df.columns if 'Cc' in i]
    
    # Fix time stamps    
    
    # Correct the data for 1) no data; 2) wrong instrument scaling coefficients; 
    # 3) range checks; 4) reversed label assignment of CO2
	
    # 1 and 2 above
    for i in CO2_names:
        df.loc[baddata_dates[0]:baddata_dates[1]] = np.nan
        df[i].loc[badcoeff_dates[0]: badcoeff_dates[1]] = df[i].loc[badcoeff_dates[0]: badcoeff_dates[1]] * coeff_correct
    		
    # 3 above
    for i in CO2_names:
        df[i]=np.where(df[i]<CO2_range[0],np.nan,np.where(df[i]>CO2_range[1],np.nan,df[i]))
    
    # 4 above (reverse the column names and reassign the heights)
    	
    # Rename the columns using the true heights
    true_heights.reverse()
    new_CO2_names=['Cc_'+str(i)+'m' for i in true_heights]
    reverse_dict = {CO2_names[i]: new_CO2_names[i] for i in range(len(CO2_names))}
    df = df.rename(columns = reverse_dict)
    
    # Reorder the columns so that they are ascending
    df = df[df.columns[::-1]]
           
    # Output the data	
    df.to_csv(outpathname, index_label = 'TIMESTAMP')
    return df       