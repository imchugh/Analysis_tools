# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 10:43:02 2016

@author: imchugh
"""

import os
import DataIO as io
import calendar
import numpy as np

input_path = '/media/Data/Dropbox/Data/Site data processing/HowardSprings/Advanced_v12a'
input_file = 'Advanced_processed_data_HowardSprings_v12a.df'
f = os.path.join(input_path, input_file)

var_list = ['Fsd_Con', 'Fc', 'ustar']
df = io.DINGO_df_to_data_structure(f, var_list = var_list, 
                                   output_structure = 'pandas')

years_list = list(set(df.index.year))
years_list.sort()
ustar_dict = {2001: 0.25,
              2002: 0.25,
              2003: 0.25,
              2004: 0.2,
              2005: 0.23,
              2006: 0.25,
              2007: 0.23,
              2008: 0.25,
              2009: 0.23,
              2010: 0.26,
              2011: 0.27,
              2012: 0.24,
              2013: 0.26,
              2014: 0.29,
              2015: 0.25,
              2016: 0.25}

for year in years_list:

    if calendar.isleap(year):
        recs = 366 * 48 / 2

    this_df = df.loc[str(year)]
    this_df = this_df[this_df['Fsd_Con'] < 5]
    Fc_recs = len(this_df.Fc.dropna())
    filt_Fc_recs = len(this_df.Fc[this_df.ustar > ustar_dict[year]].dropna())
    
    print 'For ' + str(year) + ':'
    print ' - Total available nocturnal records: ' + str(len(this_df))
    print ' - Total possible nocturnal records: ' + str(recs)
    print ' - Total available nocturnal Fc records: ' + str(Fc_recs)
    print ' - Total available u*-filtered nocturnal Fc records: ' + str(filt_Fc_recs)
    print ' - u*-filtered Fc percentage of all possible records: ' + str(np.round(filt_Fc_recs / float(recs) * 100, 1))