# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 00:22:00 2015

@author: imchugh
"""
import os
import pandas as pd

path='/media/Data/Dropbox/Data_sites non flux/Site data plots and non-flux/Sites/Whroo/Data/Species composition'
name='Whroo understorey species Apr 2013.xls'

filename=os.path.join(path,name)

# Get list of species onsite
species_df=pd.read_excel(filename,sheetname='Species ID')
#species_dict={species_df['Species_code'].iloc[k]: species_df['Species'].iloc[k] for k in range(len(species_df['Species_code']))}
species_df.index=species_df['Species_code']
species_df.drop('Species_code', axis=1, inplace=True)
species_df['Count']=0

# Get list of ground cover classes used
cover_df=pd.read_excel(filename,sheetname='Ground classes ID')
#cover_dict={cover_df['Symbol'].iloc[k]: cover_df['Description'].iloc[k] for k in range(len(cover_df['Symbol']))}
cover_df.index=cover_df['Symbol']
cover_df.drop('Symbol', axis=1, inplace=True)
cover_df['Count']=0

# Get data
data_df=pd.read_excel(filename,sheetname='Observations',skiprows=[0,1])
cover_list=[i for i in data_df.columns if 'Ground class' in i]
cover_series=pd.concat([data_df[i] for i in cover_list])
species_list=[i for i in data_df.columns if 'Species (uppermost intercept height - cm)' in i]
species_series=pd.concat([data_df[i] for i in species_list]).reset_index(drop=True)
data_df=pd.DataFrame({'Species and height':,'':})
species_trunc=species_series.dropna().reset_index()
species_trunc.columns=['index','species']

# Find instances of multiple species intercepts
multiples_list=[obs for obs in range(len(species_trunc)) if len(species_trunc['species'].iloc[obs].split(',')) >1]
for ind in multiples_list

#split_list=[obs.split('(') for obs in species_trunc]
#code_list=[obs[0] for obs in split_list]
#heights_list=[obs[1].strip(')') for obs in split_list]