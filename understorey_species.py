# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 00:22:00 2015

@author: imchugh
"""
import numpy as np
import os
import pandas as pd
import pdb


path='/media/Data/Dropbox/Data_sites non flux/Site data plots and non-flux/Sites/Whroo/Data/Species composition'
name='Whroo understorey species Apr 2013.xls'

filename=os.path.join(path,name)

# Get list of species onsite
speciesID_df=pd.read_excel(filename,sheetname='Species ID')
speciesID_df.index=speciesID_df['Species_code']
speciesID_df.drop('Species_code', axis=1, inplace=True)
speciesID_df['Count']=0
speciesID_df['Mean_height']=0

# Get list of ground cover classes used
coverID_df=pd.read_excel(filename,sheetname='Ground classes ID')
coverID_df.index=coverID_df['Symbol']
coverID_df.drop('Symbol', axis=1, inplace=True)
coverID_df['Count']=0

# Get data
data_df=pd.read_excel(filename,sheetname='Observations',skiprows=[0,1])
cover_headers_list=[i for i in data_df.columns if 'Ground class' in i]
cover_series=pd.concat([data_df[i] for i in cover_headers_list]).reset_index(drop=True)
species_headers_list=[i for i in data_df.columns if 'Species' in i]
heights_headers_list=[i for i in data_df.columns if 'intercept height' in i]
species_series=pd.concat([data_df[i] for i in species_headers_list]).reset_index(drop=True)
heights_series=pd.concat([data_df[i] for i in heights_headers_list]).reset_index(drop=True)
species_heights_df=pd.concat([species_series,heights_series],axis=1)
species_heights_df.columns=['species','heights']
trunc_df=species_heights_df.dropna().reset_index()
trunc_df.columns=['index','species','heights']

# Find instances of multiple species intercepts
multiples_list=[obs for obs in range(len(trunc_df)) if len(trunc_df['species'].iloc[obs].split(',')) >1]

# Rewrite all multiple observations as additional observations
for i in xrange(len(multiples_list)):
    multiple_species_list=trunc_df['species'].iloc[multiples_list[i]].split(',')
    multiple_species_list=[species.strip() for species in multiple_species_list]
    multiple_heights_list=trunc_df['heights'].iloc[multiples_list[i]].split(',')
    temp_df=pd.DataFrame({'index':np.ones(len(multiple_species_list))*trunc_df['index'].iloc[multiples_list[i]],
                          'species':multiple_species_list,
                          'heights':multiple_heights_list})
    trunc_df.drop(multiples_list[i],inplace=True)
    trunc_df=trunc_df.append(temp_df)
trunc_df.sort('index',inplace=True)
trunc_df.reset_index(drop=True,inplace=True)

# Do stats for species composition and abundance
for code in speciesID_df.index:
    speciesID_df['Count'].loc[code]=len(trunc_df[trunc_df.species==code])
    speciesID_df['Mean_height'].loc[code]=trunc_df['heights'][trunc_df.species==code].astype(float).mean()
    speciesID_df['Mean_height']=speciesID_df['Mean_height'].round(1)
speciesID_df.sort('Count',inplace=True,ascending=False)
speciesID_df['Percent']=(speciesID_df['Count']/float(len(trunc_df))*100).round(1)

# Do stats for ground cover
for code in coverID_df.index:
    coverID_df['Count'].loc[code]=len(cover_series[cover_series==code])
coverID_df.sort('Count',inplace=True,ascending=False)
coverID_df['Percent']=(coverID_df['Count']/float(len(cover_series))*100).round(1)

# Ouput data
filename=os.path.join(path,'species.csv')
speciesID_df.to_csv(filename,index_label='species_code')
filename=os.path.join(path,'cover.csv')
coverID_df.to_csv(filename,index_label='cover_code')

temp_list=[]
for code in coverID_df.index:
    temp_list=temp_list+list(cover_series[cover_series==code].index)