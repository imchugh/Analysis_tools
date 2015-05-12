# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import os
import numpy as np
from scipy.optimize import curve_fit
import pdb

file_in=os.path.join('/home/imchugh/Temp','GSI.xlsx')

df=pd.read_excel(file_in)

def sma(data,window):
    
    weights = np.repeat(1.0, window)/window
    smas = np.convolve(data, weights, 'valid')
    return smas

def GSI(data_array,Ta_lo,Ta_hi,DL_lo,DL_hi):

    Ta_array = data_array[:,0]
    DL_array = data_array[:,1]
    VPD_array = data_array[:,2]
    snow_array = data_array[:,3]
    
    Ta_ind_array = np.where(Ta_array < Ta_lo, 0, np.where(Ta_array > Ta_hi, 1, (Ta_array - Ta_lo) / (Ta_hi - Ta_lo)))
    
    DL_ind_array = np.where(DL_array < DL_lo, 0, np.where(DL_array > DL_hi, 1, (DL_array - DL_lo) / (DL_hi - DL_lo)))
    
    VPD_ind_array = np.where(VPD_array < 1, 1, np.where(VPD_array > 4.1, 0, 1 - (VPD_array - 1) / (4.1 - 1)))
    
    GSI_array = Ta_ind_array * DL_ind_array * VPD_ind_array * snow_array
    
    GSI_array_ext = np.tile(GSI_array, 3)
    
    GSI_array_ext_run = sma(GSI_array_ext, 21)
    
    GSI_array_run = GSI_array_ext_run[721:1452]
    
    return GSI_array_run

in_array = np.array(df[['Ta','Daylength','VPD','Snow']][df.Aopt_norm!=np.nan])
out_array = np.array(df['Aopt_norm'][df.Aopt_norm!=np.nan])

test_rslt, test_cov = curve_fit(GSI, in_array, out_array, p0=[-2,5,8,10])

