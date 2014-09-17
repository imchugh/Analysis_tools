# -*- coding: utf-8 -*-
"""
Created on Mon Sep 15 09:29:43 2014

@author: imchugh
"""

import Tkinter, tkFileDialog
import netCDF4
import pandas as pd
import datetime as dt
import xlrd
import pdb

def netCDF_input():
    
    # Prompt user for configuration file and get it
    root = Tkinter.Tk(); root.withdraw()
    file_in = tkFileDialog.askopenfilename(initialdir='')
    root.destroy()
    
    nc_obj=netCDF4.Dataset(file_in)
    
    dates_list=[dt.datetime(*xlrd.xldate_as_tuple(elem,0)) for elem in nc_obj.variables['xlDateTime']]    
    
    d_data={}
    d_attr=nc_obj.__dict__
    for i in nc_obj.variables.keys():
        ndims=len(nc_obj.variables[i].shape)
        if ndims==3:
            d_data[i]=nc_obj.variables[i][:,0,0]
        elif ndims==1:    
            d_data[i]=nc_obj.variables[i][:]
    nc_obj.close()
    
    df=pd.DataFrame(d_data,index=dates_list)    
    
    return df, d_attr