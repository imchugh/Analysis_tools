# -*- coding: utf-8 -*-
"""
Created on Mon Sep 15 09:29:43 2014

@author: imchugh
"""

import Tkinter, tkFileDialog
from configobj import ConfigObj
import netCDF4
import numpy as np
import pandas as pd
import datetime as dt
import xlrd
import pdb

def file_select_dialog():
    """ Open a file select dialog to get path for file retrieval"""
    
    root = Tkinter.Tk(); root.withdraw()
    file_in = tkFileDialog.askopenfilename(initialdir='')
    root.destroy()
    return file_in

def read_config_file(file_in):
    
    return ConfigObj(file_in)

def OzFluxQCnc_to_pandasDF(file_in, var_list = None):
    
    nc_obj=netCDF4.Dataset(file_in)

    dates_list=[dt.datetime(*xlrd.xldate_as_tuple(elem,0)) for elem in nc_obj.variables['xlDateTime']]    
    
    d_data={}
    d_attr=nc_obj.__dict__

    if var_list == None: 
        var_list = nc_obj.variables.keys()

    for var in var_list:
        ndims = len(nc_obj.variables[var].shape)
        if ndims == 3:
            d_data[var] = nc_obj.variables[var][:, 0, 0]
        elif ndims == 1:
            d_data[var] = nc_obj.variables[var][:]
    nc_obj.close()
    
    df = pd.DataFrame(d_data, index = dates_list)
    
    df.replace(-9999, np.nan, inplace = True)    

    return df, d_attr

def pandas_to_nc(file_in, file_out):
    
    return

def xlsx_to_pandas(file_in,header=True,header_row=0,skiprows_after_header=0,date_col=True,regularise=True,worksheets=[]):

    xl_book=xlrd.open_workbook(file_in)
    
    d={}
    
    start_date='1900-01-01'
    end_date='2100-01-01'    
    
    if not worksheets:
        get_sheets=xl_book.sheet_names()
    else:
        get_sheets=worksheets
    
    for sheet_name in get_sheets:

        sheet=xl_book.sheet_by_name(sheet_name)
       
        rows=sheet.nrows
        cols=sheet.ncols
        
        if rows==0: print 'Could not find any valid rows'
        if cols==0: print 'Could not find any valid columns'

        if rows!=0 and cols!=0:
            if header==True:
                if date_col==True:
                    column_names=[str(sheet.cell_value(header_row,i)) for i in range(cols)]
                    index=[]
                    for i in xrange(header_row+skiprows_after_header+1,sheet.nrows):
                        try:
                            index.append(dt.datetime(*xlrd.xldate_as_tuple(sheet.cell_value(i,0), xl_book.datemode)))
                        except ValueError:
                            index.append('')
                            print 'Error in sheet '+sheet_name+' at row '+str(i)+'; missing or invalid datetime stamp! Skipping...'
                    df=pd.DataFrame(columns=column_names[1:],index=index)
                    for i in range(1,cols):
                        arr=np.array(sheet.col_values(i)[header_row+skiprows_after_header+1:])
                        arr[arr=='']='-9999'
                        df[column_names[i]]=arr.astype(np.float)
                    if regularise==True:
                        df_freq=pd.infer_freq(df.index)
                        df_ind=pd.date_range(start=df.index[0],end=df.index[-1],freq=df_freq)
                        df=df.reindex(df_ind)
                    d[sheet_name]=df
        else:
            d[sheet_name]=pd.DataFrame()
    return d          
    
    # Multiple sheets are returned as dictionary (pandas dataframe) objects
    # Note XLRD cell type codes:
    #    XL_CELL_EMPTY: 0
    #    XL_CELL_TEXT: 1 (STRING)
    #    XL_CELL_NUMBER: 2 (FLOAT)
    #    XL_CELL_DATE: 3 (FLOAT)
    #    XL_CELL_BOOLEAN: 4 (INT)
    #    XL_CELL_ERROR: 5 (INTERNAL EXCEL CODE)
    #    XL_CELL_BLANK: 6 (EMPTY STRING)