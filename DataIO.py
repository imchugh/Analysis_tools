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

def OzFluxQCnc_to_pandasDF():
    
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

# Multiple sheets are returned as dictionary (pandas dataframe) objects
def xlsx_to_pandas(header=True,header_row=0,skiprows_after_header=0,date_col=True):
    
    # Prompt user for configuration file and get it
    root = Tkinter.Tk(); root.withdraw()
    file_in = tkFileDialog.askopenfilename(initialdir='')
    root.destroy()

    xl_book=xlrd.open_workbook(file_in)
    
    d={}
    
    start_date='1900-01-01'
    end_date='2100-01-01'    
    
    for sheet_name in xl_book.sheet_names():

        sheet=xl_book.sheet_by_name(sheet_name)

        rows=sheet.nrows
        cols=sheet.ncols
        
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
                    df[column_names[i]]=[sheet.cell_value(j,i) for j in xrange(header_row+skiprows_after_header+1,sheet.nrows)]
                df_ind=pd.date_range(start=df.index[0],end=df.index[-1],freq='30Min')
                df=df.reindex(df_ind)
                d[sheet_name]=df
                
    return d

        