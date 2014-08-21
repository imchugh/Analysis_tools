# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 16:47:37 2014

@author: imchugh
"""

import pandas as pd
import os

import plots_generic as pg

reload(pg)

def calc_cml_sums(path,name,C_name):
    
    convert_const_1=12
    convert_const_2=0.0864
        
    df=pd.read_pickle(os.path.join(path,name))
    
    daily_df=df[C_name].groupby([lambda x: x.year,lambda y: y.dayofyear]).mean()*12*0.0864
    yr_comp_df=pd.DataFrame(index=daily_df.index.levels[1],columns=['2012','2013'])
    
    yr_comp_df['2012']=daily_df.ix[2012]
    yr_comp_df['2013']=daily_df.ix[2013]
    
    yr_comp_df['2012_cml']=yr_comp_df['2012'].cumsum()
    yr_comp_df['2013_cml']=yr_comp_df['2013'].cumsum()
    
    pg.line_plot(yr_comp_df[['2012_cml','2013_cml']],'Test \n')    
    
    return yr_comp_df
    
    



