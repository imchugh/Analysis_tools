# -*- coding: utf-8 -*-
"""
Created on Wed Sep 17 13:52:53 2014

@author: imchugh
"""
import datetime as dt
import pandas as pd
import pdb

def daily_mean(df,expr):
    if isinstance(df,pd.Series):
        cols=[df.name]
    else:
        cols=df.columns
    daily_df=df.groupby([lambda x: x.year, lambda y: y.dayofyear]).mean()*eval(expr)
    daily_df=daily_df.reset_index()
    daily_df.index=(daily_df['level_0'].apply(lambda x: dt.datetime(x,1,1))+
                    daily_df['level_1'].apply(lambda x: dt.timedelta(int(x)-1)))
    daily_df.drop(['level_0','level_1'],axis=1,inplace=True)
    daily_df.columns=cols
    return daily_df
        
def cuml_sum(df):
    if isinstance(df,pd.Series):
        cols=[df.name]
    else:
        cols=df.columns
    cml_df=pd.DataFrame()
    for i in df.columns:
        cml_df[i]=df[i].cumsum()
    cml_df.columns=[i+'_cml' for i in cml_df.columns]
    return cml_df