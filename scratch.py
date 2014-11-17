# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 11:32:07 2014

@author: imchugh
"""
import pandas as pd
import pdb

def align_timestamps(d):
    for i,df_name in enumerate(d):
        df=d[df_name]
        if i==0:
            start_date=df.index[0]
            end_date=df.index[-1]
        else:
            if df.index[0]<start_date: start_date=df.index[0] 
            if df.index[-1]>end_date: end_date=df.index[-1]
    dt_index=pd.date_range(start=start_date,end=end_date,freq='30Min')
    for df_name in d:
        df=d[df_name]
        df=df.reindex(dt_index)
        d[df_name]=df
    return d
        
def collate_dfs(d):
    for i,df_name in enumerate(d):
        if i==0:
            cols=d[df_name].columns
            Tout_df=pd.DataFrame(index=d[df_name].index,columns=d.keys())
            Dout_df=pd.DataFrame(index=d[df_name].index,columns=d.keys())
        Tout_df[df_name]=d[df_name][cols[0]]
        Dout_df[df_name]=d[df_name][cols[1]]
    return Tout_df,Dout_df
        
        