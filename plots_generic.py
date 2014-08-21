# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 17:00:07 2014

@author: imchugh
"""

import matplotlib.pyplot as plt

def line_plot(df,title):
    fig=plt.figure(figsize=(12,8))
    fig.patch.set_facecolor('white')
    for i in df.columns:
        plt.plot(df.index,df[i],label=i)
    plt.title(title)
    plt.xlabel('Day of year',fontsize=16)
    plt.ylabel('Cumulative daily NEE $(gC\/m^{-2})$',fontsize=16)
    plt.axvline(y=0,color='black',linestyle='-')
    plt.legend(loc='upper right')
    plt.show()