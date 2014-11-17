# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 17:00:07 2014

Pass dictionary containing the following parameters:
    - 'xlab' - xlabel
    - 'ylab' - ylabel
    - 'title'

@author: imchugh
"""

import matplotlib.pyplot as plt
import pdb

def line_plot(df,d):
    fig=plt.figure(figsize=(12,8))
    fig.patch.set_facecolor('white')
    for i,var in enumerate(df.columns):
        if 'colors' in d:
            plt.plot(df.index,df[var],label=var,color=d['colors'][i])
        else:
            plt.plot(df.index,df[var],label=var,linewidth=1)
    plt.title(d['title'],fontsize=24)
    plt.xlabel(d['xlab'],fontsize=18)
    plt.ylabel(d['ylab'],fontsize=18)
    plt.tick_params(labelsize=14)
    plt.xticks(rotation=45)
    if 'vert_line' in d:
        plt.axvline(x=d['vert_line'],color='black',linestyle='--')
    if 'hor_line' in d:
        plt.axhline(y=d['hor_line'],color='black',linestyle='-')
    plt.legend(loc='upper left')
    plt.show()

def one_ax_dict_set():
    """
    Returns dictionary required to run the two y axis algorithm above
    """
    d={}

    d['title']='Dendrometer stem increment\n'
    d['xlab']='\nDate'
    d['ylab']='Increment (mm)'
    d['colors']=['green','red','blue','purple','orange','brown','yellow']
    #d['vert_line']=False
    #d['hor_line']=False
    
    return d

def fill_line_plot(df,d):
    fig=plt.figure(figsize=(12,8))
    fig.patch.set_facecolor('white')
    x=df.index    
    y1=df[df.columns[0]]
    y2=df[df.columns[1]]
    plt.xlim((0,24))
    plt.plot(x,y1,color='black',linestyle='-.',linewidth=2,label=df.columns[0])
    plt.plot(x,y2,color='black',linestyle='-',linewidth=2,label=df.columns[1])
    plt.fill_between(x, y1, y2, where=y2>=y1, facecolor='blue', edgecolor='None',interpolate=True)
    plt.fill_between(x, y1, y2, where=y1>=y2, facecolor='red', edgecolor='None',interpolate=True)
    plt.xlabel(d['xlab'],fontsize=18)
    plt.ylabel(d['ylab'],fontsize=18)
    plt.title(d['title'],fontsize=24)
    plt.xticks([0,4,8,12,16,20,24])
    if 'vert_line' in d:
        plt.axvline(x=d['vert_line'],color='black',linestyle='-')
    if 'hor_line' in d:
        plt.axhline(y=d['hor_line'],color='black',linestyle='-')
    plt.legend(loc='lower right')
    
def two_ax_line_plot(df,d):
    """
    Pass dataframe (df) and dictionary (d) specified by function 
    'two_ax_dict_set' - see below;
    returns open plot
    """
    fig=plt.figure(figsize=(12,8))
    fig.patch.set_facecolor('white')
    ax = plt.gca()
    ax2 = ax.twinx()
    y1=d['y1']
    y2=d['y2']
    for i,var in enumerate(y1['vars']):
        if 'colors' in y1:
            ax.plot(df.index,df[var],label=var,color=y1['colors'][i],linewidth=1)
        else:
            ax.plot(df.index,df[var],label=var,linewidth=1)
    for i,var in enumerate(y2['vars']):
        if 'colors' in y2:
            ax2.plot(df.index,df[var],label=var,color=y2['colors'][i],linewidth=1)
    plt.title(d['globals']['title'],fontsize=24)
    ax.set_xlabel(d['globals']['xlab'],fontsize=18)
    ax.set_ylabel(y1['lab'],fontsize=18)
    ax2.set_ylabel(y2['lab'],fontsize=18)
    plt.tick_params(labelsize=14)
    ax.legend(loc='upper left',fontsize=14)
    ax2.legend(loc='upper right',fontsize=14)
    plt.show()

def two_ax_dict_set():
    """
    Returns dictionary required to run the two y axis algorithm above
    """
    d={}

    d['globals']={}
    d['globals']['title']='Fc and ET\n'
    d['globals']['xlab']='\nDate'

    d['y1']={}
    d['y1']['vars']=['GPP_SOLO_rm']
    d['y1']['colors']=['green']
    d['y1']['lab']='$GPP\/(gC\/m^{2}d^{-1})$'

    d['y2']={}
    d['y2']['vars']=['Sws']
    d['y2']['colors']=['blue']
    d['y2']['lab']='$Sws\/(m^{3}m^{-3})$'
        
    return d
    