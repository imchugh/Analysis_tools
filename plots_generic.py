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
    fig=plt.figure(figsize=(16,8))
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
        for i in d['vert_line']:
            plt.axvline(x=i,color='black',linestyle='dotted',linewidth=0.5)
    if 'hor_line' in d:
        for i in d['hor_line']:
            plt.axhline(y=i,color='black',linestyle='-')
    plt.legend(loc='lower left')
    plt.show()

def line_plot_w_points(df,d):
    
    fig=plt.figure(figsize=(16,8))
    fig.patch.set_facecolor('white')
    mark_count=0
    for i,var in enumerate(df.columns):
        if not var in d['as_markers']:
            plt.plot(df.index,df[var],label=var,color=d['colors'][i],linewidth=1)
        else:
            plt.plot(df.index,df[var],marker=d['marker_style'][mark_count],label=var,color=d['colors'][i],markersize=10)
            mark_count +=1
    plt.title(d['title'],fontsize=24)
    plt.xlabel(d['xlab'],fontsize=18)
    plt.ylabel(d['ylab'],fontsize=18)
    plt.tick_params(labelsize=14)
    plt.xticks(rotation=45)
    for i in d['vert_line']:
            plt.axvline(x=i,color='black',linestyle='dotted',linewidth=0.5)
    for i in d['hor_line']:
            plt.axhline(y=i,color='black',linestyle='-')
    plt.legend(loc='lower right')
    plt.show()    

def one_ax_dict_set():
    """
    Returns dictionary required to run the two y axis algorithm above
    """
    d={}

    d['title']='Whroo LAI (LAI2200, Canon hemi and MODIS)\n'
    d['xlab']='\nDate'
    d['ylab']= '$LAI\/(m^{2}m^{-2})$'# '$NEE\/(gC\/m^{-2}d^{-1})$'
    d['colors'] = ['b','r','g','m','c','y','k'] # ['b','g','r','c','m','y','k'] # matplotlib default
    d['as_markers']=['Lai_Canon_hemi','Lai_LAI2000']
    d['marker_style']=['o','s']
    d['vert_line']=[]#['2012-04-01','2012-07-01','2012-10-01','2013-01-01','2013-04-01','2013-07-01','2013-10-01']
    d['hor_line']=[]
    
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
    if 'vert_line' in d['globals']:
        for i in d['globals']['vert_line']:
            plt.axvline(x=i,color='black',linestyle='dotted',linewidth=0.5)
    if 'hor_line' in d['globals']:
        for i in d['globals']['hor_line']:
            plt.axhline(y=i,color='black',linestyle='-')        
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
    d['globals']['title']='GPP and ET Riggs\n'
    d['globals']['xlab']='\nDate'
    d['globals']['vert_line']=['2012-04-01','2012-07-01','2012-10-01','2013-01-01','2013-04-01','2013-07-01','2013-10-01']

    d['y1']={}
    d['y1']['vars']=['r_GPP_rm']
    d['y1']['colors']=['green']
    d['y1']['lab']='$GPP\/(gC\/m^{2}d^{-1})$'

    d['y2']={}
    d['y2']['vars']=['r_ET_rm']
    d['y2']['colors']=['blue']
    d['y2']['lab']='$ET\/(gH_{2}O\/m^{-2}d^{-1})$'
        
    return d

def stacked_bar_plot(df,d):
    """    
    Pass dataframe (df) and dictionary (d) specified by function 
    'stacked_bar_dict_set' - see below;
    returns open plot
    """
    fig=plt.figure(figsize=(12,8))
    fig.patch.set_facecolor('white')
    ax = plt.gca()
    width=0.35
    plt.xlim(0.85,1.5)
    sum_vars=[]
    for i,var in enumerate(df.columns):
        sum_vars==[] if i==0 else sum_vars.append(df.ix[0,i-1])
        if 'colors' in d:
            plt.bar(df.index+1, df[var], width, color=d['colors'][i], label=var,bottom=sum(sum_vars))
    ax.axes.get_xaxis().set_ticks([])        
    plt.title(d['title'],fontsize=24)
    plt.ylabel(d['ylab'],fontsize=18)
    if 'xlab' in d: plt.ylabel(d['xlab'],fontsize=18)
    if 'ylab' in d: plt.ylabel(d['ylab'],fontsize=18)    
#    plt.legend(loc='upper right')
    plt.text(1.175,36,'Aboveground Vegetation:\n 37.75', fontsize=18,horizontalalignment='center',verticalalignment='center')
    plt.text(1.175,15,'Litter: 5.8', fontsize=18,horizontalalignment='center',verticalalignment='center')
    plt.text(1.175,5,'Belowground vegetation:\n 10.74', fontsize=18,horizontalalignment='center',verticalalignment='center')
    bbox_props = dict(boxstyle="larrow,pad=0.3", ec="black", fc='white', lw=2)
    plt.text(1.38,11,'Soil: 1.69',fontsize=18,bbox=bbox_props)
    fig.show()
    
def stacked_bar_dict_set():
    """
    Returns dictionary required to run the two y axis algorithm above
    """
    d={}
    
    d['title']='Whroo carbon pools\n'
    d['colors']=['blue','brown','grey','green']
    d['ylab']='$C\/storage\/(tC\/ha^{-1})$'
    
    return d