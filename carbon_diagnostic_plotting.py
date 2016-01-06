# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 09:28:46 2015

@author: imchugh
"""

import pandas as pd
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../Partitioning')
import partn_wrapper as partn
import pdb

import DataIO as io

reload (partn)

def get_data():
    
    reload(io)    
    
    file_in = '/home/imchugh/Ozflux/Sites/Whroo/Data/Processed/all/Whroo_2011_to_2014_L6.nc'
    
    return io.OzFluxQCnc_to_pandasDF(file_in)

def plot_storage_all_levels_funct_ustar(correct_storage = False):
    """
    This script plots all levels of the profile system as a function of ustar;
    the decline in storage of the levels below 8m with decreasing u* is thought
    to be associated with horizontal advection due to drainage currents, so 
    there is a kwarg option to output additional series where the low level storage 
    estimates are increased by regressing those series on the 8-16m level 
    (between u* = 0.4 and u* = 0.2m s-1) and extrapolating the regression to 
    0.2-0m s-1;
    """
    
    num_cats = 30   
    
    # Get data
    df, attr = get_data()
    lt_df = partn.main()[1]
    df['Fre_lt'] = lt_df['Nocturnally derived Re']
    
    # Make variable lists
    storage_vars = ['Fc_storage', 'Fc_storage_1', 'Fc_storage_2', 
                    'Fc_storage_3', 'Fc_storage_4', 'Fc_storage_5', 
                    'Fc_storage_6',]
    anc_vars = ['ustar', 'ustar_QCFlag', 'Fsd', 'Ta', 'Fc', 'Fre_lt']    
    var_names = ['0-36m', '0-0.5m', '0.5-2m', '2-4m', '4-8m', '8-16m', '16-36m']    
    
    # Remove daytime, missing or filled data where relevant
    sub_df = df[storage_vars + anc_vars]
    sub_df = sub_df[sub_df.ustar_QCFlag == 0]    
    sub_df = sub_df[sub_df.Fsd < 5]
    sub_df.drop(['ustar_QCFlag', 'Fsd'], axis = 1, inplace = True)
    sub_df.dropna(inplace = True)    
    
    # Add a variable
    sub_df['Re_take_Fc'] = sub_df['Fre_lt'] - sub_df['Fc']
    
    # Categorise data
    sub_df['ustar_cat'] = pd.qcut(sub_df.ustar, num_cats, 
                                  labels = np.linspace(1, num_cats, num_cats))
    means_df = sub_df.groupby('ustar_cat').mean()

    # Calculate uncertainty
    error_df = (sub_df[['Fc_storage', 'Re_take_Fc', 'ustar_cat']]
                .groupby('ustar_cat').std() / 
                np.sqrt(sub_df[['Fc_storage', 'Re_take_Fc', 'ustar_cat']]
                        .groupby('ustar_cat').count())) * 2

    # Generate regression statistics and apply to low u* data for low levels 
    # (0-0.5, 0.5-2, 2-4, 4-8) if correct_storage = True
    if correct_storage:
        stats_df = pd.DataFrame(columns=storage_vars[1:5], index = ['a', 'b'])
        for var in stats_df.columns:
            coeffs = np.polyfit(means_df['Fc_storage_5'][(means_df.ustar < 0.4) & 
                                                         (means_df.ustar > 0.2)],
                                means_df[var][(means_df.ustar < 0.4) & 
                                              (means_df.ustar > 0.2)],
                                1)
            stats_df.loc['a', var] = coeffs[0]
            stats_df.loc['b', var] = coeffs[1]
    
        corr_df = means_df.copy()
        for var in stats_df.columns:
            corr_df[var][corr_df['ustar'] < 0.2] = (corr_df['Fc_storage_5']
                                                    [corr_df['ustar'] < 0.2] * 
                                                    stats_df.loc['a', var] + 
                                                    stats_df.loc['b', var])
        corr_df['Fc_storage'] = corr_df[storage_vars[1:]].sum(axis = 1)
        error_df['Fc_storage_corrected'] = error_df['Fc_storage']
        error_df['Fc_storage_corrected'][corr_df['ustar'] < 0.2] = (
            error_df['Fc_storage'] * 
            corr_df[storage_vars[0]][corr_df['ustar'] < 0.2] /
            means_df[storage_vars[0]][means_df['ustar'] < 0.2])
        corr_df = corr_df[corr_df.ustar < 0.25]
        means_df['Fc_storage_corrected'] = pd.concat(
            [corr_df[storage_vars[0]][means_df.ustar < 0.25],
             means_df[storage_vars[0]][means_df.ustar > 0.25]])
    
    # Create plot
    fig = plt.figure(figsize = (12, 8))
    fig.patch.set_facecolor('white')
    ax1 = plt.gca()
    colour_idx = np.linspace(0, 1, 6)
    vars_dict = {'Fc_storage_corrected': 'blue', 'Re_take_Fc': 'grey'}
    ax1.plot(means_df.ustar, means_df[storage_vars[0]], label = var_names[0], 
             color = 'blue')
    for i, var in enumerate(storage_vars[1:]):
        ax1.plot(means_df.ustar, means_df[var], color = plt.cm.cool(colour_idx[i]), 
                 label = var_names[i + 1])
    if correct_storage:
        ax1.plot(corr_df.ustar, corr_df[storage_vars[0]], color = 'blue', 
                 linestyle = '--')
        ax1.plot(means_df.ustar, means_df.Re_take_Fc, label = '$\widehat{ER}\/-\/F_{c}$',
                 linestyle = '-', color = 'black')
        for var in vars_dict.keys():
            x = means_df.ustar         
            y1 = means_df[var]
            y2 = y1 + error_df[var]
            y3 = y1 - error_df[var]
            ax1.fill_between(x, y1, y2, where = y2 >= y1, facecolor=vars_dict[var], 
                             edgecolor='None', alpha=0.3, interpolate=True)
            ax1.fill_between(x, y1, y3, where = y3 <= y1, facecolor=vars_dict[var], 
                             edgecolor='None', alpha=0.3, interpolate=True)                
        x = corr_df.ustar         
        y1 = means_df[storage_vars[0]][means_df.ustar < 0.25]
        y2 = corr_df[storage_vars[0]]
        ax1.fill_between(x, y1, y2, where = y2 >= y1, hatch = '.', alpha = 0.5,
                         interpolate=True, edgecolor = 'black', color='None',
                         linewidth = 0.0)
        for i, var in enumerate(storage_vars[1:]):
            ax1.plot(corr_df.ustar, corr_df[var], color = plt.cm.cool(colour_idx[i]), 
                     linestyle = '--')
    ax1.axhline(y = 0, color  = 'black', linestyle = '-')
    ax1.set_ylabel(r'$C\/flux\/(\mu mol C\/m^{-2} s^{-1})$', fontsize = 22)
    ax1.set_xlabel('$u_{*}\/(m\/s^{-1})$', fontsize = 22)
    ax1.tick_params(axis = 'x', labelsize = 14)
    ax1.tick_params(axis = 'y', labelsize = 14)
    plt.setp(ax1.get_yticklabels()[0], visible = False)
    ax1.legend(fontsize = 16, loc = [0.78,0.5], numpoints = 1)
    plt.tight_layout()   
    fig.savefig('/media/Data/Dropbox/Work/Manuscripts in progress/Writing/Whroo ' \
                'basic C paper/Images/ustar_vs_storage.png',
                bbox_inches='tight',
                dpi = 300) 
    plt.show()
    
    return

def plot_storage_Fc_funct_ustar(correct_storage = False):
    
    num_cats = 30   
    
    # Get data
    df, attr = get_data()

    # Make variable lists
    use_vars = ['Fc_storage', 'Fc_storage_1', 'Fc_storage_2', 'Fc_storage_3', 
                'Fc_storage_4', 'Fc_storage_5', 'Fc_storage_6', 
                'Fc', 'Fc_QCFlag', 'ustar', 'ustar_QCFlag', 'Fsd', 'Ta', 'Sws']
    
    # Remove daytime, missing or filled data where relevant
    sub_df = df[use_vars]
    sub_df = sub_df[sub_df.Fc_QCFlag == 0]
    sub_df = sub_df[sub_df.ustar_QCFlag == 0]    
    sub_df = sub_df[sub_df.Fsd < 5]
    sub_df.drop(['Fc_QCFlag','ustar_QCFlag'], axis = 1, inplace = True)
    sub_df.dropna(inplace = True)    

    # Categorise data
    sub_df['ustar_cat'] = pd.qcut(sub_df.ustar, num_cats, labels = np.linspace(1, num_cats, num_cats))
    new_df = sub_df.groupby('ustar_cat').mean()
    
    # Generate regression statistics and apply to low u* data for low levels 
    # (0-0.5, 0.5-2, 2-4, 4-8) if correct_storage = True
    if correct_storage:                 
        stats_df = pd.DataFrame(columns=use_vars[1:5], index = ['a', 'b'])
        for var in stats_df.columns:
            coeffs = np.polyfit(new_df['Fc_storage_5'][(new_df.ustar < 0.4) & (new_df.ustar > 0.2)],
                                new_df[var][(new_df.ustar < 0.4) & (new_df.ustar > 0.2)],
                                1)
            stats_df.loc['a', var] = coeffs[0]
            stats_df.loc['b', var] = coeffs[1]
                         
        for var in stats_df.columns:
            new_df[var][new_df['ustar'] < 0.2] = (new_df['Fc_storage_5'][new_df['ustar'] < 0.2] * 
                                                  stats_df.loc['a', var] + stats_df.loc['b', var])
        
        new_df['Fc_storage'] = new_df[use_vars[1:7]].sum(axis = 1)
    
    # Create plot
    fig = plt.figure(figsize = (12, 8))
    fig.patch.set_facecolor('white')
    ax1 = plt.gca()
    ax2 = ax1.twinx()
    ser_1 = ax1.plot(new_df.ustar, new_df.Fc, 's', label = '$F_{c}$',
                     markersize = 10, color = 'grey')
    ser_2 = ax1.plot(new_df.ustar, new_df.Fc_storage, 'o', label = '$S_{c}$', 
                     markersize = 10, markeredgecolor = 'black', 
                     markerfacecolor = 'none', mew = 1)
    ser_3 = ax1.plot(new_df.ustar, new_df.Fc + new_df.Fc_storage, '^', 
                     label = '$F_{c}\/+\/S_{c}$',
                     markersize = 10, color = 'black')
    ser_4 = ax2.plot(new_df.ustar, new_df.Ta, color = 'black', linestyle = ':',
                     label = '$T_{a}$')
    ser_5 = ax2.plot(new_df.ustar, new_df.Sws * 100, color = 'black', 
                     linestyle = '-.', label = '$VWC$')                     
    ax1.axvline(x = 0.42, color  = 'black', linestyle = '--')
    ax1.axhline(y = 0, color  = 'black', linestyle = '-')
    ax1.set_ylabel(r'$C\/flux\/(\mu mol C\/m^{-2} s^{-1})$', fontsize = 22)
    ax2.set_ylabel('$T_{a}\/(^{o}C)\//\/VWC\/(m^{3}m^{-3}\/$'+'x'+'$\/10^{2})$', 
                   fontsize = 20)
    ax2.set_ylim([-5,25])
    ax1.set_xlabel('$u_{*}\/(m\/s^{-1})$', fontsize = 22)
    ax1.tick_params(axis = 'x', labelsize = 14)
    ax1.tick_params(axis = 'y', labelsize = 14)
    ax2.tick_params(axis = 'y', labelsize = 14)
    [plt.setp(ax.get_yticklabels()[0], visible = False) for ax in [ax1, ax2]]
    all_ser = ser_1 + ser_2 + ser_3 + ser_4 + ser_5
    labs = [ser.get_label() for ser in all_ser]
    ax1.legend(all_ser, labs, fontsize = 16, loc = [0.75,0.28], numpoints = 1)
    plt.tight_layout()   
    fig.savefig('/media/Data/Dropbox/Work/Manuscripts in progress/Writing/Whroo ' \
                'basic C paper/Images/ustar_vs_Fc_and_storage.png',
                bbox_inches='tight',
                dpi = 300) 
    plt.show()
    
    return

def plot_storage_Fc_advection_funct_ustar():
    
    num_cats = 50   
    
    # Get data
    df, attr = get_data()
    lt_df = partn.main()[1]
    df['Fre_lt'] = lt_df['Nocturnally derived Re']

    # Make variable lists
    use_vars = ['Fc_storage', 'Fc', 'Fre_lt', 'ustar', 'ustar_QCFlag', 
                'Fsd', 'Ta', 'Fc_QCFlag', 'Sws']
    
    # Remove daytime, missing or filled data where relevant
    sub_df = df[use_vars]
    sub_df = sub_df[sub_df.ustar_QCFlag == 0]    
    sub_df = sub_df[sub_df.Fc_QCFlag == 0]  
    sub_df = sub_df[sub_df.Fsd < 5]
    sub_df.dropna(inplace = True)    

    # Generate advection estimate
    sub_df['advection'] = sub_df['Fre_lt'] - sub_df['Fc'] - sub_df['Fc_storage']

    # Categorise data into ustar bins then do means and Sds grouped by categories
    sub_df['ustar_cat'] = pd.qcut(sub_df.ustar, num_cats, labels = np.linspace(1, num_cats, num_cats))
    means_df = sub_df.groupby('ustar_cat').mean()
    CI_df = (sub_df[['Fc','Fc_storage','Fre_lt','advection', 'ustar_cat']]
             .groupby('ustar_cat').std() / 
             np.sqrt(sub_df[['Fc','Fc_storage','Fre_lt','advection', 'ustar_cat']]
             .groupby('ustar_cat').count()) * 2)
    
    # Create plot
    fig = plt.figure(figsize = (12, 8))
    fig.patch.set_facecolor('white')
    ax1 = plt.gca()
    ser_1 = ax1.plot(means_df.ustar, means_df.Fre_lt, linestyle = '--', 
                     label = '$\hat{ER}$', color = 'black')
    ser_2 = ax1.plot(means_df.ustar, means_df.Fc, linestyle = '-', 
                     label = '$F_{c}$', color = 'black')
    ser_3 = ax1.plot(means_df.ustar, means_df.Fc_storage, linestyle = ':', 
                     label = '$S_{c}$', color = 'black')                     
    ser_4 = ax1.plot(means_df.ustar, means_df.advection, 
                     linestyle = '-', label = '$Av_{c}\/+\/Ah_{c}$', color = 'grey')
    x = means_df.ustar
    y1 = means_df.advection
    y2 = means_df.advection + CI_df.advection
    y3 = means_df.advection - CI_df.advection
    ax1.fill_between(x, y1, y2, where=y2>=y1, facecolor='0.8', edgecolor='None',
                     interpolate=True)
    ax1.fill_between(x, y1, y3, where=y3<=y1, facecolor='0.8', edgecolor='None',
                     interpolate=True)                     
    ax1.axvline(x = 0.42, color  = 'black', linestyle = '-.')
    ax1.axhline(y = 0, color  = 'black', linestyle = '-')
    ax1.set_ylabel(r'$C\/flux\/(\mu mol C\/m^{-2} s^{-1})$', fontsize = 22)
    ax1.set_xlabel('$u_{*}\/(m\/s^{-1})$', fontsize = 22)
    ax1.tick_params(axis = 'x', labelsize = 14)
    ax1.tick_params(axis = 'y', labelsize = 14)
    plt.setp(ax1.get_yticklabels()[0], visible = False)
    ax1.legend(fontsize = 16, loc = [0.7,0.28], numpoints = 1) #all_ser, labs, 
    plt.tight_layout()   
    fig.savefig('/media/Data/Dropbox/Work/Manuscripts in progress/Writing/Whroo ' \
                'basic C paper/Images/ustar_vs_Fc_and_storage_advection.png',
                bbox_inches='tight',
                dpi = 300) 
    plt.show()
    
    return

def plot_estimated_storage_and_Fc_funct_ustar():
    
    num_cats = 50   
    
    # Get data
    df, attr = get_data()
    lt_df = partn.main()[1]
    df['Fre_lt'] = lt_df['Nocturnally derived Re']

    # Make variable lists
    storage_vars = ['Fc_storage', 'Fc_storage_1', 'Fc_storage_2', 
                    'Fc_storage_3', 'Fc_storage_4', 'Fc_storage_5', 
                    'Fc_storage_6']
    anc_vars = ['ustar', 'ustar_QCFlag', 'Fsd', 'Fsd_QCFlag', 'Ta', 'Ta_QCFlag',
                'Fc', 'Fc_QCFlag', 'Fre_lt']
    var_names = ['0-32m', '0-0.5m', '0.5-2m', '2-4m', '4-8m', '8-16m', '16-32m']    
    
    # Remove daytime, missing or filled data where relevant
    sub_df = df[storage_vars + anc_vars]
    sub_df = sub_df[sub_df.ustar_QCFlag == 0]    
    sub_df = sub_df[sub_df.Fsd_QCFlag == 0]    
    sub_df = sub_df[sub_df.Fc_QCFlag == 0]  
    sub_df = sub_df[sub_df.Fsd < 5]
    sub_df.dropna(inplace = True)    

    # Categorise data
    sub_df['ustar_cat'] = pd.qcut(sub_df.ustar, num_cats, labels = np.linspace(1, num_cats, num_cats))
    new_df = sub_df[['ustar', 'Fc_storage', 'Fc_storage_1', 'Fc_storage_2', 
                     'Fc_storage_3', 'Fc_storage_4', 'Fc_storage_5', 
                     'Fc_storage_6', 'ustar_cat', 'Ta', 'Fc', 'Fre_lt']].groupby('ustar_cat').mean()
    new_df['Fc_storage_std'] = sub_df[['Fc_storage','ustar_cat']].groupby('ustar_cat').std()

    # Create plot
    fig = plt.figure(figsize = (12, 8))
    fig.patch.set_facecolor('white')
    ax1 = plt.gca()
    ax2 = ax1.twinx()
    ser_1 = ax1.plot(new_df.ustar, new_df.Fc_storage, linestyle = ':', 
                     label = '$S_{c}$', color = 'black')
    ser_2 = ax1.plot(new_df.ustar, new_df.Fc, linestyle = '-', 
                     label = '$F_{c}$', color = 'black')
    ser_3 = ax1.plot(new_df.ustar, new_df.Fre_lt - new_df.Fc, linestyle = '-', 
                     label = '$F_{c}\/u_{*}\/-\/F_{c}\/$', color = 'grey')
    ser_4 = ax1.plot(new_df.ustar, new_df.Fre_lt, linestyle = '--', 
                     label = '$F_{c}\/u_{*}$', color = 'black')
    ser_5 = ax2.plot(new_df.ustar, new_df.Ta, color = 'red', linestyle = ':',
                     label = '$Temperature$')
    x = new_df.ustar
    y1 = new_df.Fc_storage
    y2 = new_df.Fre_lt - new_df.Fc
    ax1.fill_between(x, y1, y2, where=y2>=y1, facecolor='0.8', edgecolor='None',
                     interpolate=True)
    ax1.axvline(x = 0.42, color  = 'black', linestyle = '-.')
    ax1.axhline(y = 0, color  = 'black', linestyle = '-')
    ax1.set_ylabel(r'$C\/source\/(\mu mol C\/m^{-2} s^{-1})$', fontsize = 22)
    ax2.set_ylabel('$Temperature\/(^{o}C)$', fontsize = 20)
    ax2.set_ylim([-4,20])
    ax1.set_xlabel('$u_{*}\/(m\/s^{-1})$', fontsize = 22)
    ax1.tick_params(axis = 'x', labelsize = 14)
    ax1.tick_params(axis = 'y', labelsize = 14)
    ax2.tick_params(axis = 'y', labelsize = 14)
    plt.setp(ax1.get_yticklabels()[0], visible = False)
    all_ser = ser_1 + ser_2 + ser_3 + ser_4 + ser_5
    labs = [ser.get_label() for ser in all_ser]
    ax1.legend(all_ser, labs, fontsize = 16, loc = [0.7,0.28], numpoints = 1)
    plt.tight_layout()   
    fig.savefig('/media/Data/Dropbox/Work/Manuscripts in progress/Writing/Whroo ' \
                'basic C paper/Images/ustar_vs_Fc_and_storage_advection1.png',
                bbox_inches='tight',
                dpi = 300) 
    plt.show()
    
    return

def plot_test():
    
    num_cats = 50   
    
    # Get data
    df, attr = get_data()
    lt_df = partn.main()[1]
    df['Fre_lt'] = lt_df['Nocturnally derived Re']

    # Make variable lists
    storage_vars = ['Fc_storage', 'Fc_storage_1', 'Fc_storage_2', 
                    'Fc_storage_3', 'Fc_storage_4', 'Fc_storage_5', 
                    'Fc_storage_6']
    anc_vars = ['ustar', 'ustar_QCFlag', 'Fsd', 'Fsd_QCFlag', 'Ta', 'Ta_QCFlag',
                'Fc', 'Fc_QCFlag', 'Fre_lt']
    var_names = ['0-32m', '0-0.5m', '0.5-2m', '2-4m', '4-8m', '8-16m', '16-32m']    
    
    # Remove daytime, missing or filled data where relevant
    sub_df = df[storage_vars + anc_vars]
    sub_df = sub_df[sub_df.ustar_QCFlag == 0]    
    sub_df = sub_df[sub_df.Fsd_QCFlag == 0]    
    sub_df = sub_df[sub_df.Fc_QCFlag == 0]  
    sub_df = sub_df[sub_df.Fsd < 5]
    sub_df.dropna(inplace = True)    

    # Categorise data
    sub_df['ustar_cat'] = pd.qcut(sub_df.ustar, num_cats, labels = np.linspace(1, num_cats, num_cats))
    new_df = sub_df[['ustar', 'Fc_storage', 'Fc_storage_1', 'Fc_storage_2', 
                     'Fc_storage_3', 'Fc_storage_4', 'Fc_storage_5', 
                     'Fc_storage_6', 'ustar_cat', 'Ta', 'Fc', 'Fre_lt']].groupby('ustar_cat').mean()
    new_df['Fc_storage_std'] = sub_df[['Fc_storage','ustar_cat']].groupby('ustar_cat').std()

    # Create plot
    fig = plt.figure(figsize = (8, 8))
    fig.patch.set_facecolor('white')
    ax1 = plt.gca()
    ax1.plot(new_df.ustar, new_df.Fc_storage + new_df.Fc, linestyle = ':', 
             label = '$S_{c}$', color = 'black')
    ax1.plot(new_df.ustar, new_df.Fc, linestyle = '-', 
             label = '$F_{c}$', color = 'black')
    ax1.plot(new_df.ustar, new_df.Fre_lt, linestyle = '--', 
             label = '$\widehat{ER}$', color = 'black')

    x = new_df.ustar
    y1 = new_df.Fc
    y2 = new_df.Fc_storage + new_df.Fc
    ax1.fill_between(x, y1, y2, where=y2>=y1, facecolor='0.8', edgecolor='None',
                     interpolate=True)
    y1 = new_df.Fc_storage + new_df.Fc
    y2 = new_df.Fre_lt
    ax1.fill_between(x, y1, y2, where=y2>=y1, facecolor='0.6', edgecolor='None',
                 interpolate=True)
#    ax1.axvline(x = 0.42, color  = 'black', linestyle = '-.')
    ax1.axhline(y = 0, color  = 'black', linestyle = '-')
    ax1.xaxis.set_ticks_position('bottom')
    ax1.yaxis.set_ticks_position('left')
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.set_ylabel(r'$C\/source\/(\mu mol C\/m^{-2} s^{-1})$', fontsize = 22)
    ax1.set_xlabel('$u_{*}\/(m\/s^{-1})$', fontsize = 22)
    ax1.set_xlim([0, 0.42])
    ax1.tick_params(axis = 'x', labelsize = 14)
    ax1.tick_params(axis = 'y', labelsize = 14)
    ax1.text(0.08, 1.73, '$Av_c\/+\/Ah_c$\n$=\/0.37$', fontsize = 22, 
             horizontalalignment = 'center')
    ax1.text(0.1, 1.1, '$S_c\/=\/0.63$', fontsize = 22)
    ax1.text(0.3, 0.5, '$F_c$', fontsize = 22)
    plt.setp(ax1.get_yticklabels()[0], visible = False)
    plt.tight_layout()   
    fig.savefig('/media/Data/Dropbox/Work/Manuscripts in progress/Writing/Whroo ' \
                'basic C paper/Images/ustar_vs_Fc_and_storage_advection1.png',
                bbox_inches='tight',
                dpi = 300) 
    plt.show()
    
    print new_df.Fc[new_df.ustar<0.42].mean()
    print new_df.Fre_lt[new_df.ustar<0.42].mean()
    print new_df.Fc_storage[new_df.ustar<0.42].mean()
    
    return new_df

def storage_and_ustar_example_time_series():
    
    df, attr = get_data()
    
    dates_list = ['2012-02-11 12:00:00','2012-02-19 12:00:00']    
    vars_list = ['Cc_LI840_1m', 'Cc_LI840_2m', 'Cc_LI840_4m', 
                 'Cc_LI840_8m', 'Cc_LI840_16m', 'Cc_LI840_32m']
    new_list = ['36m', '16m', '8m', '4m', '2m', '0.5m']
    tick_locs = [i for i in 
                 df.loc[dates_list[0]: dates_list[1]].index
                 if i.hour == 0 and i.minute == 0]
    tick_labs = [dt.datetime.strftime(i.date(), '%Y-%m-%d') for i in 
                 df.loc[dates_list[0]: dates_list[1]].index
                 if i.hour == 0 and i.minute == 0]

    # Create plot
    fig = plt.figure(figsize = (12, 6))
    fig.patch.set_facecolor('white')
    ax1 = plt.gca()
    ax2 = ax1.twinx()
    colour_idx = np.linspace(0, 1, 6)
    ax1.set_ylim([360,520])
    ax1.set_ylabel('$CO_{2}\/(ppm)$', fontsize = 22)

    ax1.set_xticks(tick_locs)
    ax1.set_xticklabels(tick_labs, rotation = 'vertical', fontsize = 14)
    ax1.tick_params(axis = 'y', labelsize = 14)
    
    ax2.set_ylabel('$u_{*}\/(ms^{-1})$', fontsize = 22)
    ax2.tick_params(axis = 'y', labelsize = 14)
    for i, var in enumerate(vars_list):
        ax1.plot(df.loc[dates_list[0]: dates_list[1]].index, 
                 df.loc[dates_list[0]: dates_list[1], var], 
                 label = new_list[i], color = plt.cm.cool(colour_idx[i]))
    ax2.plot(df.loc[dates_list[0]: dates_list[1]].index, 
             df.loc[dates_list[0]: dates_list[1], 'ustar'], label = 'ustar',
             linestyle = ':', color = '0.5')
    ax2.axhline(0.42, color = 'grey', linestyle = '--')
    ax1.legend(loc = [0.7, 0.72], ncol = 2)
    plt.tight_layout()
    plt.show()    

def plot_storage_diurnal_with_ustar():
    
    file_in = '/home/imchugh/Ozflux/Sites/Whroo/Data/Processed/all/Whroo_2011_to_2014_L6_stor.nc'
    
    storage_vars = ['Fc_storage', 'Fc_storage_1', 'Fc_storage_2', 
                    'Fc_storage_3', 'Fc_storage_4', 'Fc_storage_5', 
                    'Fc_storage_6',]    
    
    df = io.OzFluxQCnc_to_data_structure(file_in, 
                                         var_list = (storage_vars + 
                                                     ['ustar','Ta', 'Fc',
                                                     'Fc_storage', 'Flu',
                                                     'Fsd']), 
                                         output_structure='pandas')

    diurnal_df = df.groupby([lambda x: x.hour, lambda y: y.minute]).mean()
    diurnal_df['Fc_storage_std'] = df['Fc_storage'].groupby([lambda x: x.hour, 
                                                             lambda y: y.minute]).std()
    diurnal_df.index = np.linspace(0, 23.5, 48)
    
    storage_mean = diurnal_df.Fc_storage.mean()
    var_names = ['0-36m', '0-0.5m', '0.5-2m', '2-4m', '4-8m', '8-16m', '16-36m']
    # Create plot
    fig = plt.figure(figsize = (12, 8))
    fig.patch.set_facecolor('white')
    colour_idx = np.linspace(0, 1, 6)
    ax1 = plt.gca()
    ax2 = ax1.twinx()
    ax1.set_xlim([0, 24])
    ax1.set_xticks([0,4,8,12,16,20,24])
    ax1.tick_params(axis = 'x', labelsize = 14)
    ax1.tick_params(axis = 'y', labelsize = 14)
    ax2.tick_params(axis = 'y', labelsize = 14)
    ax1.set_xlabel('$Time (hours)$', fontsize = 22)
    ax1.set_ylabel('$S_c\/(\mu mol C\/m^{-2} s^{-1})$', fontsize = 22)
    ax2.set_ylabel('$u_{*}\/(ms^{-1})$', fontsize = 22)
    series_list = []
    for i, var in enumerate(storage_vars[1:]):
        series_list.append(ax1.plot(diurnal_df.index, diurnal_df[var], 
                                    color = plt.cm.cool(colour_idx[i]), 
                                    label = var_names[i + 1]))
    series_list.append(ax1.plot(diurnal_df.index, diurnal_df.Fc_storage, 
                                color = '0.5', label = var_names[i + 1]))
    series_list.append(ax2.plot(diurnal_df.index, diurnal_df.ustar, 
                                color = 'black', linestyle = ':', 
                                label = '$u_*$'))
    ax1.axhline(storage_mean, color = 'black')
    ax2.axhline(0.42, linestyle = '--', color = 'black')    
    plt.setp(ax1.get_yticklabels()[0], visible = False)
    plt.setp(ax2.get_yticklabels()[0], visible = False)
    labs = [ser[0].get_label() for ser in series_list]
    lst = [i[0] for i in series_list]
    ax1.legend(lst, labs, fontsize = 16, loc = [0.035,0.75], 
               numpoints = 1, ncol = 2)
    plt.tight_layout()
    fig.savefig('/media/Data/Dropbox/Work/Manuscripts in progress/Writing/Whroo ' \
                'basic C paper/Images/diurnal_storage.png',
                bbox_inches='tight',
                dpi = 300) 
    plt.show()

def plot_Fc_storage_diurnal_with_ustar_Fsd():
    
    file_in = '/home/imchugh/Ozflux/Sites/Whroo/Data/Processed/all/Whroo_2011_to_2014_L6.nc'
        
    df = io.OzFluxQCnc_to_data_structure(file_in, 
                                         var_list = ['ustar','Ta', 'Fc',
                                                     'Fc_storage', 'Flu',
                                                     'Fsd'], 
                                         output_structure='pandas')

    df.Fsd = df.Fsd * 4.6 * 0.46
    diurnal_df = df.groupby([lambda x: x.hour, lambda y: y.minute]).mean()
    diurnal_df.index = np.linspace(0, 23.5, 48)
    
    storage_mean = diurnal_df.Fc_storage.mean()

    # Create plot
    fig = plt.figure(figsize = (12, 8))
    fig.patch.set_facecolor('white')
    ax1 = plt.gca()
    ax2 = ax1.twinx()
    ax1.set_xlim([8, 16])
    ax1.set_ylim([-0.01,-0.002])
    ax2.set_ylim([0,1])
    ax1.set_xticks([8,9,10,11,12,13,14,15,16])
    ax1.tick_params(axis = 'x', labelsize = 14)
    ax1.tick_params(axis = 'y', labelsize = 14)
    ax2.tick_params(axis = 'y', labelsize = 14)
    ax1.set_xlabel('$Time (hours)$', fontsize = 22)
    ax1.set_ylabel('$RUE\/(\mu mol C\/\mu mol\/photon^{-1}\/m^{-2} s^{-1})$', fontsize = 22)
    ax2.set_ylabel('$S_c\//\/F_c$', fontsize = 22)
    series_list = []
    series_1 = ax1.plot(diurnal_df.index, (diurnal_df.Fc + 
                                           diurnal_df.Fc_storage) / diurnal_df.Fsd, 
                        color = '0.5', label = '$S_c$')
    series_2 = ax1.plot(diurnal_df.index, diurnal_df.Fc / diurnal_df.Fsd, 
                        color = 'black', label = '$F_c\/+\/S_c$')
    series_3 = ax2.plot(diurnal_df.index, diurnal_df.Fc_storage / diurnal_df.Fc, 
                        color = 'black', label = '$S_c\//\/F_c$', linestyle = ':')
    series_list = series_1 + series_2 + series_3                        
    labs = [ser.get_label() for ser in series_list]
    ax1.legend(series_list, labs, fontsize = 18, loc = [0.75,0.77], 
               numpoints = 1, frameon = False)    
    plt.setp(ax1.get_yticklabels()[0], visible = False)
    plt.setp(ax2.get_yticklabels()[0], visible = False)
    plt.tight_layout()
    fig.savefig('/media/Data/Dropbox/Work/Manuscripts in progress/Writing/Whroo ' \
                'basic C paper/Images/RUE_and_Sc_on_Fc.png',
                bbox_inches='tight',
                dpi = 300) 
    plt.show()
    
def storage_and_T_by_wind_sector():

    # Assign storage and met variables
    stor_var = 'Fc_storage_6'
    T_var = 'Ta_HMP_32m'
    ws_var = 'Ws_RMY_32m'
    wd_var = 'Wd_RMY_32m'
    
    # Wind sectors
    wind_sectors_dict = {'NNE': [0,45], 'ENE': [45,90], 'ESE': [90,135],
                         'SSE': [135,180], 'SSW': [180,225], 'WSW': [225,270],
                         'WNW': [270,315], 'NNW': [315,360]}    
    
    # Get data, then subset to exclude extraneous / bad data
    df, attr = get_data()
    df = df[df.Fsd < 10]
    df = df[df[ws_var] != 0]
    df = df[[stor_var, T_var, ws_var, wd_var]]
    df.dropna(inplace = True)

    # Separate out by wind sector
    results_df = pd.DataFrame(index = wind_sectors_dict.keys(), 
                              columns = ['count','stor_mean', 'T_mean'])
    for sector in wind_sectors_dict:
        results_df.loc[sector, 'count'] = len(df[(df[wd_var] > wind_sectors_dict[sector][0]) & 
                                                     (df[wd_var] < wind_sectors_dict[sector][1])])
        results_df.loc[sector, 'stor_mean'] = df[stor_var][(df[wd_var] > wind_sectors_dict[sector][0]) & 
                                                         (df[wd_var] < wind_sectors_dict[sector][1])].mean()
        results_df.loc[sector, 'T_mean'] = df[T_var][(df[wd_var] > wind_sectors_dict[sector][0]) & 
                                                     (df[wd_var] < wind_sectors_dict[sector][1])].mean()                                                 

    return results_df                                                 

def temp():
    
    num_cats = 50   
    
    # Get data
    df, attr = get_data()

    # Make variable lists
    storage_vars = ['Fc_storage', 'Fc_storage_1', 'Fc_storage_2', 
                    'Fc_storage_3', 'Fc_storage_4', 'Fc_storage_5', 
                    'Fc_storage_6']
    anc_vars = ['ustar', 'ustar_QCFlag', 'Fsd', 'Fsd_QCFlag', 'Ta', 'Ta_QCFlag',
                'Fc', 'Fc_QCFlag']
    var_names = ['0-32m', '0-0.5m', '0.5-2m', '2-4m', '4-8m', '8-16m', '16-32m']    
    
    # Remove daytime, missing or filled data where relevant
    sub_df = df[storage_vars + anc_vars]
    sub_df = sub_df[sub_df.ustar_QCFlag == 0]    
    sub_df = sub_df[sub_df.Fsd_QCFlag == 0]    
    sub_df = sub_df[sub_df.Fc_QCFlag == 0]  
    sub_df = sub_df[sub_df.Fsd < 5]
    sub_df.dropna(inplace = True)    

    # Categorise data
    sub_df['ustar_cat'] = pd.qcut(sub_df.ustar, num_cats, labels = np.linspace(1, num_cats, num_cats))
    new_df = sub_df[['ustar', 'Fc_storage', 'Fc_storage_1', 'Fc_storage_2', 
                     'Fc_storage_3', 'Fc_storage_4', 'Fc_storage_5', 
                     'Fc_storage_6', 'ustar_cat', 'Ta', 'Fc']].groupby('ustar_cat').mean()
#    new_df['Fc_storage_std'] = sub_df[['Fc_storage','ustar_cat']].groupby('ustar_cat').std()

    return new_df