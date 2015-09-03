# -*- coding: utf-8 -*-
"""
Created on Mon May 25 14:27:04 2015

@author: imchugh
"""

import pandas as pd
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import windrose
import matplotlib.cm as cm
import pdb

import DataIO as io

def get_data():
    
    reload(io)    
    
    file_in = '/home/imchugh/Ozflux/Sites/Whroo/Data/Processed/all/Whroo_2011_to_2014_L6.nc'
    
    return io.OzFluxQCnc_to_pandasDF(file_in)

def calculate_plot_NEE_diurnal():
    
    file1_in = '/home/imchugh/Ozflux/Sites/Whroo/Data/Processed/all/Whroo_2011_to_2014_L6.nc'
    file2_in = '/home/imchugh/Ozflux/Sites/Whroo/Data/Processed/all/Whroo_2011_to_2014_L6_stor.nc'
    
    df, attr = io.OzFluxQCnc_to_pandasDF(file1_in)
    
    df_stor, attr_stor = io.OzFluxQCnc_to_pandasDF(file2_in)    
    
    df = df[['Fc', 'NEE_SOLO', 'Fc_storage']]
    df.columns = ['Fc', 'NEE_SOLO_no_storage', 'Fc_storage']
    df['NEE_SOLO_incl_storage'] = df_stor['NEE_SOLO']
    
    df['NEE_SOLO_no_storage'] = df['NEE_SOLO_no_storage'] * (1/(12.0 * 1800 / 10**6))
    df['NEE_SOLO_incl_storage'] = df['NEE_SOLO_incl_storage'] * (1/(12.0 * 1800 / 10**6))
    df['Fc_plus_storage'] = df['Fc'] + df['Fc_storage']

    df.dropna(inplace = True)    

    # Do calculations of daily totals for u* uncorrected and corrected Fc and storage
    diurnal_df = df.groupby([lambda x: x.hour, lambda y: y.minute]).mean()
    diurnal_df.reset_index(inplace = True)    
    diurnal_df.drop(['level_0','level_1'], axis = 1, inplace = True)
    diurnal_df.index = np.linspace(0, 23.5, 48)

    vars_dict = {1: ['Fc', 'Fc_plus_storage'],
                 2: ['Fc', 'NEE_SOLO_no_storage'],
                 3: ['Fc_plus_storage', 'NEE_SOLO_no_storage'],
                 4: ['NEE_SOLO_no_storage', 'NEE_SOLO_incl_storage']}
                 
    names_dict = {1: ['$F$', '$F\/+\/S$'],
                  2: ['$F$', '$F\/+\/u_{*}$'],
                  3: ['$F\/+\/S$', '$F\/+\/u_{*}$'],
                  4: ['$F\/+\/u_{*}$', '$F\/+\/S\/+\/u_{*}$'] }

    lines_dict = {'Fc': '-.',
                  'Fc_plus_storage': ':',
                  'NEE_SOLO_no_storage':'-',
                  'NEE_SOLO_incl_storage': '--'}

    # Instantiate plot
    fig = plt.figure(figsize = (12, 8))
    fig.patch.set_facecolor('white')                              
    fig_labels = ['a)', 'b)', 'c)', 'd)']

    for i, var in enumerate(vars_dict.keys()):

        sbplt = 220 + i + 1
        ax = fig.add_subplot(sbplt)
        ax.set_xlim([0, 24])
        ax.set_ylim([-10, 4])
        ax.set_xticks([0,4,8,12,16,20,24])
        x = diurnal_df.index
        var1 = vars_dict[var][0]
        var2 = vars_dict[var][1]
        y1 = diurnal_df[var1]
        y2 = diurnal_df[var2]
        ax.plot(x, y1, color = 'black', linestyle = lines_dict[var1], 
                linewidth = 2, label = names_dict[var][0])
        ax.plot(x, y2, color = 'black', linestyle = lines_dict[var2], 
                linewidth = 2, label = names_dict[var][1])
        ax.fill_between(x, y1, y2, where=y2>=y1, facecolor='blue', edgecolor='None',interpolate=True)
        plt.fill_between(x, y1, y2, where=y1>=y2, facecolor='red', edgecolor='None',interpolate=True)
        if i % 2 == 0:        
            ax.set_ylabel(r'$F_{C}\/(\mu mol C\/m^{-2} s^{-1})$', fontsize = 18)
        if i > 1:
            ax.set_xlabel(r'$Time\/(hours)$', fontsize = 18)
        ax.legend(fontsize = 14, loc = 'lower right', frameon = False)
        ax.axhline(y = 0, color = 'black', linestyle = '-')
        ax.text(-2, 3.9, fig_labels[i], fontsize = 12)
    
    plt.tight_layout()   
    fig.savefig('/media/Data/Dropbox/Work/Manuscripts in progress/Writing/Whroo ' \
                'basic C paper/Images/diurnal_Fc_effects_of_storage_and_ustar.png',
                bbox_inches='tight',
                dpi = 300)     
    plt.show()
    
    return

def calculate_plot_NEE_diurnal_2():
    
    file1_in = '/home/imchugh/Ozflux/Sites/Whroo/Data/Processed/all/Whroo_2011_to_2014_L6.nc'
    file2_in = '/home/imchugh/Ozflux/Sites/Whroo/Data/Processed/all/Whroo_2011_to_2014_L6_stor.nc'
    
    df, attr = io.OzFluxQCnc_to_pandasDF(file1_in)
    
    df_stor, attr_stor = io.OzFluxQCnc_to_pandasDF(file2_in)    
    
    df = df[['Fc', 'NEE_SOLO', 'Fc_storage']]
    df.columns = ['Fc', 'NEE_SOLO_no_storage', 'Fc_storage']
    df['NEE_SOLO_incl_storage'] = df_stor['NEE_SOLO']
    df['Fc_plus_storage'] = df['Fc'] + df['Fc_storage']
    
    df.dropna(inplace = True)    

    # Do calculations of daily totals for u* uncorrected and corrected Fc and storage
    diurnal_df = df.groupby([lambda x: x.hour, lambda y: y.minute]).mean()
    diurnal_df.reset_index(inplace = True)    
    diurnal_df.drop(['level_0','level_1'], axis = 1, inplace = True)
    diurnal_df.index = np.linspace(0, 23.5, 48)

    vars_dict = {1: ['Fc', 'Fc_plus_storage'],
                 2: ['Fc_plus_storage', 'NEE_SOLO_no_storage'],
                 3: ['NEE_SOLO_no_storage', 'NEE_SOLO_incl_storage']}
                 
    names_dict = {1: ['$F$', '$F\/+\/S$'],
                  2: ['$F\/+\/S$', '$F\/+\/u_{*}$'],
                  3: ['$F\/+\/u_{*}$', '$F\/+\/S\/+\/u_{*}$'] }

    lines_dict = {'Fc': '-.',
                  'Fc_plus_storage': ':',
                  'NEE_SOLO_no_storage':'-',
                  'NEE_SOLO_incl_storage': '--'}

    # Instantiate plot
    fig = plt.figure(figsize = (15, 5))
    fig.patch.set_facecolor('white')                              
    fig_labels = ['a)', 'b)', 'c)']

    for i, var in enumerate(vars_dict.keys()):

        sbplt = 130 + i + 1
        ax = fig.add_subplot(sbplt)
        ax.set_xlim([0, 24])
        ax.set_ylim([-10, 4])
        ax.set_xticks([0,4,8,12,16,20,24])
        x = diurnal_df.index
        var1 = vars_dict[var][0]
        var2 = vars_dict[var][1]
        y1 = diurnal_df[var1]
        y2 = diurnal_df[var2]
        ax.plot(x, y1, color = 'black', linestyle = lines_dict[var1], 
                linewidth = 2, label = names_dict[var][0])
        ax.plot(x, y2, color = 'black', linestyle = lines_dict[var2], 
                linewidth = 2, label = names_dict[var][1])
        ax.fill_between(x, y1, y2, where=y2>=y1, facecolor='0.8', edgecolor='None',interpolate=True)
        plt.fill_between(x, y1, y2, where=y1>=y2, facecolor='0.8', edgecolor='None',interpolate=True)
        if i == 0:        
            ax.set_ylabel(r'$F_{C}\/(\mu mol C\/m^{-2} s^{-1})$', fontsize = 18)
        ax.set_xlabel(r'$Time\/(hours)$', fontsize = 18)
        ax.legend(fontsize = 14, loc = 'lower right', frameon = False)
        ax.axhline(y = 0, color = 'black', linestyle = '-')
        ax.text(-2.3, 3.9, fig_labels[i], fontsize = 12)
    
    plt.tight_layout()   
    fig.savefig('/media/Data/Dropbox/Work/Manuscripts in progress/Writing/Whroo ' \
                'basic C paper/Images/diurnal_Fc_effects_of_storage_and_ustar_2.png',
                bbox_inches='tight',
                dpi = 300)     
    plt.show()
    
    return

def calculate_plot_NEE_diurnal_3():
    
    file1_in = '/home/imchugh/Ozflux/Sites/Whroo/Data/Processed/all/Whroo_2011_to_2014_L6.nc'
    file2_in = '/home/imchugh/Ozflux/Sites/Whroo/Data/Processed/all/Whroo_2011_to_2014_L6_stor.nc'
    
    df, attr = io.OzFluxQCnc_to_pandasDF(file1_in)
    
    df_stor, attr_stor = io.OzFluxQCnc_to_pandasDF(file2_in)    
    
    df = df[['Fc', 'NEE_SOLO', 'Fc_storage']]
    df.columns = ['Fc', 'NEE_SOLO_no_storage', 'Fc_storage']
    df['NEE_SOLO_incl_storage'] = df_stor['NEE_SOLO']
    df['Fc_plus_storage'] = df['Fc'] + df['Fc_storage']
    
    df.dropna(inplace = True)    

    # Do calculations of daily totals for u* uncorrected and corrected Fc and storage
    diurnal_df = df.groupby([lambda x: x.hour, lambda y: y.minute]).mean()
    diurnal_df.reset_index(inplace = True)    
    diurnal_df.drop(['level_0','level_1'], axis = 1, inplace = True)
    diurnal_df.index = np.linspace(0, 23.5, 48)
    corr_stor_df = calculate_storage_levels_diurnal()
    diurnal_df['Fc_storage_corr'] = corr_stor_df['Fc_storage']
    diurnal_df['Fc_plus_storage_corr'] = diurnal_df['Fc'] + diurnal_df['Fc_storage_corr']

    vars_dict = {1: ['Fc', 'Fc_plus_storage'],
                 2: ['Fc_plus_storage', 'Fc_plus_storage_corr'],
                 3: ['Fc_plus_storage', 'NEE_SOLO_no_storage'],
                 4: ['Fc_plus_storage_corr', 'NEE_SOLO_no_storage']}
                 
    names_dict = {1: ['$F$', '$F\/+\/S$'],
                  2: ['$F\/+\/S$', '$F\/+\/S_{corr}$'],
                  3: ['$F\/+\/S$', '$F\/+\/u_{*}$'],
                  4: ['$F\/+\/S_{corr}$', '$F\/+\/u_{*}$']}

    lines_dict = {'Fc': '-.',
                  'Fc_plus_storage': ':',
                  'Fc_plus_storage_corr': '-',
                  'NEE_SOLO_no_storage': '--'}

    # Instantiate plot
    fig = plt.figure(figsize = (12, 8))
    fig.patch.set_facecolor('white')                              
    fig_labels = ['a)', 'b)', 'c)', 'd)']

    for i, var in enumerate(vars_dict.keys()):

        sbplt = 220 + i + 1
        ax = fig.add_subplot(sbplt)
        ax.set_xlim([0, 24])
        ax.set_ylim([-10, 4])
        ax.set_xticks([0,4,8,12,16,20,24])
        x = diurnal_df.index
        var1 = vars_dict[var][0]
        var2 = vars_dict[var][1]
        y1 = diurnal_df[var1]
        y2 = diurnal_df[var2]
        ax.plot(x, y1, color = 'black', linestyle = lines_dict[var1], 
                linewidth = 2, label = names_dict[var][0])
        ax.plot(x, y2, color = 'black', linestyle = lines_dict[var2], 
                linewidth = 2, label = names_dict[var][1])
        ax.fill_between(x, y1, y2, where=y2>=y1, facecolor='0.8', edgecolor='None',interpolate=True)
        plt.fill_between(x, y1, y2, where=y1>=y2, facecolor='0.8', edgecolor='None',interpolate=True)
        if i == 0:        
            ax.set_ylabel(r'$F_{C}\/(\mu mol C\/m^{-2} s^{-1})$', fontsize = 18)
        ax.set_xlabel(r'$Time\/(hours)$', fontsize = 18)
        ax.legend(fontsize = 14, loc = 'lower right', frameon = False)
        ax.axhline(y = 0, color = 'black', linestyle = '-')
        ax.text(-2.3, 3.9, fig_labels[i], fontsize = 12)
    
    plt.tight_layout()   
    fig.savefig('/media/Data/Dropbox/Work/Manuscripts in progress/Writing/Whroo ' \
                'basic C paper/Images/diurnal_Fc_effects_of_storage_and_ustar_3.png',
                bbox_inches='tight',
                dpi = 300)     
    plt.show()
    
    return
    
def calculate_plot_NEE_cml_annual():
    
    # Get data    
    df, attr = get_data()
    
    # Do calculations of daily totals
    daily_count_S = df['NEE_SOLO'].groupby([lambda x: x.year, lambda y: y.dayofyear]).count()
    daily_mean_S = df['NEE_SOLO'].groupby([lambda x: x.year, lambda y: y.dayofyear]).mean() * 86.4
    daily_mean_S[daily_count_S < 48] = np.nan
    
    # Split into years and align
    years = list(set(daily_mean_S.index.levels[0]))
    new_index = np.arange(1,367)
    daily_df = pd.concat([daily_mean_S.loc[yr].reindex(new_index) for yr in years], axis = 1)
    years_str = [str(yr) for yr in years]
    daily_df.columns = years_str
    
    # Do running mean
    new_df = pd.concat([pd.rolling_mean(daily_df[yr], 14, center = True) 
                        for yr in years_str], axis = 1)
    new_df.columns = [yr + '_rm' for yr in years_str]
    daily_df = daily_df.join(new_df)
    
    # Do cumulative plot
    new_df = pd.concat([daily_df[yr].cumsum() for yr in years_str], axis = 1)
    new_df.columns = [yr + '_cml' for yr in years_str]
    daily_df = daily_df.join(new_df)
    
    # Instantiate plot
    fig = plt.figure(figsize=(12, 8))
    fig.patch.set_facecolor('white')
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
        
    colours = ['b', 'r', 'g']
    labels = ['2012', '2013', '2014']
    
    # Do plotting
    for i, series in enumerate(daily_df.columns[1:4]):
        ax1.plot(daily_df.index, daily_df[series], linewidth = 0.3, color = colours[i])
    for i, series in enumerate(daily_df.columns[5:8]):
        ax1.plot(daily_df.index, daily_df[series], linewidth = 3, color = colours[i], label = labels[i])
    for i, series in enumerate(daily_df.columns[9:]):
        ax2.plot(daily_df.index, daily_df[series], linewidth = 3, color = colours[i], label = labels[i])
    
    # Shared properties
    tick_locs = [int(dt.datetime.strftime(dt.datetime(2010,i,1),'%j')) for i in np.arange(1,13)]
    tick_labs = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
    for ax in ax1, ax2:
        ax.axhline(y = 0, color = 'black', linestyle = '-')
        ax.set_xlim([1,366])
        ax.set_xticks(tick_locs)
        ax.set_xticklabels(tick_labs, fontsize = 12)
        ax.tick_params(axis = 'y', labelsize = 12)
    
    # ax1 specific
    ax1.set_ylabel('$NEE\/(gC\/m^{-2}\/d^{-1})$', fontsize = 18)
    ax1.legend(fontsize = 14, loc = [0.85,0.6], frameon = False)
    ax1.set_ylim([-4,4])
    ax1.text(-19, 3.9, 'a)', fontsize = '16')
    
    # ax2 specific
    ax2.set_ylabel('$NEE\/(gC\/m^{-2})$', fontsize = 18)
    ax2.set_xlabel('$Month$', fontsize = 18)
    ax2.text(-26, 90, 'b)', fontsize = '16')
       
    plt.tight_layout()   
    fig.savefig('/media/Data/Dropbox/Work/Manuscripts in progress/Writing/Whroo ' \
                'basic C paper/Images/Cumulative and daily NEE.png',
                bbox_inches='tight',
                dpi = 300) 
    plt.show()    
    
    return

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

    num_cats = 50   
    
    # Get data
    df, attr = get_data()
    
    # Make variable lists
    storage_vars = ['Fc_storage', 'Fc_storage_1', 'Fc_storage_2', 
                    'Fc_storage_3', 'Fc_storage_4', 'Fc_storage_5', 
                    'Fc_storage_6',]
    anc_vars = ['ustar', 'ustar_QCFlag', 'Fsd', 'Fsd_QCFlag', 'Ta', 'Ta_QCFlag',
                'Ta_HMP_2m', 'Ta_HMP_2m_QCFlag']    
    var_names = ['0-32m', '0-0.5m', '0.5-2m', '2-4m', '4-8m', '8-16m', '16-32m']    
    
    # Remove daytime, missing or filled data where relevant
    sub_df = df[storage_vars + anc_vars]
    sub_df = sub_df[sub_df.ustar_QCFlag == 0]    
    sub_df = sub_df[sub_df.Fsd_QCFlag == 0]    
    sub_df = sub_df[sub_df.Fsd < 5]
    sub_df.dropna(inplace = True)    
    
    # Categorise data
    sub_df['ustar_cat'] = pd.qcut(sub_df.ustar, num_cats, 
                                  labels = np.linspace(1, num_cats, num_cats))
    new_df = sub_df[['ustar', 'Fc_storage', 'Fc_storage_1', 'Fc_storage_2', 
                     'Fc_storage_3', 'Fc_storage_4', 'Fc_storage_5', 'Fc_storage_6', 
                     'ustar_cat', 'Ta', 'Ta_HMP_2m']].groupby('ustar_cat').mean()

    # Generate regression statistics and apply to low u* data for low levels 
    # (0-0.5, 0.5-2, 2-4, 4-8) if correct_storage = True
    if correct_storage:
        stats_df = pd.DataFrame(columns=storage_vars[1:5], index = ['a', 'b'])
        for var in stats_df.columns:
            coeffs = np.polyfit(new_df['Fc_storage_5'][(new_df.ustar < 0.4) & 
                                                       (new_df.ustar > 0.2)],
                                new_df[var][(new_df.ustar < 0.4) & 
                                            (new_df.ustar > 0.2)],
                                1)
            stats_df.loc['a', var] = coeffs[0]
            stats_df.loc['b', var] = coeffs[1]
    
        corr_df = new_df.copy()
        for var in stats_df.columns:
            corr_df[var][corr_df['ustar'] < 0.2] = (corr_df['Fc_storage_5']
                                                    [corr_df['ustar'] < 0.2] * 
                                                    stats_df.loc['a', var] + 
                                                    stats_df.loc['b', var])
        corr_df['Fc_storage'] = corr_df[storage_vars[1:]].sum(axis = 1)
        corr_df = corr_df[corr_df.ustar < 0.25]
        
    # Create plot
    fig = plt.figure(figsize = (12, 8))
    fig.patch.set_facecolor('white')
    ax1 = plt.gca()
    colour_idx = np.linspace(0, 1, 6)
    ax1.plot(new_df.ustar, new_df[storage_vars[0]], label = var_names[0], 
             color = 'grey')
    for i, var in enumerate(storage_vars[1:]):
        ax1.plot(new_df.ustar, new_df[var], color = plt.cm.cool(colour_idx[i]), 
                 label = var_names[i + 1])
    if correct_storage:
        ax1.plot(corr_df.ustar, corr_df[storage_vars[0]], color = 'grey', 
                 linestyle = '--')
        x = corr_df.ustar         
        y1 = new_df[storage_vars[0]][new_df.ustar < 0.25]
        y2 = corr_df[storage_vars[0]]
        ax1.fill_between(x, y1, y2, where = y2 >= y1, facecolor='0.8', 
                         edgecolor='None', interpolate=True)
        for i, var in enumerate(storage_vars[1:]):
            ax1.plot(corr_df.ustar, corr_df[var], color = plt.cm.cool(colour_idx[i]), 
                     linestyle = '--')
    ax1.axhline(y = 0, color  = 'black', linestyle = '-')
    ax1.set_ylabel(r'$C\/storage\/(\mu mol C\/m^{-2} s^{-1})$', fontsize = 22)
    ax1.set_xlabel('$u_{*}\/(m\/s^{-1})$', fontsize = 22)
    ax1.tick_params(axis = 'x', labelsize = 14)
    ax1.tick_params(axis = 'y', labelsize = 14)
    plt.setp(ax1.get_yticklabels()[0], visible = False)
    ax1.legend(fontsize = 16, loc = [0.78,0.55], numpoints = 1)
    plt.tight_layout()   
    fig.savefig('/media/Data/Dropbox/Work/Manuscripts in progress/Writing/Whroo ' \
                'basic C paper/Images/ustar_vs_storage.png',
                bbox_inches='tight',
                dpi = 300) 
    plt.show()
    
    return
    
def plot_storage_and_Fc_funct_ustar(correct_storage = False):
    
    num_cats = 30   
    
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
    
    # Generate regression statistics and apply to low u* data for low levels 
    # (0-0.5, 0.5-2, 2-4, 4-8) if correct_storage = True
    if correct_storage:                 
        stats_df = pd.DataFrame(columns=storage_vars[1:5], index = ['a', 'b'])
        for var in stats_df.columns:
            coeffs = np.polyfit(new_df['Fc_storage_5'][(new_df.ustar < 0.4) & (new_df.ustar > 0.2)],
                                new_df[var][(new_df.ustar < 0.4) & (new_df.ustar > 0.2)],
                                1)
            stats_df.loc['a', var] = coeffs[0]
            stats_df.loc['b', var] = coeffs[1]
                         
        for var in stats_df.columns:
            new_df[var][new_df['ustar'] < 0.2] = (new_df['Fc_storage_5'][new_df['ustar'] < 0.2] * 
                                                  stats_df.loc['a', var] + stats_df.loc['b', var])
        
        new_df['Fc_storage'] = new_df[storage_vars[1:]].sum(axis = 1)
    
    # Create plot
    fig = plt.figure(figsize = (12, 8))
    fig.patch.set_facecolor('white')
    ax1 = plt.gca()
    ax2 = ax1.twinx()
    ser_1 = ax1.plot(new_df.ustar, new_df.Fc_storage, 'o', label = '$S_{c}$', 
                     markersize = 10, markeredgecolor = 'black', 
                     markerfacecolor = 'none', mew = 1)
    ser_2 = ax1.plot(new_df.ustar, new_df.Fc, 's', label = '$F_{c}$',
                     markersize = 10, color = 'grey')
    ser_3 = ax1.plot(new_df.ustar, new_df.Fc + new_df.Fc_storage, '^', 
                     label = '$F_{c}\/+\/S_{c}$',
                     markersize = 10, color = 'black')
    ser_4 = ax2.plot(new_df.ustar, new_df.Ta, color = 'black', linestyle = ':',
                     label = '$Temperature$')
    ax1.axvline(x = 0.42, color  = 'black', linestyle = '-.')
    ax1.axvline(x = 0.36, color  = 'black', linestyle = '--')
    ax1.axhline(y = 0, color  = 'black', linestyle = '-')
    ax1.set_ylabel(r'$C\/source\/(\mu mol C\/m^{-2} s^{-1})$', fontsize = 22)
    ax2.set_ylabel('$Temperature\/(^{o}C)$', fontsize = 20)
    ax2.set_ylim([-4,20])
    ax1.set_xlabel('$u_{*}\/(m\/s^{-1})$', fontsize = 22)
    ax1.tick_params(axis = 'x', labelsize = 14)
    ax1.tick_params(axis = 'y', labelsize = 14)
    ax2.tick_params(axis = 'y', labelsize = 14)
    plt.setp(ax1.get_yticklabels()[0], visible = False)
    all_ser = ser_1 + ser_2 + ser_3 + ser_4
    labs = [ser.get_label() for ser in all_ser]
    ax1.legend(all_ser, labs, fontsize = 16, loc = [0.7,0.28], numpoints = 1)
    plt.tight_layout()   
    fig.savefig('/media/Data/Dropbox/Work/Manuscripts in progress/Writing/Whroo ' \
                'basic C paper/Images/ustar_vs_Fc_and_storage.png',
                bbox_inches='tight',
                dpi = 300) 
    plt.show()
    
    return

def calculate_storage_levels_diurnal():

    num_cats = 30   
    
    # Get data
    df, attr = get_data()
    
    # Make variable lists
    storage_vars = ['Fc_storage', 'Fc_storage_1', 'Fc_storage_2', 
                    'Fc_storage_3', 'Fc_storage_4', 'Fc_storage_5', 
                    'Fc_storage_6',]
    anc_vars = ['ustar', 'ustar_QCFlag', 'Fsd']    
    var_names = ['0-32m', '0-0.5m', '0.5-2m', '2-4m', '4-8m', '8-16m', '16-32m']    
    
    # Remove daytime, missing or filled data where relevant
    sub_df = df[storage_vars + anc_vars]
    sub_df = sub_df[sub_df.ustar_QCFlag == 0]    
    sub_df.dropna(inplace = True)    
    
    # Categorise data
    sub_df['ustar_cat'] = pd.qcut(sub_df.ustar, num_cats, 
                                  labels = np.linspace(1, num_cats, num_cats))
    group_df = sub_df[['ustar', 'Fc_storage', 'Fc_storage_1', 'Fc_storage_2', 
                     'Fc_storage_3', 'Fc_storage_4', 'Fc_storage_5', 'Fc_storage_6', 
                     'ustar_cat']].groupby('ustar_cat').mean()

    # Generate regression statistics and apply to low u* data for low levels 
    # (0-0.5, 0.5-2, 2-4, 4-8) if correct_storage = True
    stats_df = pd.DataFrame(columns=storage_vars[1:5], index = ['a', 'b'])
    for var in stats_df.columns:
        coeffs = np.polyfit(group_df['Fc_storage_5'][(group_df.ustar < 0.4) & 
                                                     (group_df.ustar > 0.2)],
                            group_df[var][(group_df.ustar < 0.4) & 
                                          (group_df.ustar > 0.2)],
                            1)
        stats_df.loc['a', var] = coeffs[0]
        stats_df.loc['b', var] = coeffs[1]

    corr_df = group_df.copy()
    for var in stats_df.columns:
        corr_df[var][corr_df['ustar'] < 0.2] = (corr_df['Fc_storage_5']
                                                [corr_df['ustar'] < 0.2] * 
                                                stats_df.loc['a', var] + 
                                                stats_df.loc['b', var])
    corr_df['Fc_storage'] = corr_df[storage_vars[1:]].sum(axis = 1)
    corr_df = corr_df[corr_df.ustar < 0.25]

    # Make a 3rd order polynomial to provide continuous estimation as function
    # of ustar
    params_df = pd.DataFrame(index = ['0', '1', '2', '3'], 
                             columns = storage_vars[1:5])
    for var in params_df.columns:
        params_df[var] = np.polyfit(corr_df.ustar, corr_df[var], 3)
        
    # Apply to all data on basis of nocturnal ustar
    result_df = sub_df[storage_vars + ['ustar', 'Fsd']].copy()
    for var in storage_vars[1:5]:
        result_df['series'] = (params_df.loc['0', var] * df.ustar ** 3 + 
                               params_df.loc['1', var] * df.ustar ** 2 +   
                               params_df.loc['2', var] * df.ustar + 
                               params_df.loc['3', var])
        result_df[var][(result_df.ustar < 0.22) & (result_df.Fsd < 10)] = result_df['series']
    result_df = result_df[storage_vars]    
    result_df[storage_vars[0]] = result_df[storage_vars[1:]].sum(axis = 1)
    
    # Now do diurnal averaging
    diurnal_df = (result_df[storage_vars].groupby([lambda x: x.hour, 
                                                   lambda y: y.minute]).mean())
    diurnal_df.index = np.linspace(0, 23.5, 48)

    return diurnal_df
    
#    df, attr = get_data()
#    
#    df = df[['Fc', 'Fc_QCFlag', 'Fc_storage','Fc_storage_1', 'Fc_storage_2', 
#             'Fc_storage_3', 'Fc_storage_4', 'Fc_storage_5', 'Fc_storage_6']]
#    
#    # Remove daytime, missing or filled data where relevant
#    df = df[df.Fc_QCFlag == 0]  
#    df.dropna(inplace = True)        
#    
#    # Create 
#    diurnal_df = df.groupby([lambda x: x.hour, lambda y: y.minute]).mean()
#    plt.plot(diurnal_df[diurnal_df.columns[3:]])
#    plt.show()
    
def plot_storage_hist():
    
    # Get data
    df, attr = get_data()
    
    fig = plt.figure(figsize = (12, 8))
    fig.patch.set_facecolor('white')
    
    avg = str(round(df.Fc_storage.mean(), 3))
    plt.hist(df.Fc_storage, 200, [-10, 10], color = 'grey')
    plt.axvline(x = 0, color  = 'black')
    plt.tick_params(axis = 'x', labelsize = 14)
    plt.tick_params(axis = 'y', labelsize = 14)
    plt.text(5, 3500,  
             '$Average = $' + avg + '$\mu mol C\/m^{-2} s^{-1}$', 
             fontsize = 18)    
    
    plt.show()

def calc_plot_vars():
    
    # Get data
    df, attr = get_data()
    
    monthly_mean_df = df[['Ta', 'Fsd', 'Sws']].groupby([lambda x: x.year, 
                                                        lambda y: y.month]).mean()
    
    years = list(set(monthly_mean_df.index.levels[0]))
    years_str = [str(yr) for yr in years]
    new_index = np.linspace(1, 12, 12)
    Ta_df = pd.concat([monthly_mean_df.loc[yr, 'Ta'].reindex(new_index) 
                       for yr in years], axis = 1)
    Fsd_df = pd.concat([monthly_mean_df.loc[yr, 'Fsd'].reindex(new_index) 
                        for yr in years], axis = 1)
    VWC_df = pd.concat([monthly_mean_df.loc[yr, 'Sws'].reindex(new_index) 
                        for yr in years], axis = 1)
    Fsd_df = Fsd_df * 0.0864
    Ta_df.columns = years_str
    Fsd_df.columns = years_str
    VWC_df.columns = years_str

    fig = plt.figure(figsize = (16, 6))
    fig.patch.set_facecolor('white')
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)

    lines = ['-', '--', ':']
    colours = ['0.5', '0', '0']

    for ax in [ax1, ax2, ax3]:
        ax.set_xlim([0, 11])
        ax.set_xlabel('$Month$', fontsize = 18, labelpad = 10)
        ax.set_xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
        ax.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 
                            'N', 'D'], fontsize = 14)
        ax.tick_params(axis = 'y', labelsize = 14)
    
    ax1.set_ylabel('$Temperature\/(^{o}C)$', fontsize = 18)
    ax2.set_ylabel('$Insolation\/(MJ\/m^{-2}\/d^{-1})$', fontsize = 18)
    ax3.set_ylabel('$Volumetric\/soil\/H_{2}O\/content\/(m^{2}\/m^{-2})$', 
                   fontsize = 18)
    
    for i, col in enumerate(Ta_df.columns[1:]):
        ax1.plot(Ta_df[col], label = col, linestyle = lines[i], 
                 color = colours[i], linewidth = 2)
        ax2.plot(Fsd_df[col], label = col, linestyle = lines[i], 
                 color = colours[i], linewidth = 2)
        ax3.plot(VWC_df[col], label = col, linestyle = lines[i], 
                 color = colours[i], linewidth = 2)
    ax1.legend(loc='lower right', frameon = False)
    ax2.legend(loc='lower right', frameon = False)
    ax3.legend(loc='upper right', frameon = False)
    plt.tight_layout()
    plt.show()
    
def calc_WUE():
    
    df, attr = get_data()
    
    df['ET'] = df.Fe * 1800 / 2260 / 10**3
    df.GPP_SOLO = df.GPP_SOLO
    
    daily_sums_df = df[['GPP_SOLO', 'ET', 'Precip']]. \
                    groupby([lambda x: x.year, lambda y: y.dayofyear]).sum()
    
    daily_sums_df['GPP_SOLO'] = daily_sums_df['GPP_SOLO'] * 1800 * 12 * 10**-6
    daily_sums_df['WUE'] = daily_sums_df['GPP_SOLO'] / daily_sums_df['ET']
    daily_sums_df['WUE_rm'] = pd.rolling_mean(daily_sums_df['WUE'], 14, 
                                              center = True)
                                              
    daily_sums_df['WUE_no_rain'] = daily_sums_df['WUE']
    daily_sums_df['WUE_no_rain'][daily_sums_df.Precip!=0] = np.nan
    daily_sums_df['WUE_no_rain'] = daily_sums_df['WUE_no_rain'].interpolate()
    daily_sums_df['WUE_no_rain_rm'] = pd.rolling_mean(daily_sums_df['WUE_no_rain'], 30, 
                                                      center = True)
    daily_sums_df = daily_sums_df.reset_index()
    daily_sums_df.index=(daily_sums_df['level_0'].apply(lambda x: dt.datetime(x,1,1))+
                         daily_sums_df['level_1'].apply(lambda x: dt.timedelta(int(x)-1)))
    daily_sums_df.drop(['level_0','level_1'], axis = 1, inplace = True)
    
    fig = plt.figure(figsize = (16, 6))
    fig.patch.set_facecolor('white')
    plt.xlabel('$Month$', fontsize = 18, labelpad = 10)
    plt.ylabel('$WUE\/(gC\/kg^{-1}H_{2}O)$', fontsize = 18)
    plt.tick_params(axis = 'x', labelsize = 14)
    plt.tick_params(axis = 'y', labelsize = 14)
    plt.plot(daily_sums_df.index, daily_sums_df['WUE_no_rain'], color = '0.5',
             linewidth = 0.5)
    plt.plot(daily_sums_df.index, daily_sums_df['WUE_no_rain_rm'], color = '0.3',
             linewidth = 2)
    plt.xlim(['2012-01-01','2014-12-31'])
    plt.ylim([0.5, 3.5])
    plt.xticks(['2012-01-01', '2012-04-01', '2012-07-01', '2012-10-01',
                '2013-01-01', '2013-04-01', '2013-07-01', '2013-10-01',
                '2014-01-01', '2014-04-01', '2014-07-01', '2014-10-01'],
               ['Jan', 'Apr', 'Jul', 'Oct',
                'Jan', 'Apr', 'Jul', 'Oct',
                'Jan', 'Apr', 'Jul', 'Oct',
                'Jan', 'Apr'])
    [plt.axvline(line, color = 'black') for line in ['2013-01-01', '2014-01-01']]
    plt.axhline(daily_sums_df['WUE_no_rain'].mean(), color = 'black', linestyle = ':')
    plt.text('2012-06-05', 3.2, '2012', fontsize = 22)
    plt.text('2013-06-05', 3.2, '2013', fontsize = 22)
    plt.text('2014-06-05', 3.2, '2014', fontsize = 22)
    plt.tight_layout()
    plt.show()

    return
    
def plot_BOM_rainfall():

    avg_int = [1971,2000]
    
    years_compare = [2011,2012,2013,2014]
    
    plot_compare = ['climatol']
    
    exclude_noQC = False
    
    df = pd.read_csv('/home/imchugh/Analysis/Whroo/Data/External/BOM_081043_precip.csv')
    
    df.index = [dt.datetime(df.loc[i, 'Year'],df.loc[i, 'Month'],df.loc[i, 'Day'])
                for i in df.index]
    
    cols_list = df.columns
    
    # Remove data which has not been QC'd
    if exclude_noQC:
        df[cols_list[-3]] = np.where(df[cols_list[-1]]=='Y', df[cols_list[-3]], np.nan)
    
    # Calculate amount of missing data and monthly sums for each requested year
    yrs_data_list = []
    for yr in years_compare:
    
        num_days_obs = df.loc[str(yr), cols_list[-3]].groupby(lambda x: x.year).count().loc[yr]
        num_days_yr = 365 if yr % 4 != 0 else 366
        print str(num_days_yr - num_days_obs) + ' days of observations missing from year ' + str(yr)
        temp_df = df.loc[str(yr), cols_list[-3]].groupby(lambda x: x.month).sum()
        temp_df.name = str(yr)
        yrs_data_list.append(temp_df)
    
    monthly_df = pd.concat(yrs_data_list, axis = 1)
    
    
    # Calculate standard climatology for site using specified averaging interval
    climatol = df.loc[str(avg_int[0]): str(avg_int[1]), cols_list[-3]].groupby([lambda x: x.dayofyear]).mean()
    climatol.index = [(dt.datetime(2010,1,1) + dt.timedelta(i - 1)).month for i in climatol.index]
    monthly_df['climatol_' + str(avg_int[0]) + '-' + str(avg_int[1])] = climatol.groupby(level=0).sum().round(1)
    
    # Calculate mean for all available years    
    num_records = df[cols_list[-3]].groupby([lambda x: x.year]).count()
    full_years = list(num_records[num_records >= data_min].index)
    full_df = pd.concat([df.loc[str(i)] for i in full_years])
    monthly_df['all_avail_data'] = (full_df[cols_list[-3]].groupby([lambda x: x.month]).sum() / len(full_years)).round(1)
    
    print 'The following years had required minimum number of days of obs (' + str(data_min) + '):'
    print full_years
    
    # Do plotting
    fig = plt.figure(figsize=(12,8))
    fig.patch.set_facecolor('white')
    
    width = 0.3
    x_series = np.linspace(0.5, 11.5, 12)
    var_name = [i for i in df.columns if 'climatol in i'] if plot_compare == 'climatol' else 'all_avail_data'
    LT_annual_mean = monthly_df[var_name].sum().round(1)
    
    for i, yr in enumerate(years_compare):
        
        sbplt = i + 1 + 220
        ax = fig.add_subplot(sbplt)
        
        annual_mean = monthly_df[str(yr)].sum().round(1)
        
        LT = plt.bar(x_series, monthly_df[var_name], width, color = '0.2')
        annual = plt.bar(x_series + width, monthly_df[str(yr)], width, color = '0.6')
        
        ax.set_title(str(yr), fontsize=20, y=1.03)
        
        if i == 0:
            ax.legend((LT, annual), ('1971-2000', 'year'), bbox_to_anchor=(0.5, 0.99), 
                      frameon=False, fontsize=14)
        if i == 0:
            ax.text(7, 159, 'Climatology: ' + str(LT_annual_mean) + 'mm', fontsize=14)
            ax.text(7, 143, 'Annual: ' + str(annual_mean) + 'mm', fontsize=14)
        else:
            ax.text(8, 159, 'Annual: ' + str(annual_mean) + 'mm', fontsize=14)        
                
        ax.set_xlim([0, 12.6])
        ax.set_xticks(x_series+width)
        if i > 1:    
            ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
        else:
            ax.xaxis.set_visible(False)
        ax.set_ylim([0, 180])    
        if i%2 == 0:    
            ax.set_ylabel('Rainfall (mm)', fontsize=16)
    
        plt.setp(ax.get_xticklabels(), fontsize=14)
        plt.setp(ax.get_yticklabels(), fontsize=12)
        
    #plt.figlegend((LT, annual), ('1971-2000', 'year'), 'center', ncol=2, frameon = False, fontsize=10)
    plt.tight_layout()
    plt.savefig('/home/imchugh/Analysis/Whroo/Images/mangalore_annual_precip_2011_2014.png', bbox_inches='tight')
    
def plot_LAI():

    file_in = '/media/Data/Dropbox/Data_sites non flux/MODIS_cutout_timeseries/' \
               'Whroo.MOD15A2.Lai_1km.dat'
    
    hemi_LAI_dict = {'2012-04-04': 0.82,
                     '2012-07-17': 0.90,
                     '2012-10-16': 0.97,
                     '2013-01-29': 0.94,
                     '2013-04-24': 0.82,
    #                 '2013-11-28': 0.52,
                     '2014-07-23': 1.01,
                     '2015-03-18': 1.01,
    #                 '2015-06-25': 0.57
                     }
                     
    LAI2000_dict = {'2012-07-17': 0.90,
                    '2012-10-16': 0.97,
                    '2013-11-28': 0.81,
                    '2014-05-02': 0.95,
                    '2014-06-11': 0.96,
    #                '2014-07-23': 0.85,
    #                '2014-11-12': 0.71
                    }
    
    df = pd.read_csv(file_in, skiprows = [0, 2, 3], parse_dates = ['TIMESTAMP'], 
                     index_col = 'TIMESTAMP', na_values=-9999)
                     
    short_df = df.loc['2012':].copy()
    short_df[short_df.N_AVE < 9] = np.nan
    new_index = pd.date_range(short_df.index[0], short_df.index[-1], freq = 'D')
    short_df = short_df.reindex(new_index)
    short_df['Lai_1km'] = short_df['Lai_1km'].interpolate()
    short_df['Lai_1km_run'] = pd.rolling_mean(short_df.Lai_1km, 30, center = True)
    short_df['hemi_cam'] = np.nan
    
    for date in hemi_LAI_dict.keys():
        short_df.loc[date, 'hemi_cam'] = hemi_LAI_dict[date]
    short_df['LAI2000'] = np.nan
    for date in LAI2000_dict.keys():
        short_df.loc[date, 'LAI2000'] = LAI2000_dict[date]
    
    fig = plt.figure(figsize = (16, 6))
    fig.patch.set_facecolor('white')
    
    plt.plot(short_df.index, short_df.Lai_1km, color = '0.6', linewidth = 0.5)
    plt.plot(short_df.index, short_df.Lai_1km_run, color = '0.4', linewidth = 2, 
             label = 'MOD15A2')
    plt.plot(short_df.index, short_df.hemi_cam, 'o', mew = 1.5,
             markersize = 12, markeredgecolor = 'black', markerfacecolor = 'none',
             label = 'DHP')
    plt.plot(short_df.index, short_df.LAI2000, '^', markeredgewidth = 1.5,
             markersize = 12, markeredgecolor = 'black', markerfacecolor = 'none',
             label = 'LAI2200')
    plt.xlim(['2012-01-01',short_df.index[-1]])
    plt.tick_params(axis = 'y', labelsize = 14)
    plt.tick_params(axis = 'x', labelsize = 14)
    plt.xlabel('$Month$', fontsize = 22, labelpad = 10)
    plt.ylabel('$LAI\/(m^{2}m^{-2})$', fontsize = 22, labelpad = 10)
    plt.xticks(['2012-01-01', '2012-04-01', '2012-07-01', '2012-10-01',
                '2013-01-01', '2013-04-01', '2013-07-01', '2013-10-01',
                '2014-01-01', '2014-04-01', '2014-07-01', '2014-10-01',
                '2015-01-01', '2015-04-01'], ['Jan', 'Apr', 'Jul', 'Oct',
                                              'Jan', 'Apr', 'Jul', 'Oct',
                                              'Jan', 'Apr', 'Jul', 'Oct',
                                              'Jan', 'Apr'])
    [plt.axvline(line, color = 'black') for line in ['2013-01-01', '2014-01-01', '2015-01-01']]
    plt.legend(fontsize = 18, loc = [0.14, 0.74], numpoints = 1, frameon = False)
    plt.text('2012-01-20', 1.24, '2012', fontsize = 22)
    plt.text('2013-01-20', 1.24, '2013', fontsize = 22)
    plt.text('2014-01-20', 1.24, '2014', fontsize = 22)
    plt.text('2015-01-20', 1.24, '2015', fontsize = 22)
    yticks = plt.gca().yaxis.get_major_ticks()    
    yticks[0].label1.set_visible(False)    
    plt.tight_layout()
    plt.show()
    
def wind_roses():
    
    # Program error corrections
    sonic_az = 282
    program_az = 214.5
    correction = sonic_az - program_az
    
    # Get data
    df, attr = get_data()

    # Create variables lists
    wind_spd_list = ['Ws_RMY_1m', 
                     'Ws_RMY_2m', 
                     'Ws_RMY_4m', 
                     'Ws_RMY_8m',
                     'Ws_RMY_16m', 
                     'Ws_RMY_32m',
                     'Ws_CSAT']
    wind_dir_list = ['Wd_RMY_1m', 
                     'Wd_RMY_2m', 
                     'Wd_RMY_4m', 
                     'Wd_RMY_8m',
                     'Wd_RMY_16m', 
                     'Wd_RMY_32m',
                     'Wd_CSAT']

    # Remove bad data
    for var in wind_spd_list:
        df[var] = np.where(df[var + '_QCFlag'] == 0, df[var], np.nan)
        df.drop(var + '_QCFlag', axis = 1, inplace = True)
    for var in wind_dir_list:
        df[var] = np.where(df[var + '_QCFlag']==0, df[var], np.nan)    
        df.drop(var + '_QCFlag', axis = 1, inplace = True)

    # Subset df columns    
    df = df[wind_spd_list + wind_dir_list + ['ustar', 'Fsd']]

    # Subset to night    
    df = df[df.Fsd < 10]
    df.drop('Fsd', axis = 1, inplace = True)
    df.dropna(axis = 0, inplace = True)
    
    # Correct CSAT wind direction for incorrect angle in program
    df['Wd_CSAT'] = df['Wd_CSAT'] + correction
    df['Wd_CSAT'] = np.where(df['Wd_CSAT'] > 360, df['Wd_CSAT'] - 360, 
                             df['Wd_CSAT'])

    # Correct 1m windspeed for error

    ustar_quant = df['ustar'].quantile(0.1)
    print '20% of u_star values below ' + str(round(ustar_quant, 3))

    low_ustar_df = df[df.ustar < ustar_quant]

    # Calculate frequencies for wind directions
    low_rslt_df = pd.DataFrame(columns = wind_dir_list, index = range(0, 16))
    for var in low_rslt_df.columns:
        for sector in xrange(0, 16):
            lower_bound = sector * 22.5 - 12.25
            upper_bound = lower_bound + 22.5
            if sector == 0: lower_bound = 360 + lower_bound
            if sector == 0:
                low_rslt_df[var].iloc[sector] = len(low_ustar_df[var][(low_ustar_df[var] >= lower_bound) | 
                                               (low_ustar_df[var] < upper_bound)])
            else:
                low_rslt_df[var].iloc[sector] = len(low_ustar_df[var][(low_ustar_df[var] >= lower_bound) & 
                                               (low_ustar_df[var] < upper_bound)])

    low_rslt_df = low_rslt_df / low_rslt_df.sum() * 100
    
    theta_1 = np.array([348.75])
    theta_2 = np.linspace(11.25, 326.25, 15)
    theta = np.radians(np.append(theta_1, theta_2))
    width = np.pi/ 8
    
    fig = plt.figure(figsize = (3, 21))
    fig.patch.set_facecolor('white')
    sbplt = 711
    for i, var in enumerate(low_rslt_df.columns):
        ax = plt.subplot(sbplt + i, polar = True)
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        ax.bar(theta, low_rslt_df[var], width = width, alpha=0.5,
               facecolor='0.5', edgecolor = 'None')
        plt.setp(ax.get_xticklabels(), visible = False)
        [plt.axvline(line, color = '0.5') for line in np.radians([0, 90, 180, 270])]
        ax.set_ylim([0, 25])
    
    plt.tight_layout()
#    # Plot
#    fig = plt.figure(figsize=(8, 8), dpi=80, facecolor='w', edgecolor='w')
#    rect = [0.1, 0.1, 0.8, 0.8]
#    ax = windrose.WindroseAxes(fig, rect, axisbg='w')
#    print type(ax)
#    fig.add_axes(ax)
##    ax.contourf(df['Wd_CSAT'], df['Ws_CSAT'], bins = np.arange(0,8,1), cmap=cm.hot)
#    ax.bar(sub_df['Wd_CSAT'], sub_df['Ws_CSAT'], normed=True, opening=0.8, 
#           edgecolor='white')
##    l = ax.legend(borderaxespad=-0.10)
##    plt.setp(l.get_texts(), fontsize=8)
    
    return df
    
