# -*- coding: utf-8 -*-
"""
Created on Mon May 25 14:27:04 2015

@author: imchugh
"""

import pandas as pd
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
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

def plot_storage_funct_ustar():
    
    num_cats = 30   
    
    # Get data
    df, attr = get_data()
    
    # Remove daytime, missing or filled data where relevant
    sub_df = df[['Fc', 'Fc_QCFlag', 'ustar', 'ustar_QCFlag', 'Fc_storage', 
                 'Fsd', 'Fsd_QCFlag', 'Ta']]
    sub_df = sub_df[sub_df.Fc_QCFlag == 0]    
    sub_df = sub_df[sub_df.ustar_QCFlag == 0]
    sub_df = sub_df[sub_df.Fsd_QCFlag == 0]    
    sub_df = sub_df[sub_df.Fsd < 5]
    sub_df.dropna(inplace = True)    
    
    # Categorise data
    sub_df['ustar_cat'] = pd.qcut(sub_df.ustar, num_cats, labels = np.linspace(1, num_cats, num_cats))
    new_df = sub_df[['Fc', 'ustar', 'Fc_storage', 'ustar_cat', 'Ta']].groupby('ustar_cat').mean()

    # Create plot
    fig = plt.figure(figsize = (12, 8))
    fig.patch.set_facecolor('white')
    ax1 = plt.gca()
    ax2 = ax1.twinx()
    ser_1 = ax1.plot(new_df.ustar, new_df.Fc_storage, 'o', label = '$Storage$', 
                     markersize = 10, markeredgecolor = 'black', 
                     markerfacecolor = 'none', mew = 1)
    ser_2 = ax1.plot(new_df.ustar, new_df.Fc, 's', label = '$F_{c}$',
                     markersize = 10, color = 'grey')
    ser_3 = ax1.plot(new_df.ustar, new_df.Fc + new_df.Fc_storage, '^', 
                     label = '$F_{c}\/+\/Storage$',
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

def calc_plot_solar():
    
    # Get data
    df, attr = get_data()
    
    # Do calculations of daily totals
    daily_count_S = df['Ta'].groupby([lambda x: x.year, lambda y: y.month]).count()
    daily_mean_S = df['Ta'].groupby([lambda x: x.year, lambda y: y.month]).mean() * 0.0864
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

    return daily_df
    
    data_min = 364

def plot_BOM_rainfall():

    avg_int = [1971,2000]
    
    years_compare = [2011,2012,2013,2014]
    
    plot_compare = ['climatol']
    
    exclude_noQC = False
    
    df = pd.read_csv('/home/imchugh/Analysis/Whroo/Data/External/BOM_088109_precip_all.csv')
    
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
    
    fig = plt.figure(figsize = (16, 8))
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
    plt.xlabel('$Date$', fontsize = 22, labelpad = 10)
    plt.ylabel('$LAI\/(m^{2}m^{-2})$', fontsize = 22, labelpad = 10)
    plt.xticks(['2012-01-01', '2012-04-01', '2012-07-01', '2012-10-01',
                '2013-01-01', '2013-04-01', '2013-07-01', '2013-10-01',
                '2014-01-01', '2014-04-01', '2014-07-01', '2014-10-01',
                '2015-01-01', '2015-04-01'], ['Jan', 'Apr', 'Jul', 'Oct',
                                              'Jan', 'Apr', 'Jul', 'Oct',
                                              'Jan', 'Apr', 'Jul', 'Oct',
                                              'Jan', 'Apr'])
    [plt.axvline(line, color = 'black') for line in ['2013-01-01', '2014-01-01', '2015-01-01']]
    plt.legend(fontsize = 14, loc = [0.03, 0.8], numpoints = 1)
    plt.tight_layout()
    plt.show()        