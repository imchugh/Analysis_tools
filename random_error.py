# -*- coding: utf-8 -*-
import os
import math
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib import gridspec
import pdb
import DataIO as io

#-----------------------------------------------------------------------------#

def main():
    """
    Pass the following arguments: 1) df containing non-gap filled QC'd Fc, Fsd, Ta, ws
                                  2) config file (dictionary) containing options as 
                                      specified below
    
    Returns the linear regression statistics for daytime and nocturnal data
    
    Dictionary should contain the following:
        'results_out_path' (path for results to be written to)
        'flux_freq' (frequency of flux data)
        'averaging_bins' (for calculation of random error variance as function of flux magnitude)    
        'radiation_difference_threshold' (the maximum allowable Fsd difference between Fc pairs;
                                          see below for threshold values)
        'temperature_difference_threshold' (as above)
        'windspeed_difference_threshold' (as above)
                                     
    Algorithm from reference:
        Hollinger, D.Y., Richardson, A.D., 2005. Uncertainty in eddy covariance measurements 
        and its application to physiological models. Tree Physiol. 25, 873–885.
    
    Uses daily differencing procedure to estimate random error (note this will overestimate by factor of up to 2)
    
    Fc pairs must pass difference constraints as follows (as suggested in ref above):
        Fsd:35W m^-2
        Ta: 3C
        ws: 1m s^-1
    """    
    print 'Opening data file...'
    
    # Get config file name then go get it
    config_file_in = io.file_select_dialog()
    configs = io.read_config_file(config_file_in)
    
    # Open the file and convert from .nc to pandas df
    data_file_in = configs['files']['data_file_in']
    results_out_path = configs['files']['results_out_path']
    df, d_attr = io.OzFluxQCnc_to_pandasDF(data_file_in)    

    # Calculate the sigma_delta regression coefficients
    stats_df, error_df = calculate_sigma_delta(configs, df)
    
    # Report and output results
    print '\nRegression results (also saved at ' + results_out_path+'):'
    print stats_df
    stats_df.to_csv(os.path.join(results_out_path,'random_error_stats.csv'))        
    
    # Check whether user wants error propagation and single instance error estimate
    propagate_error = configs['options']['propagate_error']=='True'
    error_generator = configs['options']['error_generator']=='True'

    if propagate_error:
        years_df = propagate_random_error(error_df, configs)

    # Report and output results    
    print '\nAnnual uncertainty in gC m-2 (also saved at ' + results_out_path + '):'
    print years_df
    years_df.to_csv(os.path.join(results_out_path, 'annual_uncertainty.csv'))
    
    if error_generator:
        random_error_generator(error_df)

    # Output data
    error_df.to_csv(os.path.join(results_out_path, 'timestep_error.csv'), header = True, index_label = 'DateTime')

    print 'Finished Processing!'

    return stats_df, error_df    
 
#-----------------------------------------------------------------------------#

def calculate_sigma_delta(configs, df):

    """

    """    
    
    print '\nCalculating random error regression coefficients'
    print '------------------------------------------------\n'
    
    # Unpack:
    # Files
    results_out_path = configs['files']['results_out_path']
    # Variables
    Fc = configs['variables']['Fc']
    Fsd = configs['variables']['Fsd']    
    Ta = configs['variables']['Ta']
    ws = configs['variables']['ws']
    ustar = configs['variables']['ustar']
    # Options    
    Ta_threshold = int(configs['options']['temperature_difference_threshold'])
    ws_threshold = int(configs['options']['windspeed_difference_threshold'])
    rad_threshold = int(configs['options']['radiation_difference_threshold'])
    num_classes = int(configs['options']['averaging_bins'])
    noct_threshold = int(configs['options']['noct_threshold'])
    records_per_day = 1440/int(configs['options']['flux_period'])
    ustar_threshold = float(configs['options']['ustar_threshold'])
        
    # Remove low ustar obs
    if not ustar_threshold == 0:
        df[Fc] = np.where(df[ustar] < ustar_threshold, np.nan, df[Fc])
    
    # Do the paired difference analysis (drop all cases containing any NaN)
    diff_df=pd.DataFrame({'Fc_mean':(df[Fc]+df[Fc].shift(records_per_day))/2,
                          'Fc_diff':df[Fc]-df[Fc].shift(records_per_day),
                		 'Ta_diff':abs(df[Ta]-df[Ta].shift(records_per_day)),
            			 'ws_diff':abs(df[ws]-df[ws].shift(records_per_day)),
            			 'Fsd_diff':abs(df[Fsd]-df[Fsd].shift(records_per_day)),
                          'Day_ind':(df[Fsd]+df[Fsd].shift(records_per_day))/2 > noct_threshold}).reset_index()
    diff_df.dropna(axis=0,how='any',inplace=True)
    
    # Find number of passed tuples, report, then drop everything except mean and difference of Fc tuples
    diff_df['pass_constraints'] = ((diff_df['Ta_diff'] < Ta_threshold) & 
                                   (diff_df['ws_diff'] < ws_threshold) &
                                   (diff_df['Fsd_diff'] < rad_threshold))
    total_tuples = len(diff_df)
    passed_tuples = len(diff_df[diff_df['pass_constraints']])	
    print (str(passed_tuples) +' of ' + str(total_tuples) + ' available tuples passed difference constraints (Fsd = ' +
  	     str(rad_threshold) + 'Wm^-2, Ta = ' + str(Ta_threshold) + 'C, ws = ' + str(ws_threshold) + 'ms^-1)\n')
    diff_df = diff_df[['Fc_mean', 'Fc_diff', 'Day_ind']][diff_df['pass_constraints']]
    diff_df['Class'] = np.nan    
    
    # Calculate and report n for each bin
    n_per_class = int(len(diff_df) / num_classes)
    print 'Number of observations per bin = ' + str(n_per_class)
    
    # Calculate and report nocturnal and daytime share of data
    day_classes = int(len(diff_df[diff_df['Day_ind']]) / float(len(diff_df)) * num_classes)
    noct_classes = num_classes - day_classes
    print 'Total bins = ' + str(num_classes) + '; day bins = ' + str(day_classes) + '; nocturnal bins = ' + str(noct_classes)

    # Cut into desired number of quantiles and get the flux and flux error means for each category
    diff_df.loc[diff_df['Day_ind'], 'Class'] = pd.qcut(diff_df.loc[diff_df['Day_ind'], 'Fc_mean'], day_classes, labels = False)
    diff_df.loc[~diff_df['Day_ind'], 'Class'] = pd.qcut(diff_df.loc[~diff_df['Day_ind'], 'Fc_mean'], noct_classes, labels = False) + day_classes
    diff_df.index = diff_df['Class'].astype('int64')
    cats_df = pd.DataFrame({'Fc_mean':diff_df['Fc_mean'].groupby(diff_df['Class']).mean()})
    cats_df['sig_del'] = np.nan
    cats_df['Day_ind'] = np.nan
    cats_df.loc[:day_classes, 'Day_ind'] = True
    cats_df.loc[day_classes:, 'Day_ind'] = False
    for i in cats_df.index:
        i = int(i)
        cats_df.loc[i, 'sig_del'] = (abs(diff_df.loc[i, 'Fc_diff'] - diff_df.loc[i, 'Fc_diff'].mean())).mean() * np.sqrt(2)
    
    # Remove daytime positive bin averages and nocturnal negative bin averages
    cats_df['exclude'] = ((cats_df.Fc_mean > 0) & (cats_df.Day_ind)) | ((cats_df.Fc_mean < 0) & (cats_df.Day_ind==False))
    cats_df = cats_df[['Fc_mean', 'sig_del', 'Day_ind']][~cats_df.exclude]
    
    # Calculate linear fit for day and night values...
    linreg_stats_df = pd.DataFrame(columns = ['slope','intcpt','r_val','p_val','SE'], index = ['day','noct'])
    linreg_stats_df.loc['day'] = stats.linregress(cats_df['Fc_mean'][cats_df['Day_ind']], cats_df['sig_del'][cats_df['Day_ind']])
    linreg_stats_df.loc['noct'] = stats.linregress(cats_df['Fc_mean'][cats_df['Day_ind'] == False], cats_df['sig_del'][cats_df['Day_ind']==False])

    # Calculate the estimated sigma_delta for each datum
    error_df = pd.DataFrame(df[Fc])
    error_df['sig_del'] = np.where(df[Fsd] >= noct_threshold,
                                   (abs(df[Fc] * linreg_stats_df.loc['noct', 'slope']) + 
                                    linreg_stats_df.loc['noct', 'intcpt']) / np.sqrt(2),
                    		     (abs(df[Fc] * linreg_stats_df.loc['day', 'slope']) +
                                    linreg_stats_df.loc['day', 'intcpt']) / np.sqrt(2))
    
    # Output plots
    print '\nPlotting: 1) PDF of random error'
    print '          2) sigma_delta (variance of random error) as a function of flux magnitude'
    fig=main_plot(diff_df['Fc_diff'], cats_df, linreg_stats_df, num_classes)
    fig.savefig(os.path.join(results_out_path, 'Random_error_plots.jpg'))    
    
    return linreg_stats_df, error_df

#-----------------------------------------------------------------------------#

def propagate_random_error(error_df, configs):
    """
        
    """

    # Unpack:
    # Variables
    Fc = configs['variables']['propagation_series']
    # Options    
    num_trials = int(configs['options']['num_trials'])
    flux_period = int(configs['options']['flux_period'])

    print '\nRunning Monte-Carlo propagation of random error to annual NEE uncertainty (n trials = ' + str(num_trials) + ')...'
    print '-----------------------------------------------------------------------------------------------'
    
    # Calculate critical t
    crit_t = stats.t.isf(0.025, num_trials) # Calculate critical t-value for p=0.095
                
    # Calculate uncertainty due to random error for all measured data at annual time step
    years_df = pd.DataFrame({'Valid_n': error_df[Fc].groupby([lambda x: x.year]).count()})
    years_df['Random error (Fc)'] = np.nan
    for i in years_df.index:
        if not len(error_df.loc[str(i), 'sig_del']) > 1:
            print '\nWarning! Insufficient data for year ' + str(i) + '! Skipping...'
            continue
        temp_df = error_df.dropna(axis=0, how='any')
        temp_arr = np.array([np.random.laplace(0, temp_df.loc[str(i), 'sig_del']).sum() for j in xrange(num_trials)])
        years_df.loc[i, 'Random error (Fc)'] = temp_arr.std() * crit_t * 12 * 10 **-6 * flux_period * 60
  
    return years_df
    
#-----------------------------------------------------------------------------#

def random_error_generator(error_df):

    print '\nCalculating single realisation of random error...'
    print '-------------------------------------------------'    
    error_df['random_error'] = np.random.laplace(0, error_df.loc[:, 'sig_del'])
    
#-----------------------------------------------------------------------------#


#-----------------------------------------------------------------------------#

def myround(x,base=10):
    return int(base*round(x/base))

# Create two plot axes
def main_plot(Fc_diff,cats_df,linreg_stats_df,num_classes):

    fig = plt.figure(figsize=(12, 6))
    fig.patch.set_facecolor('white')
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 2])
    ax1 = plt.subplot(gs[0])    
    ax2 = plt.subplot(gs[1])
    ax3=ax2.twinx()

    ### Histogram ###    
    
    # Calculate scaling parameter for Laplace (sigma / sqrt(2)) and Gaussian 
    # (sigma) distributions over entire dataset 
    beta=(abs(Fc_diff-Fc_diff.mean())).mean()
    sig=Fc_diff.std()

    # Get edge quantiles and range, then calculate Laplace pdf over range
    x_low=myround(Fc_diff.quantile(0.005))
    x_high=myround(Fc_diff.quantile(0.995))
    x_range=(x_high-x_low)
    x=np.arange(x_low,x_high,1/(x_range*10.))
    pdf_laplace=np.exp(-abs(x/beta))/(2.*beta)
	    	
    # Plot normalised histogram with Laplacian and Gaussian pdfs
    ax1.hist(np.array(Fc_diff), bins=200, range=[x_low,x_high], normed=True,
            color='0.7', edgecolor='none')
    ax1.plot(x,pdf_laplace,color='black', label='Laplacian')
    ax1.plot(x,mlab.normpdf(x,0,sig),color='black', linestyle = '--', 
            label='Gaussian')
    ax1.set_xlabel(r'$\delta\/(\mu mol\/m^{-2} s^{-1}$)',fontsize=22)
    ax1.set_ylabel('$Percent$', fontsize=22)
    ax1.axvline(x=0, color = 'black', linestyle = ':')
    ax1.tick_params(axis = 'x', labelsize = 14)
    ax1.set_yticks([0, 0.1, 0.2, 0.3, 0.4])
    ax1.tick_params(axis = 'y', labelsize = 14)

    ### Line plot ###

    # Create series for regression lines
    influx_x=cats_df['Fc_mean'][cats_df['Fc_mean']<=0]
    influx_y=np.polyval([linreg_stats_df.ix['day'][0],linreg_stats_df.ix['day'][1]],influx_x)
    efflux_x=cats_df['Fc_mean'][cats_df['Fc_mean']>=0]
    efflux_y=np.polyval([linreg_stats_df.ix['noct'][0],linreg_stats_df.ix['noct'][1]],efflux_x)
    
    # Do formatting
    ax2.set_xlim(round(cats_df['Fc_mean'].iloc[0]),math.ceil(cats_df['Fc_mean'].iloc[-1]))
    ax2.set_xlabel(r'$C\/flux\/(\mu molC\/m^{-2} s^{-1}$)',fontsize=22)
    ax2.set_ylabel('$\sigma(\delta)\/(\mu molC\/m^{-2} s^{-1})$',fontsize=22)
    ax2.yaxis.set_ticklabels([])
    ax2.tick_params(axis='y', which='both', left='off', right='off')
    ax3.spines['right'].set_position('zero')
    ax3.set_yticks([1, 2, 3, 4, 5, 6, 7])
    plt.setp(ax3.get_yticklabels()[-1], visible = False)


    # Do plotting
    ax2.plot(cats_df['Fc_mean'], cats_df['sig_del'], 'o', markerfacecolor='0.8',
             markeredgecolor='black', markersize=6)
    ax2.plot(influx_x,influx_y, linestyle=':', color='black')
    ax2.plot(efflux_x,efflux_y, linestyle=':', color='black')

#    ax.set_title('Random error SD binned over flux magnitude quantiles (n='+str(num_classes)+')\n')
    
    # Move axis and relabel
    str_influx=('a = '+str(round(linreg_stats_df.ix['day'][0],2))+
                '\nb = '+str(round(linreg_stats_df.ix['day'][1],2))+
                '\nr = '+str(round(linreg_stats_df.ix['day'][2],2))+
                '\np = '+str(round(linreg_stats_df.ix['day'][3],2))+
                '\nSE = '+str(round(linreg_stats_df.ix['day'][4],2)))
    str_efflux=('a = '+str(round(linreg_stats_df.ix['noct'][0],2))+
                '\nb = '+str(round(linreg_stats_df.ix['noct'][1],2))+
                '\nr = '+str(round(linreg_stats_df.ix['noct'][2],2))+
                '\np = '+str(round(linreg_stats_df.ix['noct'][3],2))+
                '\nSE = '+str(round(linreg_stats_df.ix['noct'][4],2)))
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax2.text(0.055, 0.28, str_influx, transform=ax2.transAxes, fontsize=12,
            verticalalignment='top',bbox=props)
    ax2.text(0.83, 0.28, str_efflux, transform=ax2.transAxes, fontsize=12,
            verticalalignment='top',bbox=props)        

    fig.tight_layout()
    fig.show()
        
#-----------------------------------------------------------------------------#