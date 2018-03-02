#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 12:29:02 2018

@author: ian
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pdb
from scipy import stats

import DataIO as io
reload(io)


class random_error(object):
    
    def __init__(self, dataframe, t_threshold = 3, ws_threshold = 1,
                 k_threshold = 35, configs_dict = False):
        
        if not configs_dict:
            configs_dict = {'flux_name': 'Fc',
                            'mean_flux_name': 'Fc_SOLO',
                            'windspeed_name': 'Ws',
                            'temperature_name': 'Ta',
                            'insolation_name': 'Fsd',
                            'QC_name': 'Fc_QCFlag',
                            'QC_code': 0}
        
        # Get and check the interval
        interval = int(filter(lambda x: x.isdigit(), 
                              pd.infer_freq(dataframe.index)))
        assert interval % 30 == 0
        recs_per_day = 1440 / interval
        self.recs_per_day = recs_per_day
        self.df = dataframe
        self.configs_dict = configs_dict
        self.t_threshold = t_threshold
        self.ws_threshold = ws_threshold
        self.k_threshold = k_threshold

    #------------------------------------------------------------------------------
    # Calculate regression parameters for random error
    #------------------------------------------------------------------------------
    def get_flux_binned_sigma_delta(self):    

        # Set stuff
        num_cats = 60
        nocturnal_threshold = 10
    
        #----------------------------------------------------------------------
        # Internal functions
        #----------------------------------------------------------------------
        # Split day and night
        def bin_series():
            
            convert_names_and_qc
            
            
            
            
            
            
            
            
            
            
            # Function that actually does the binning
            def get_sigmas(df):
                def calc(s):
                    return abs(s - s.mean()).mean() * np.sqrt(2)
                return pd.DataFrame({'sigma_delta': 
                                      map(lambda x: calc(df.loc[df['quantile_label'] == x, 
                                                                'flux_diff']), 
                                          df['quantile_label'].unique().categories),
                                     'mean': 
                                      map(lambda x: df.loc[df['quantile_label'] == x,
                                                           'flux_mean'].mean(),
                                          df['quantile_label'].unique().categories)})
            
            # Separate into night and day 
            noct_df = filter_df.loc[filter_df.Fsd_mean < nocturnal_threshold, 
                                    ['flux_mean', 'flux_diff']]
            day_df = filter_df.loc[filter_df.Fsd_mean > nocturnal_threshold, 
                                   ['flux_mean', 'flux_diff']]
                    
            # Calculate n bins for day and night
            nocturnal_propn = float(len(noct_df)) / len(filter_df)
            num_cats_night = int(round(num_cats * nocturnal_propn))
            num_cats_day = num_cats - num_cats_night
            
            # Do the binned flux mean and sigma_delta calculation for night
            noct_df['quantile_label'] = pd.qcut(noct_df.flux_mean, num_cats_night, 
                                                labels = np.arange(num_cats_night))
            noct_group_df = get_sigmas(noct_df)
            noct_group_df = noct_group_df.loc[noct_group_df['mean'] > 0]
    
            # Do the binned flux mean and sigma_delta calculation for day
            day_df['quantile_label'] = pd.qcut(day_df.flux_mean, num_cats_day, 
                                               labels = np.arange(num_cats_day))
            day_group_df = get_sigmas(day_df)
            day_group_df = day_group_df.loc[day_group_df['mean'] < 0]

            return noct_group_df, day_group_df
        
        #----------------------------------------------------------------------
    # Convert external to internal names and drop filled flux time series 
    # values
    def convert_names_and_QC(self):
        
        new_names_dict = {'flux_name': 'flux',
                          'mean_flux_name': 'flux_mean',
                          'QC_name': 'QC',
                          'windspeed_name': 'Ws',
                          'temperature_name': 'Ta',
                          'insolation_name': 'Fsd'}

        internal_dict = self.configs_dict.copy()
        try:
            QC_code = internal_dict.pop('QC_code')
            assert 'QC_name' in internal_dict.keys()
        except (KeyError, AssertionError):
            QC_code = None
            new_names_dict.pop('QC_name')
        old_names = [internal_dict[name] for name in 
                     sorted(internal_dict.keys())]
        new_names = [new_names_dict[name] for name in 
                     sorted(new_names_dict.keys())]
        sub_df = self.df[old_names].copy()
        sub_df.columns = new_names
        if not QC_code is None:
            sub_df.loc[sub_df.QC != QC_code, 'flux'] = np.nan
            sub_df.drop('QC', axis = 1, inplace = True)
        return sub_df
        
        #--------------------------------------------------------------------------    
        # Do differencing
        def difference_time_series():
            diff_df = pd.DataFrame(index = work_df.index)
            for var in ['flux', 'Ta', 'Fsd', 'Ws']:
                var_name = var + '_diff'
                temp = work_df[var] - work_df[var].shift(recs_per_day) 
                diff_df[var_name] = temp if var == 'flux' else abs(temp)
            diff_df['flux_mean'] = (work_df['flux_mean'] + 
                                    work_df['flux_mean'].shift(recs_per_day)) / 2
            diff_df['Fsd_mean'] = (work_df['Fsd'] + 
                                   work_df['Fsd'].shift(recs_per_day)) / 2
            return diff_df
        
        #--------------------------------------------------------------------------
        # Filter on criteria
        def filter_time_series():
            bool_s = ((diff_df['Ws_diff'] < ws_threshold) & 
                      (diff_df['Ta_diff'] < t_threshold) & 
                      (diff_df['Fsd_diff'] < k_threshold))
            return pd.DataFrame({var: diff_df[var][bool_s] for var in 
                                 ['flux_diff', 'flux_mean', 'Fsd_mean']}).dropna()
        #--------------------------------------------------------------------------
        # Do plotting
        def plot_data():
            
            fig, ax1 = plt.subplots(1, 1, figsize = (14, 8))
            fig.patch.set_facecolor('white')
            ax1.xaxis.set_ticks_position('bottom')
            ax1.set_xlim([round(combined_df['mean'].min() * 1.05), 
                          round(combined_df['mean'].max() * 1.05)])
            ax1.set_ylim([0, round(combined_df['sigma_delta'].max() * 1.05)])
            ax1.yaxis.set_ticks_position('left')
            ax1.spines['right'].set_visible(False)
            ax1.spines['top'].set_visible(False)
            ax1.tick_params(axis = 'y', labelsize = 14)
            ax1.tick_params(axis = 'x', labelsize = 14)
            
            ax1.set_xlabel('$flux\/(\mu mol\/CO_2\/m^{-2}\/s^{-1})$', fontsize = 18)
            ax1.set_ylabel('$\sigma[\delta]\/(\mu mol\/CO_2\/m^{-2}\/s^{-1})$', fontsize = 18)
            ax2 = ax1.twinx()
            ax2.spines['right'].set_position('zero')
            ax2.spines['left'].set_visible(False)
            ax2.spines['top'].set_visible(False)
            ax2.set_ylim(ax1.get_ylim())
            ax2.tick_params(axis = 'y', labelsize = 14)
            plt.setp(ax2.get_yticklabels()[0], visible = False)
            
            ax1.plot(combined_df['mean'], combined_df['sigma_delta'], 
                     'o', mfc = 'grey')
            
            day_x = np.linspace(ax1.get_xlim()[0], 0, 5)
            day_y = day_x * stats_dict['day'].slope + stats_dict['day'].intercept
            night_x = np.linspace(0, ax1.get_xlim()[-1], 5)
            night_y = night_x * stats_dict['night'].slope + stats_dict['night'].intercept
            ax1.plot(day_x, day_y, color = 'grey')
            ax1.plot(night_x, night_y, color = 'grey')
            
            return fig
        #--------------------------------------------------------------------------
        # Main routine
        #--------------------------------------------------------------------------
        

        
        # Convert names
        work_df = convert_names_and_QC()
        
        # Do the differencing
        diff_df = difference_time_series()

        # Filter
        filter_df = filter_time_series()
            
        # Do calculations for nocturnal and daytime
        noct_group_df, day_group_df = bin_series()
        
        return pd.concat([day_group_df, noct_group_df]).reset_index(drop = True)
#        
#        # Do the stats
#        stats_dict = do_stats()
#        
#        # Organise the outputs
#        output_dict = {'Stats': stats_dict}
#        combined_df = 
#        if return_data: output_dict['Data'] = combined_df
#        if return_plot: output_dict['Figure'] = plot_data()
#        
#        return output_dict
    #------------------------------------------------------------------------------
    
    #--------------------------------------------------------------------------
    # Calculate basic regression statistics
    def regression_statistics(binned_df):
            return {'night': stats.linregress(noct_group_df['mean'], 
                                              noct_group_df['sigma_delta']),
                    'day': stats.linregress(day_group_df['mean'], 
                                            day_group_df['sigma_delta'])}

#------------------------------------------------------------------------------
# Estimate sigma_delta for a time series using regression coefficients
#------------------------------------------------------------------------------
def estimate_sigma_delta(flux_series, stats_dict):
    """
    Calculates sigma_delta value for each member of a time series
    
    Args:
        * flux_series (array-like): the series for which sigma-delta is to be \
        calculated
        * stats_dict (dictionary): regression parameters (slope and intercept) \
        generated from *regress_sigma_delta()*

    Returns:
        * ret1 (array-like): time series of sigma_delta for each valid value in \
        the originally passed array-like
    """
    sig_del_series = np.where(flux_series > 0, flux_series * 
                                               stats_dict['night'].slope + 
                                               stats_dict['night'].intercept,
                                               flux_series * 
                                               stats_dict['day'].slope + 
                                               stats_dict['day'].intercept)
    if any(sig_del_series < 0):
        n_below = len(sig_del_series[sig_del_series < 0])
        print ('Warning: approximately {0} estimates of sigma_delta have value ' 
               'less than 0 - setting to mean of all other values'
               .format(str(n_below)))
        sig_del_series = np.where(sig_del_series > 0, 
                                  sig_del_series, 
                                  sig_del_series[sig_del_series > 0].mean())
    return sig_del_series
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Generate scaled noise realisation for time series
#------------------------------------------------------------------------------
def estimate_random_error(sigma_delta_series):
    """
    Generates a random error realisation for a single point
    Pass the following arguments: 1) numpy array or pandas series 
                                     containing sigma_delta estimate for NEE
    Returns a series of random error estimate drawn from the Laplace 
    distribution with sigma_delta as the scaling parameter
    (location parameter is 0)
    """
    return np.random.laplace(0, sigma_delta_series / np.sqrt(2))
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Propagate random error
#------------------------------------------------------------------------------
def propagate_random_error(sigma_delta_series, n_trials, scaling_coefficient = 1):
    '''
    Function to run Monte Carlo-style trials to assess uncertainty due to 
    random error.
    
    Args:
        * sigma_delta_series (array-like): sigma_delta estimates for each valid \
        data point in the flux series.
        * n_trials (int):  number of trials over which to compound the sum.
        * scaling_coefficient (int or float): scales summed value to required \
        units.
    
    Returns:
        * float: scaled estimate of 2-sigma bounds of all random error\
        trial sums.
    '''
    # Calculate critical t-statistic for p = 0.095
    crit_t = stats.t.isf(0.025, n_trials)    

    results_list = []
    for this_trial in xrange(n_trials):
        results_list.append(sum(estimate_random_error(sigma_delta_series) *
                                scaling_coefficient))
    return float(pd.DataFrame(results_list).std() * crit_t)
#------------------------------------------------------------------------------


    '''
    Calculates regression for sigma_delta on flux magnitude 
    
    Args:
        * df (pandas dataframe): dataframe containing the required data, with \
        names defined in the configuration dictionary (see below).
        * config_dict (dictionary): dictionary for specifying naming convention \
        for dataset passed to the function; keys must use the default names \
        (specified below) expected by the script, and the values specify the \
        name the relevant variable takes in the dataset; default names are as \
        follows:\n
            - flux_name: the turbulent flux for which to calculate random \
            error 
            - mean_flux_name: the flux series to use for calculating the \
            bin average for the flux when estimating sigma_delta \
            (since the turbulent flux already contains random error, a \
            model series is generally recommended for this purpose)
            - windspeed_name
            - temperature_name
            - insolation_name
            - QC_name: (optional) if passed, used as a filter variable \
            for the flux variable (if QC_code is not present, this is \
            ignored, and no warning or error is raised)
            - QC_code: (optional, int) if passed, all flux records that \
            coincide with the occurrence of this code in the QC variable \
            are retained, and all others set to NaN
    Kwargs:
        * return_data (bool): if True attaches bin averaged flux and \
        sigma_delta estimates to the returned results dictionary
        * return_plot (bool): if True attaches plot of sigma_delta as a \
        function of flux magnitude to the returned results dictionary
        * t_threshold (int): user-set temperature difference threshold \
        default = 3superscriptoC)
        * ws_threshold (int): user-set wind speed difference threshold \
        (default = 1m s\^{-1})
        * k_threshold (int): user set insolation difference threshold \
        (default = 35Wm-2)
    '''
    


##------------------------------------------------------------------------------
## Main function for stringing the above together
##------------------------------------------------------------------------------
#def main(df, scaling_coefficient, config_dict = None, n_trials = 10**4):
#    '''
#    Function that combines the individual functions of the random_error module
#    together for convenience.
#    
#    Args:
#        * df (pandas dataframe): dataframe containing the required data, with \
#        names defined in the configuration dictionary (see below).
#        * scaling_coefficient (int or float): converts the flux data (an \
#        instantaneous average quantity) to the appropriate units for \
#        integrating over the time series length.
#        
#    Kwargs:
#        * config_dict (dictionary): sictionary for specifying naming convention \
#        for dataset passed to the function; keys must use the default names \
#        (specified below) expected by the script, and the values specify the \
#        name the relevant variable takes in the dataset; default names are as \
#        follows:\n
#            - flux_name: the turbulent flux for which to calculate random \
#            error 
#            - mean_flux_name: the flux series to use for calculating the \
#            bin average for the flux when estimating sigma_delta \
#            (since the turbulent flux already contains random error, a \
#            model series is generally recommended for this purpose)
#            - windspeed_name
#            - temperature_name
#            - insolation_name
#            - QC_name: (optional) if passed, used as a filter variable \
#            for the flux variable (if QC_code is not present, this is \
#            ignored, and no warning or error is raised)
#            - QC_code: (optional, int) if passed, all flux records that \
#            coincide with the occurrence of this code in the QC variable \
#            are retained, and all others set to NaN
#        * n_trials (int): the number of trials over which to compound to \
#        calculate the uncertainty due to random error.
#    
#    Returns:
#        * float: 2 sigma of the population (n = n_trials) of random error sums
#    '''
#    
#    if config_dict is None:
#        
#        config_dict = {'flux_name': 'Fc',
#                       'mean_flux_name': 'Fc_SOLO',
#                       'windspeed_name': 'Ws',
#                       'temperature_name': 'Ta',
#                       'insolation_name': 'Fsd',
#                       'QC_name': 'Fc_QCFlag',
#                       'QC_code': 0}
#        
#    results = regress_sigma_delta(df, config_dict)
#
#    sigma_delta_series = estimate_sigma_delta(df[config_dict ['mean_flux_name']], 
#                                                 results['Stats'])
#
#    summed_error = propagate_random_error(sigma_delta_series, 10000, 
#                                          scaling_coefficient)
#    
#    return summed_error
##------------------------------------------------------------------------------
#
#
#path = '/home/ian/OzFlux/Sites/GatumPasture/Data/Processed/All/GatumPasture_L5.nc'
#
#df = io.OzFluxQCnc_to_data_structure(path, output_structure = 'pandas')
#
#scaling_coefficient = 30 * 60 * 12 * 10**-6
#
#error = main(df, scaling_coefficient)