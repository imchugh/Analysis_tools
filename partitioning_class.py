#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 13:53:16 2018

@author: ian
"""

import datetime as dt
from lmfit import Model
import numpy as np
import operator
import pandas as pd
import pdb

import utils

#------------------------------------------------------------------------------
# Init
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
class respiration(object):
    """
    Class for fitting of respiration parameters and estimation of respiration
    
    Args:
        * dataframe (pd.dataframe): containing a minimum of temperature, solar 
          radiation and CO2 flux.
    Kwargs:
        * names_dict (dict): maps the variable names used in the dataset to 
          common names (keys must be 'air_temperature', 'soil_temperature', 
          'insolation', 'Cflux'); if None, defaults to the internal 
          specification, which works for PyFluxPro.
        * weighting (str, int or float): if str, must be either 'air' or 'soil',
          which determines which temperature series is used for the fit; if 
          int or float is supplied, the number is used as a ratio for weighting
          the air and soil series - note that the ratio is air: soil, such that
          e.g. choice of 3 would cause weighting of 3:1 in favour of air 
          temperature, or e.g. float(1/3) would result in the reverse.
    """
    def __init__(self, dataframe, names_dict = None, weighting = 'air',
                 fit_daytime_rb = False, noct_threshold = 10,
                 convert_to_photons = True):

        interval = int(filter(lambda x: x.isdigit(), 
                              pd.infer_freq(dataframe.index)))
        assert interval % 30 == 0
        assert isinstance(fit_daytime_rb, bool)
        self.interval = interval 
        if not names_dict: 
            self.external_names = self._define_default_external_names()
        else:
            self.external_names = names_dict
        self.internal_names = self._define_default_internal_names()

        self.df = self._make_formatted_df(dataframe, weighting)
        self.fit_daytime_rb = fit_daytime_rb
        self.noct_threshold = noct_threshold
        self.convert_to_photons = convert_to_photons
#------------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    # Methods
    #--------------------------------------------------------------------------
    
    #--------------------------------------------------------------------------
    def get_subset(self, date, size, mode):
        
        ops = {">": operator.gt, "<": operator.lt}
        ref_date = date + dt.timedelta(0.5)
        date_tuple = (ref_date - dt.timedelta(size / 2.0 - 
                                              self.interval / 1440.0),
                      ref_date + dt.timedelta(size / 2.0))
        sub_df = self.df.loc[date_tuple[0]: date_tuple[1], 
                             ['NEE', 'Fsd', 'TC', 'VPD']].dropna()
        return sub_df[ops[mode](sub_df.Fsd, self.noct_threshold)]
    #--------------------------------------------------------------------------
    
    #--------------------------------------------------------------------------
    def make_date_iterator(self, size, step):
        
        start_date = (self.df.index[0].to_pydatetime().date() + 
                      dt.timedelta(size / 2))
        end_date = self.df.index[-1].to_pydatetime().date()
        return pd.date_range(start_date, end_date, 
                             freq = '{}D'.format(str(step)))
    #--------------------------------------------------------------------------
    
    #--------------------------------------------------------------------------
    def estimate_Eo(self, window_size = 15, window_step = 5):
        
        Eo_list = []
        for date in self.make_date_iterator(window_size, window_step):
            df = self.get_subset(date, size = window_size, mode = '<')
            if not len(df) > 6: continue
            if not df.TC.max() - df.TC.min() >= 5: continue
            f = _Lloyd_and_Taylor
            model = Model(f, independent_vars = ['t_series'])
            params = model.make_params(rb = 1, 
                                       Eo = self.get_prior_parameter_estimates()['Eo'])
            result = model.fit(df.NEE,
                               t_series = df.TC,
                               params = params)
            if not 50 < result.params['Eo'].value < 400: continue
            se = (result.conf_interval()['Eo'][4][1] - 
                  result.conf_interval()['Eo'][2][1]) / 2
            if se > result.params['Eo'].value / 2.0: continue
            Eo_list.append([result.params['Eo'].value, se])
        if len(Eo_list) == 0: raise RuntimeError('Could not find any valid '
                                                 'estimates of Eo! Exiting...')
        print ('Found {} valid estimates of Eo'.format(str(len(Eo_list))))
        Eo_array = np.array(Eo_list)
        Eo = ((Eo_array[:, 0] / (Eo_array[:, 1])).sum() / 
              (1 / Eo_array[:, 1]).sum())
        if not 50 < Eo < 400: raise RuntimeError('Eo value {} outside '
                                                 'acceptable parameter range '
                                                 '(50-400)! Exiting...'
                                                 .format(str(round(Eo, 2))))
        return Eo
    #--------------------------------------------------------------------------
        
    #--------------------------------------------------------------------------
    def get_parameters(self, mode, Eo = None, window_size = 4, window_step = 4):
        
        priors_dict = self.get_prior_parameter_estimates()
        func = self._get_func()[mode]
        if not Eo: Eo = self.estimate_Eo()
        result_list, date_list = [], []
        print 'Processing the following dates ({} mode): '.format(mode)
        for date in self.make_date_iterator(window_size, window_step):
            print date.date(),
            try:
                result_list.append(func(date, Eo, window_size, priors_dict))
                date_list.append(date)
                print
            except RuntimeError, e:
                print '- {}'.format(e)
                continue
        out_df = pd.DataFrame(result_list, index = date_list)
        out_df = out_df.resample('D').interpolate()
        out_df = out_df.reindex(np.unique(self.df.index.date))
        out_df.fillna(method = 'bfill', inplace = True)
        out_df.fillna(method = 'ffill', inplace = True)
        return out_df
    #--------------------------------------------------------------------------
    
    #--------------------------------------------------------------------------
    def _get_func(self):
        
        return {'night': self.nocturnal_params, 'day': self.daytime_params}
    #--------------------------------------------------------------------------
    
    #--------------------------------------------------------------------------
    def nocturnal_params(self, date, Eo, window_size, priors_dict):
        
        df = self.get_subset(date, size = window_size, mode = '<')
        if not len(df) > 2: raise RuntimeError
        f = _Lloyd_and_Taylor
        model = Model(f, independent_vars = ['t_series'])
        params = model.make_params(rb = priors_dict['rb'], 
                                   Eo = Eo)
        params['Eo'].vary = False
        result = model.fit(df.NEE,
                           t_series = df.TC, 
                           params = params)
        if result.params['rb'].value < 0: raise RuntimeError
        return result.best_values
    #--------------------------------------------------------------------------    

    #--------------------------------------------------------------------------
    def daytime_params(self, date, Eo, window_size, priors_dict):
        
        def model_fit(these_params):
            return model.fit(df.NEE,
                             par_series = df.Fsd, vpd_series = df.Fsd, 
                             t_series = df.TC, params = these_params)
            
        if self.fit_daytime_rb:
            rb_prior = priors_dict['rb']
        else:
            rb_prior = self.nocturnal_params(date, Eo, window_size, 
                                             priors_dict)['rb']
        beta_prior = priors_dict['beta']
        df = self.get_subset(date, size = window_size, mode = '>')
        try:
            if not len(df) > 4: 
                raise RuntimeError('insufficient data for fit')
            f = _NEE_model
            model = Model(f, independent_vars = ['par_series', 'vpd_series',
                                                 't_series'])
            params = model.make_params(rb = rb_prior, Eo = Eo,
                                       alpha = priors_dict['alpha'],
                                       beta = beta_prior, 
                                       k = priors_dict['k'])
            rmse_list, params_list = [], []
            for this_beta in [beta_prior, beta_prior / 2, beta_prior * 2]:
                params['beta'].value = this_beta
                params['Eo'].vary = False
                params['rb'].vary = self.fit_daytime_rb
                result = model_fit(these_params = params)
                if result.params['rb'] < 0: 
                    raise RuntimeError('rb parameter out of range')
                if not 0 <= result.params['k'].value <= 10:
                    params['k'].value = priors_dict['k']
                    params['k'].vary = False
                    result = model_fit(these_params = params)
                if not -0.22 <= result.params['alpha'].value <= 0:
                    params['alpha'].value = priors_dict['alpha']
                    params['alpha'].vary = False
                    result = model_fit(these_params = params)
                if not -100 <= result.params['beta'].value <= 0: 
                    raise RuntimeError('beta parameter out of range')
                rmse_list.append(np.sqrt(((df.NEE - result.best_fit)**2).sum()))
                params_list.append(result.best_values)
            idx = rmse_list.index(min(rmse_list))
            priors_dict['alpha'] = result.best_values['alpha']
            return params_list[idx]
        except RuntimeError, e:
            priors_dict['alpha'] = 0
            raise RuntimeError(e)
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    def estimate_er(self, params_df = False):
        
        if not isinstance(params_df, pd.core.frame.DataFrame):
            params_df = self.get_parameters(mode = 'night')
        resp_series = pd.Series()
        for date in params_df.index:
            params = params_df.loc[date]
            str_date = dt.datetime.strftime(date, '%Y-%m-%d')
            data = self.df.loc[str_date, 'TC']
            resp_series = resp_series.append(_Lloyd_and_Taylor
                                             (t_series = data, 
                                              Eo = params.Eo, rb = params.rb))
        return resp_series
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    def estimate_gpp(self, params_df = False):
        
        if not isinstance(params_df, pd.core.frame.DataFrame):
            params_df = self.get_parameters(mode = 'day')
        gpp_series = pd.Series()
        for date in params_df.index:
            params = params_df.loc[date]
            str_date = dt.datetime.strftime(date, '%Y-%m-%d')
            data = self.df.loc[str_date, ['TC', 'Fsd', 'VPD']]
            gpp_series = gpp_series.append(_rectangular_hyperbola
                                           (par_series = data.Fsd,
                                            vpd_series = data.VPD,
                                            alpha = params.alpha,
                                            beta = params.beta,
                                            k = params.k))
        return gpp_series
    #--------------------------------------------------------------------------
    
    #--------------------------------------------------------------------------
    def estimate_nee(self, params_df = False):
        return self.estimate_gpp(params_df) + self.estimate_er(params_df)
    #--------------------------------------------------------------------------
    
    #--------------------------------------------------------------------------
    def get_prior_parameter_estimates(self):
        
        return {'rb': self.df.loc[self.df.Fsd < 10, 'NEE'].mean(),
                'Eo': 100,
                'alpha': -0.01,
                'beta': (self.df.loc[self.df.Fsd > 10, 
                                     'NEE'].quantile(0.03)-
                          self.df.loc[self.df.Fsd > 10, 
                                      'NEE'].quantile(0.97)),
                'k': 0}
    #--------------------------------------------------------------------------    
    
    #--------------------------------------------------------------------------
    def _define_default_internal_names(self):

        return {'Cflux': 'NEE',
                'air_temperature': 'Ta',
                'soil_temperature': 'Ts',
                'insolation': 'Fsd',
                'vapour_pressure_deficit': 'VPD'}
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    def _define_default_external_names(self):

        return {'Cflux': 'Fc',
                'air_temperature': 'Ta',
                'soil_temperature': 'Ts',
                'insolation': 'Fsd',
                'vapour_pressure_deficit': 'VPD'}
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    def _make_formatted_df(self, df, weighting):
        
        sub_df = utils.rename_df(df, self.external_names,
                                 self.internal_names)
        if weighting == 'air':
            s = sub_df['Ta'].copy()
        elif weighting == 'soil':
            s = sub_df['Ts'].copy()
        elif isinstance(weighting, (int, float)):
            s = ((sub_df['Ta'] * weighting + sub_df['Ts']) / (weighting + 1))
        s.name = 'TC'
        return sub_df.join(s)
    #--------------------------------------------------------------------------
    
#--------------------------------------------------------------------------
def _Lloyd_and_Taylor(t_series, rb, Eo):

    return rb  * np.exp(Eo * (1 / (10 + 46.02) - 1 / (t_series + 46.02)))
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def _rectangular_hyperbola(par_series, vpd_series, alpha, beta, k):
    
    beta_VPD = beta * np.exp(-k * (vpd_series - 1))
    index = vpd_series <= 1
    beta_VPD[index] = beta
    GPP = (alpha * par_series) / (1 - (par_series / 2000) + 
           (alpha * par_series / beta_VPD))    
    index = par_series < 5
    GPP[index] = 0
    return GPP
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------    
def _NEE_model(par_series, vpd_series, t_series, rb, Eo, alpha, beta, k):
    
    return (_rectangular_hyperbola(par_series, vpd_series, alpha, beta, k) + 
            _Lloyd_and_Taylor(t_series, rb, Eo))
#------------------------------------------------------------------------------