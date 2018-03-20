#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 13:53:16 2018

@author: ian
"""

import datetime as dt
from lmfit import Model
import numpy as np
import pandas as pd
import pdb

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
    def __init__(self, dataframe, names_dict = None, weighting = 'air'):
        
        df = dataframe.copy()
        if not names_dict: names_dict = self.get_default_external_names()
        self.df = self._rewrite_dataframe_for_internal(df, names_dict, 
                                                       weighting)
        
        interval = int(filter(lambda x: x.isdigit(), 
                              pd.infer_freq(dataframe.index)))
        assert interval % 30 == 0
        self.interval = interval
#------------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    # Methods
    #--------------------------------------------------------------------------
    
    #--------------------------------------------------------------------------
    def get_subset(self, date, size):
        
        noct_threshold = 10
        ref_date = date + dt.timedelta(0.5)
        date_tuple = (ref_date - dt.timedelta(size / 2.0 - 
                                              self.interval / 1440.0),
                      ref_date + dt.timedelta(size / 2.0))
        sub_df = self.df.loc[date_tuple[0]: date_tuple[1], 
                             ['NEE', 'Fsd', 'TC']].dropna()
        return sub_df.loc[sub_df.Fsd < noct_threshold].drop('Fsd', axis = 1)
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
        date_list = self.make_date_iterator(size = window_size, 
                                            step = window_step)
        for date in date_list:
            df = self.get_subset(date, size = window_size)
            if not len(df) > 6: continue
            model = Model(_LT_Eo_long, independent_vars = ['t_series'])
            params = model.make_params(rb = 1, Eo = 100)
            result = model.fit(df.NEE,
                               t_series = df.TC,
                               params = params)
            if not 50 < result.params['Eo'].value < 400: continue
            se = (result.conf_interval()['Eo'][4][1] - 
                  result.conf_interval()['Eo'][2][1]) / 2
            if se > result.params['Eo'].value / 2.0: continue
            Eo_list.append([result.params['Eo'].value, se])
        print ('Found {} valid estimates of Eo'.format(str(len(Eo_list))))
        if len(Eo_list) == 0: raise RuntimeError
        Eo_array = np.array(Eo_list)
        Eo = ((Eo_array[:, 0] / (Eo_array[:, 1])).sum() / 
              (1 / Eo_array[:, 1]).sum())
        if not 50 < Eo < 400: raise RuntimeError
        return Eo
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    def estimate_rb(self, window_size = 4, window_step = 4, Eo = None):
        
        out_df = pd.DataFrame(index = self.make_date_iterator
                              (size = window_size, 
                               step = window_step))
        if not Eo:
            try:
                Eo = self.estimate_Eo()
            except RuntimeError:
                print 'Could not find any valid values of Eo'
                return
        out_df['Eo'] = Eo
        out_df['rb'] = np.nan
        for date in out_df.index:
            df = self.get_subset(date, size = window_size)
            if not len(df) > 0: continue
            test = _LT_Eo_long
            model = Model(test, independent_vars = ['t_series'])
            params = model.make_params(rb = 1, Eo = Eo)
            params['Eo'].vary = False
            result = model.fit(df.NEE,
                               t_series = df.TC, 
                               params = params)
            if result.params['rb'].value < 0: continue
            out_df.loc[date, 'rb'] = result.params['rb'].value
        out_df = out_df.resample('D').interpolate()
        out_df = out_df.reindex(np.unique(self.df.index.date))
        out_df.fillna(method = 'bfill', inplace = True)
        out_df.fillna(method = 'ffill', inplace = True)
        return out_df
    #--------------------------------------------------------------------------            
    
    #--------------------------------------------------------------------------
    def estimate_er(self, params_df = False):
        
        if not isinstance(params_df, pd.core.frame.DataFrame):
            params_df = self.estimate_rb()
        resp_series = pd.Series()
        for date in params_df.index:
            Eo = params_df.loc[date, 'Eo']
            rb = params_df.loc[date, 'rb']
            str_date = dt.datetime.strftime(date, '%Y-%m-%d')
            t_series = self.df.loc[str_date, 'TC']
            resp_series = resp_series.append(_LT_Eo_long(t_series = t_series, 
                                                         Eo = Eo, rb = rb))
        return resp_series
    #--------------------------------------------------------------------------
    
    #--------------------------------------------------------------------------
    def get_default_external_names(self):

        return {'Cflux': 'Fc',
                'air_temperature': 'Ta',
                'soil_temperature': 'Ts',
                'insolation': 'Fsd'}
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    def _rewrite_dataframe_for_internal(self, df, external_names, weighting):

        internal_names = {'Cflux': 'NEE',
                          'air_temperature': 'Ta',
                          'soil_temperature': 'Ts',
                          'insolation': 'Fsd'}
        
        swap_dict = {external_names[key]: internal_names[key] 
                     for key in internal_names.keys()}
        sub_df = df[swap_dict.keys()]
        sub_df.columns = swap_dict.values()
        
        if weighting == 'air':
            s = sub_df[internal_names['air_temperature']].copy()
        elif weighting == 'soil':
            s = sub_df[internal_names['soil_temperature']].copy()
        elif isinstance(weighting, (int, float)):
            s = ((sub_df['Ta'] * weighting + sub_df['Ts']) / (weighting + 1))
        s.name = 'TC'
        return sub_df.join(s)
    #--------------------------------------------------------------------------
    
#------------------------------------------------------------------------------
def _LT_Eo_long(t_series, rb, Eo):

    return rb  * np.exp(Eo * (1 / (10 + 46.02) - 1 / (t_series + 46.02)))
#------------------------------------------------------------------------------