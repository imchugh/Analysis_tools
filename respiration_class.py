#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 08:21:07 2017

@author: ian
"""
import numpy as np
from scipy. optimize import curve_fit
import matplotlib.pyplot as plt

class MyClass(object):

    def __init__(self, T, ER, sws = None):
        self.T = T
        self.ER = ER
        self.sws = sws
    
    def get_respiration(self, temp, rb, Eo):
        return rb  * np.exp(Eo * (1 / (10 + 46.02) - 1 / (temp + 46.02)))

    def get_respiration_sm(self, temp, sws, rb, Eo, theta_1, theta_2):
        return (rb  * np.exp(Eo * (1 / (10 + 46.02) - 1 / (temp + 46.02))) * 
                (1 / (1 + np.exp(theta_1 - theta_2 * sws))))

    def get_fit(self, Eo = None, theta_1 = None, theta_2 = None):
        
        p0_list = [1, 100, 1, 10]
        
        if not theta_1 == None and theta_2 == None:
            if self.sws == None:
                raise RuntimeError('Series sws must be passed to class '
                                   'instance if theta parameters are passed '
                                   'to fitting function... exiting!')
        try:
            if Eo == None:
                params, cov = curve_fit(self.get_respiration, self.T, self.ER, 
                                        p0 = [1, 100])
            else:
                params, cov = curve_fit(lambda x, a: 
                                        self.get_respiration(x, a, Eo), 
                                        self.T, self.ER, 
                                        p0 = [1])
            error_state = 0
        except RuntimeError:
            params = [np.nan, np.nan]
            cov = None
            error_state = 1
        results_d = {'parameters': params,
                     'covariance_matrix': cov,
                     'error_state': error_state}
        return results_d

    def plot_respiration(self, title_str):
        
        results_dict = self.get_fit()
         
        x = self.T
        y1 = self.ER
        y2 = self.get_respiration(self.T, results_dict['parameters'][0], 
                                  results_dict['parameters'][1])
        
        # Plot
        fig = plt.figure(figsize = (12,8))
        fig.patch.set_facecolor('white')
        ax = plt.gca()
        ax.plot(x, y1, 'o' , markerfacecolor = 'none',
                 markeredgecolor = 'black', label = 'NEE_obs', color = 'black')
        ax.plot(x, y2, linestyle = ':', color = 'black', 
                 label = 'NEE_est')
        ax.set_title(title_str)
        ax.set_xlabel('$Temperature\/(^oC)$', fontsize = 18)
        ax.set_ylabel('$NEE\/(\mu mol C\/m^{-2}\/s^{-1}$)', fontsize = 18)
        ax.axhline(y = 0, color = 'black')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        
        return fig
        
Eo = 200
rb = 2.5
theta_1 = 3
theta_2 = 30

# Make some fake data
temp = np.linspace(0, 30, 100)
vwc = np.linspace(0.5, 0.1, 100)
vwc_func = 1 / (1 + np.exp(theta_1 - theta_2 * vwc))
est_resp = (rb  * np.exp(Eo * (1 / (10 + 46.02) - 1 / (temp + 46.02))) +
            np.random.randn(100)) #* vwc_func

# Instantiate class
mc = MyClass(temp, est_resp)

# Get the results dictionary
d = mc.get_fit()

# Get a respiration series based on parameters
#a = mc.get_respiration(temp, d['parameters'][0], d['parameters'][1])