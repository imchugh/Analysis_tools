#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 08:21:07 2017

@author: ian
"""
import numpy as np
import pdb
from scipy. optimize import curve_fit
import matplotlib.pyplot as plt

class MyClass(object):

    def __init__(self, T, ER, sws = None):
        self.T = T
        self.ER = ER
        self.sws = sws
        
        if sws == None:
            self.drivers = np.reshape(T, [len(T), 1])
        else:
            self.drivers = np.column_stack([T, sws])            
       
    def get_respiration(self, drivers, *params):
        T_response = params[0]  * np.exp(params[1] * (1 / (10 + 46.02) - 
                                         1 / (drivers[:, 0] + 46.02)))
        if drivers.shape[1] == 1:
            return T_response
        else:
            return T_response * (1 / (1 + np.exp(params[2] - params[3] *
                                                 drivers[:, 1])))

    def get_fit(self, rb = None, Eo = None, theta_1 = None, theta_2 = None):
        
#       Create a binary word from parameter arguments and generate a base-10
#       ID and specify a set of starting values for the parameters to be 
#       fitted
        init_est_array = np.array([1, 100, 1, 10])
        param_list = [rb, Eo, theta_1, theta_2]
        if self.sws == None:
            init_est_array = init_est_array[:-2]
            param_list = param_list[:-2]
        else:
            if not theta_1 == None and theta_2 == None:
                raise RuntimeError('Theta parameters can only be passed to '
                                   'fitting function of a class instance '
                                   'containing Series sws... exiting!')
        
        bool_array = np.array([i == None for i in param_list])
        p0_list = init_est_array[bool_array]

        if Eo == None:
            func = lambda x, a, b: self.get_respiration(x, a, b)
        else:
            func = lambda x, a: self.get_respiration(x, a, Eo)

        try:
            params, cov = curve_fit(func, self.drivers, self.ER, p0 = p0_list)
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