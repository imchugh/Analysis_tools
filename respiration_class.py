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

class ER(object):

    def __init__(self, T, ER, sws = None):
        self.T = T
        self.ER = ER
        self.sws = sws
        
        if sws == None:
            self.drivers = np.reshape(T, [len(T), 1])
        else:
            self.drivers = np.column_stack([T, sws])            
       
    def get_ER(self, drivers, *params):
        T_response = params[0]  * np.exp(params[1] * (1 / (10 + 46.02) - 
                                         1 / (drivers[:, 0] + 46.02)))
        if drivers.shape[1] == 1:
            return T_response
        else:
            return T_response * (1 / (1 + np.exp(params[2] - params[3] *
                                                 drivers[:, 1])))

    def get_fit(self, rb = None, Eo = None, theta_1 = None, theta_2 = None):
        
        # Check if sws is a valid attribute of class instance - if None, set
        # truncating index accordingly; if theta parameters have been passed,
        # warn the user and then ignore
        if self.sws == None:
            index = -2
            if not theta_1 == None and theta_2 == None:
                print ('Theta parameters can only be passed to '
                       'fitting function of a class instance '
                       'containing Series sws... ignoring theta parameters!')
        else:
            index = None

        # Set up lists for string construction to be passed into function 
        # string for evaluation
        passed_list = ['rb', 'Eo', 'theta_1', 'theta_2'][:index]
        fitted_list = ['a', 'b', 'c', 'd'][:index]
        init_est_list = [1, 100, 1, 10][:index]
        boolean_list = [i == None for i in [rb, Eo, theta_1, theta_2]][:index]
        
        # Make strings including function string
        sub_fitted_str = ','.join([fitted_list[i] for i, d 
                                   in enumerate(boolean_list) if d])
        all_str = ','.join([fitted_list[i] if d else passed_list[i] for i, d  
                            in enumerate(boolean_list)])
        func_str = ('lambda x, {0}, self = self: self.get_ER(x, {1})'
                    .format(sub_fitted_str, all_str))

        # Make list of prior estimates for curve_fit
        p0_list = [init_est_list[i] for i, d in enumerate(boolean_list) if d]

        # Make function by evaluating string (find another way before this goes
        # public!!!)      
        func = eval(func_str)

        # Now fit
        try:
            params, cov = curve_fit(func, self.drivers, self.ER, p0 = p0_list)
            error_state = 0
        except RuntimeError:
            params = [np.nan] * len(sub_fitted_str)
            cov = None
            error_state = 1
        
        results_d = {'parameters': params,
                     'covariance_matrix': cov,
                     'error_state': error_state}
        
        return results_d

    def get_series(self):
        
        results_dict = self.get_fit()
        return self.get_ER(self.drivers, results_dict['parameters'][0], 
                                             results_dict['parameters'][1])
        
    def plot_respiration(self, title_str):
        
        if not self.sws == None:
            
            print 'Watch this space! Going to make a 3D contour plot!'
            return
        
        x = self.drivers[:, 0]
        y1 = self.ER
        y2 = self.get_series()
        
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
            np.random.randn(100))

est_resp_H2O = (rb  * np.exp(Eo * (1 / (10 + 46.02) - 1 / (temp + 46.02))) +
                np.random.randn(100)) * vwc_func            
            
# Instantiate class without soil moisture
this_ER = ER(temp, est_resp)

# Get the results dictionary
params_dict = this_ER.get_fit()

# Instantiate class with soil moisture
wt_ER = ER(temp, est_resp_H2O, vwc)

this_dict = wt_ER.get_fit()