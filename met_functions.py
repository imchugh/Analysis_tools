# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 11:16:19 2014

@author: imchugh
"""

# Humidity functions and conversions
def es(t=20):
    """Pass temperature (t) in celsius; returns saturation vapour pressure in kPa"""
    return 0.611*10**(7.5*t/(237.3+t))
    
def VPD(e,t):
    """Pass temperature (t) in celsius and vapour pressure (e) in kPa;
    returns vapour pressure deficit in kPa"""    
    return 0.611*10**(7.5*t/(237.3+t)) - e
    
def q_to_e(q,p=101.3):
    """Pass specific humidity (q) in g.kg-1, barometric pressure (p) in kpa;
    returns vapour pressure in kPa"""
    return (q*28.97)/28.97*p
   
