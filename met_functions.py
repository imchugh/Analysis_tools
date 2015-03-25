# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 11:16:19 2014

@author: imchugh
"""
import pdb
import numpy as np
import datetime as dt

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
   
# Estimate clear sky radiation
def Insol_calc(date_time,k,GMT_zone,latit,longit,ALT_m):
    
    d = {}    
    
    DOY = date_time.timetuple().tm_yday
    hour = date_time.hour
    minute = date_time.minute
    
    # Calculate equation of time correction, solar noon, declination and TOA radiation
    EqofTime = 0.17 * np.sin(4 * np.pi * (DOY-80) / 373) - 0.129 * np.sin(2 * np.pi *(DOY-8) / 355) # DiLaura (1984)
    solar_noon = 12 + (GMT_zone * 15.0 - longit) / 360 * 24 - EqofTime # Me
    decl = np.radians(23.4) * np.sin((DOY + 284) / 365.0 * 2 * np.pi) # Oke (1987)
    TOArad = (1 + 0.034 * np.cos(DOY / 365.25 * 2 * np.pi)) * 1367.0 # Duffie and Beckman (1980)

    # Calculate hour angle    
    hr_angle =abs(np.radians((solar_noon - (minute/60.0 + hour)) * 15))

    # Calculate solar zenith angle
    zenith = np.arccos(np.sin(np.radians(latit)) * np.sin(decl) + 
             np.cos(np.radians(latit)) * np.cos(decl) * np.cos(hr_angle))
   
    # Check if night - if yes, return   
    if zenith > np.pi/2:
        return 0             
 
    # Calculate optical air mass term for all valid Z 
    m = (np.exp(-1 * ALT_m / 8343.5) / (np.cos(zenith) + 0.15 *
         (np.degrees(90 - zenith) + 3.855)** -1.253)) # Wunderlich (1972)
   
    # Instantaneous clear sky surface radiation in Wm-2
    Kdown = TOArad * np.exp(-k * m) * np.cos (zenith)
    
    d['solar_noon'] = 
    
    pdb.set_trace()    
  