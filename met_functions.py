# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 11:16:19 2014

@author: imchugh
"""
import pdb
import numpy as np

# Humidity functions and conversions
def es(t=20):
    """Pass temperature (t) in celsius; returns saturation vapour pressure in kPa"""
    return 0.611*10**(7.5*t/(237.3+t))
    
def VPD(e,t):
    """Pass temperature (t) in celsius and vapour pressure (e) in kPa;
       returns vapour pressure deficit in kPa"""    
    return 0.611*10**(7.5*t/(237.3+t)) - e
    
def q_to_e(q, p=101.3):
    """Pass specific humidity (q) in g.kg-1, barometric pressure (p) in kpa;
       returns vapour pressure in kPa"""
    return (q*28.97)/28.97*p
    
def Ah_to_e(Ah, Ta):
    """Pass absolute humidity (Ah) in g.m-3 and air temperature in celsius;
       returns vapour pressure in kPa"""
    return Ah / 18 * (Ta + 273.15) * 8.3143 / 10**3

def Lv(t):
    """Pass temperature in C; returns Lv in J g-1;
    note: valid for range -25 to +40C"""
   
    return 2500.8 -2.36 * t + 0.0016 * t**2 - 0.00006 * t**3  
   
# Estimate clear sky radiation
def Insol_calc(date_time, GMT_zone, latit, longit, ALT_m, k):
    """
    Pass array of date_times in datetime format, GMT_zone in decimal hours, 
    lat and long in decimal degrees and altitude in m; extinction coefficient 
    is unitless, with range 0-1 \n
    Currently does not include solar disk diameter or refraction corrections
    
    Refs:
    
    DiLaura, D. L. (1984), IES Calculation Procedures Committee Recommended
    practice for the calculation of daylight availability, J. Illuminating
    Engineering Soc. of North America, 13(4), 381-392.
    
    Duffie, J. and W. Beckman (1980). Solar Engineering of Thermal Processes. 
    New York, John Wiley and Sons.
    
    Wunderlich, W. (1972), Heat and Mass Transfer between a Water Surface
    and the Atmosphere, Report No 14, Report Publication No. 0-6803,
    Water Resources Research Laboratory, TennesseeValleyAuthority,Division
    of Water Control Planning, Engineering Laboratory, Norris, TN.
    """
    
    # Get date and time components
    DOY = np.array([i.timetuple().tm_yday for i in date_time])
    hour = np.array([i.hour for i in date_time])
    minute = np.array([i.minute for i in date_time])

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
    zenith_msk = np.ma.masked_greater_equal(zenith, np.pi / 2) # Mask night values

    # Calculate optical air mass term
    m = (np.exp(-1 * ALT_m / 8343.5) / (np.cos(zenith_msk) + 0.15 *
         (np.degrees(90 - zenith) + 3.855)** -1.253)) # Wunderlich (1972)
    # Instantaneous clear sky surface radiation in Wm-2 (Beer-Lambert variant)
    Kdown = (TOArad * np.exp(-k * m) * np.cos (zenith_msk)).filled(0)
    m = m.filled(np.nan)
    
    d = {}    
    d['solar_noon'] = solar_noon
    d['declination'] = decl
    d['TOA_radiation'] = TOArad
    d['zenith'] = zenith
    d['optical_air_mass'] = m
    d['Kdown'] = Kdown
    
    return d
  