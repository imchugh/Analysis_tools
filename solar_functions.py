# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 15:59:38 2016

@author: imchugh
"""

import ephem
import datetime as dt
import numpy as np

def get_zenith(data_dict, lat, lon, alt, GMT_zone):

    # set local variable datetime
    date_time = data_dict['date_time']

    # Create and populate local observer (lat and long must be strings)
    obs = ephem.Observer()
    obs.elev = alt
    obs.lat = str(lat)
    obs.long = str(lon)

    # Specify body
    sun = ephem.Sun(obs)

    # Convert to UTC
    UTC_datetime = date_time - dt.timedelta(hours = GMT_zone)

    # Convert to array if single value
    try: 
        iter(UTC_datetime)
    except:
        UTC_datetime = np.array([UTC_datetime, ])

    # Get zenith for each date_time
    z_list = []

    for i, this_dt in enumerate(UTC_datetime):
        if 'T' in data_dict.keys():
            obs.temp = data_dict['T'][i]
        if 'P' in data_dict.keys():
            obs.pressure = data_dict['P'][i]
        obs.date = this_dt
        sun.compute(obs)
        z_list.append(sun.alt)
        
    z_array = np.array(z_list)
        
    return np.pi / 2 - z_array
        
# Estimate clear sky radiation
def Insol_calc(data_dict, GMT_zone, latit, longit, ALT_m, k, use_ephem = False):
    """
    Pass args: 
        1) data_dict: must contain equal-length arrays, with a minimum of 
                      key / value pair of:
                          - 'date_time': array of datetimes (python datetime);
                      additionally, if using pyephem for zenith calculation and
                      accurate refraction correction required, include two 
                      additional key / value pairs of: 
                          - 'T': array of temperatures in C (float)
                          - 'P': array of pressures in hPa (float)
        2) GMT_zone: time zone in decimal hours (int or float),
        3) latit: latitude in decimal degrees (int or float),
        4) longit: longitude in decimal degrees (int or float),
        5) alt: altitude in mASL (int or float); 
        6) k: extinction coefficient - unitless with range 0-1 (int or float)

    Optional kwargs:
        1) use_ephem - use pyephem instead of algoirthms here (boolean); if 
                       valid entries for T and P are found in data_dict, the 
                       atmospheric refraction correction will be more accurate
                       (otherwise defaults to )
    
    The algorithms used are from the references below:    
   
    DiLaura, D. L. (1984), IES Calculation Procedures Committee Recommended
    practice for the calculation of daylight availability, J. Illuminating
    Engineering Soc. of North America, 13(4), 381-392.
    
    Duffie, J. and W. Beckman (1980). Solar Engineering of Thermal Processes. 
    New York, John Wiley and Sons.
    
    Wunderlich, W. (1972), Heat and Mass Transfer between a Water Surface
    and the Atmosphere, Report No 14, Report Publication No. 0-6803,
    Water Resources Research Laboratory, TennesseeValleyAuthority,Division
    of Water Control Planning, Engineering Laboratory, Norris, TN.
    
    Note that solar disk diameter or refraction corrections are not yet 
    included; if u want these, set use_ephem to True
    """
    
    date_time = data_dict['date_time']    
    
    # Get date and time components
    try: 
        iter(date_time)
    except:
        date_time = [date_time]
        
    DOY = np.array([i.timetuple().tm_yday for i in date_time])
    hour = np.array([i.hour for i in date_time])
    minute = np.array([i.minute + 15 for i in date_time])

    # Calculate equation of time correction, solar noon, declination and TOA radiation
    EqofTime = (0.17 * np.sin(4 * np.pi * (DOY-80) / 373) 
                - 0.129 * np.sin(2 * np.pi *(DOY-8) / 355)) # DiLaura (1984)
    solar_noon = 12 + (GMT_zone * 15.0 - longit) / 360 * 24 - EqofTime # Me
    decl = np.radians(23.4) * np.sin((DOY + 284) / 365.0 * 2 * np.pi) # Oke (1987)
    TOArad = (1 + 0.034 * np.cos(DOY / 365.25 * 2 * np.pi)) * 1367.0 # Duffie and Beckman (1980)

    # Calculate hour angle    
    hr_angle =abs(np.radians((solar_noon - (minute/60.0 + hour)) * 15))

    # Calculate solar zenith angle
    if use_ephem:
        date_time = np.array([this_dt - dt.timedelta(minutes = 15) 
                              for this_dt in date_time])
        zenith = get_zenith(data_dict, latit, longit, ALT_m, GMT_zone)
    else:
        zenith = np.arccos(np.sin(np.radians(latit)) * np.sin(decl) + 
                 np.cos(np.radians(latit)) * np.cos(decl) * np.cos(hr_angle))
    zenith_msk = np.ma.masked_greater_equal(zenith, np.pi / 2) # Mask night values

    # Calculate optical air mass term
    m = (np.exp(-1 * ALT_m / 8343.5) / (np.cos(zenith_msk) + 0.15 *
         (np.degrees(np.pi - zenith) + 3.855)** -1.253)) # Wunderlich (1972)
    
    # Instantaneous clear sky surface radiation in Wm-2 (Beer-Lambert variant)
    Kdown = (TOArad * np.exp(-k * m) * np.cos (zenith_msk)).filled(0)
    m = m.filled(np.nan)

    # Make result dict    
    d = {}    
    d['solar_noon'] = solar_noon
    d['declination'] = decl
    d['TOA_radiation'] = TOArad
    d['zenith'] = zenith
    d['optical_air_mass'] = m
    d['Kdown'] = Kdown
    
    return d    

