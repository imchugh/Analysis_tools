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
   
# Estimate clear sky radiation and optimise k for site obs data
def Insol_calc(date_time,k,GMT_zone,latit,longit):
    
    DOY=date_time.timetuple().tm_yday
    hour=date_time.hour
    minute=date_time.minute
    
    # Calculate equation of time correction, solar noon, declination and TOA radiation
    EqofTime=0.17*np.sin(4*np.pi*(DOY-80)/373)-0.129*np.sin(2*np.pi*(DOY-8)/355) # DiLaura (1984)
    solar_noon=12+(GMT_zone*15.0-longit)/360*24-EqofTime # Me
    decl=np.radians(23.4)*np.sin((DOY+284)/365.0*2*np.pi) # Oke (1987)
    TOArad=(1+0.034*np.cos(DOY/365.25*2*np.pi))*1367.0 # Duffie and Beckman (1980)
    pdb.set_trace()
    
    # Create an hour angle array for each minute of day and each day of year
    array_h=np.tile(np.linspace(0,1439.0/1440*24,num=1440),(len(DOY),1))
    array_h=abs(np.radians((array_solar_noon.reshape(len(DOY),1)-array_h)*15))
    
    # Duplicate declination array for each time of day
    array_decl=np.tile(array_decl,(1440,1)).T

    # Calculate zenith angles
    array_z=np.arccos(np.sin(np.radians(lat_decdeg))*np.sin(array_decl)+
            np.cos(np.radians(lat_decdeg))*np.cos(array_decl)*np.cos(array_h))
    array_z_msk=np.ma.masked_greater_equal(array_z,np.pi/2) # Mask night values    
  
    # Calculate optical air mass term for all valid Z 
    array_m=(np.exp(-1*ALT_m/8343.5)/(np.cos(array_z_msk)+0.15*
            (np.degrees(90-array_z_msk)+3.855)**-1.253)) # Wunderlich (1972)           
    
    # Instantaneous clear sky surface radiation in Wm-2 for each minute of the day
    array_Kdown_clr_mins=array_TOArad.reshape(len(array_TOArad),1)*np.exp(-k*array_m)*np.cos(array_z_msk)
    
    # Aggregate one-minute instantaneous clear sky rad to period average
    array_Kdown_clr_hr=np.empty([len(DOY),1440/rec_length])
    for i in xrange(len(DOY)):
        array_temp=array_Kdown_clr_mins[i][:].reshape(1440/rec_length,rec_length) # Temporary bins
        array_Kdown_clr_hr[i][:]=np.ma.mean(array_temp,axis=1) # Average bin content  
    
    # Aggregate to daily
    array_Kdown_clr_daily=(array_Kdown_clr_hr*(rec_length*60.0/10**6)).sum(axis=1)
        
    if boolOutput==False:
        return array_Kdown_clr_daily # Result for optimisation
    else:
        return array_Kdown_clr_daily,array_Kdown_clr_hr    # Output of final data