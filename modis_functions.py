# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 16:36:26 2016

@author: imchugh
"""

import numpy as np
from suds.client import *
import webbrowser
import datetime as dt
import sys

def get_products():
    webbrowser.open('https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table')

def get_sites():
    webbrowser.open('http://www.ozflux.org.au/monitoringsites/index.html')
    
def get_MODIS_subset(lat, lon, product, data_band, QC_band, 
                     requestStart = None, requestEnd = None):
    
    wsdlurl = 'https://modis.ornl.gov/cgi-bin/MODIS/soapservice/MODIS_soapservice.wsdl'
                
    # Get the available dates on the server
    client = Client(wsdlurl)
    
    # Get dates and set start and end dates
    server_datestr_arr = np.array(client.service.getdates(lat, lon, product))
    if requestStart or requestEnd:
        all_dates_arr = np.array([dt.datetime.strptime(i[1:], '%Y%j').date() 
                                  for i in server_datestr_arr])
        if not requestStart: 
            requestStart = all_dates_arr[0]
        if not requestEnd: 
            requestEnd = all_dates_arr[-1]
        target_dates_index = np.where((all_dates_arr >= requestStart) & 
                                      (all_dates_arr <= requestEnd))
        server_datestr_arr = server_datestr_arr[target_dates_index]
                                  
    # Iterate through dates and get data
    date_list = []
    data_list = []
    QC_list = []
    for date in server_datestr_arr:
    
        try:
            QC = client.service.getsubset(lat, lon, product, QC_band, 
                                          date, date, 
                                          0, 0)
            data = client.service.getsubset(lat, lon, product, data_band, 
                                            date, date, 
                                            0, 0)
        except Exception, e:
            print 'ORNL DAAC Server error with the following message: '
            print e[0]
            quit
        this_QC_str_list = QC.subset[0].split(',')
        this_QC = '{0:08b}'.format(int(this_QC_str_list[-1]))
        QC_list.append(this_QC)
        this_data_str_list = data.subset[0].split(',')
        this_date = this_data_str_list[2]
        date_list.append(dt.datetime.strptime(this_date[1:], '%Y%j').date())
        data_list.append(int(this_data_str_list[-1]) * data.scale)
    
    # Build dictionary
    return {'product': product,
            'band': data_band,
            'date': np.array(date_list),
            'QC': np.array(QC_list),
            'data': np.array(data_list)}