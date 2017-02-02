# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 16:20:41 2015

@author: imchugh
"""

import warnings
import logging
import os
import datetime as dt
import data_formatting as dt_fm

reload(logging)

def this_main():
    log_dir = os.path.expanduser('~')
    log_file = 'Test'
    full_fname = os.path.join(log_dir, log_file)
    logging.basicConfig(filename = full_fname, filemode = 'w', level = logging.INFO)
    cur_time = dt.datetime.now()
    time_str = dt.datetime.strftime(cur_time, '%Y-%m-%d %H:%M:%S')
    
    warnings.showwarning = dt_fm.send_warnings_to_log
    
    warnings.warn('This is a warning message, received on: {0}'.format(time_str))
    logging.debug('This message should go to the log file')
    logging.info('So should this')
    logging.warning('And this, too')
            
