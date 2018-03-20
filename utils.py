#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 11:20:56 2018

@author: ian
"""

import pandas as pd
import pdb

def rename_df(df, external_names, internal_names):
    
    assert all(sorted(external_names.keys()) == sorted(internal_names.keys))
    pdb.set_trace()