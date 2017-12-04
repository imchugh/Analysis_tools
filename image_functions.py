#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 12:32:54 2017

@author: ian
"""

import datetime as dt
import os
from PIL import Image
from PIL.ExifTags import TAGS
import pdb

def get_exif(f, item_list = None):
    ret_dict = {}
    i = Image.open(f)
    info = i._getexif()
    if not item_list is None:
        for tag, value in info.items():
            decoded = TAGS.get(tag, tag)
            try:
                assert decoded in item_list
                ret_dict[decoded] = value
            except AssertionError:
                pass
        return ret_dict
    else:
        for tag, value in info.items():
            decoded = TAGS.get(tag, tag)
            ret_dict[decoded] = value
        return ret_dict