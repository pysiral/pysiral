# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 18:04:43 2015

@author: Stefan
"""

from dateutil import parser as dtparser


def parse_datetime_str(dtstr):
    """ Converts a time string to a datetime object using dateutils """
    return dtparser.parse(dtstr)
