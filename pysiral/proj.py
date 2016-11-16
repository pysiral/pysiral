# -*- coding: utf-8 -*-
"""
Created on Thu Apr 07 20:02:58 2016

@author: Stefan
"""


class BaseProjection(object):

    def __init__(self):
        self.projection = {}

    @property
    def projection_keyw(self):
        return self.projection


class EASE2North(BaseProjection):

    def __init__(self):
        super(EASE2North, self).__init__()
        self.projection = {
            "proj": "laea",
            "lat_0": 90.0,
            "lon_0": 0.0,
            "ellps": "WGS84",
            "datum": "WGS84",
            "units": "m"}


class EASE2South(BaseProjection):

    def __init__(self):
        super(EASE2South, self).__init__()
        self.projection = {
            "proj": "laea",
            "lat_0": -90.0,
            "lon_0": 0.0,
            "ellps": "WGS84",
            "datum": "WGS84",
            "units": "m"}
