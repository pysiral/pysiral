# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 13:02:45 2015

@author: Stefan
"""

from treedict import TreeDict
import numpy as np


class ROIBase(object):
    """ Base class for Region of Interest Definitions """

    def __init__(self):
        self._options = None

    def set_options(self, **opt_dict):
        options = self.default
        options.update(opt_dict)
        self._options = TreeDict.fromdict(options, expand_nested=True)


class LowerLatLimit(ROIBase):
    """ Lower Latitude Treshold """
    def __init__(self):
        super(LowerLatLimit, self).__init__()
        self.default = {"latitude_threshold": 60.0}

    def get_roi_list(self, longitude, latitude):
        if self._options.latitude_threshold > 0.0:  # Northern Hemisphere
            return np.where(latitude >= self._options.latitude_threshold)[0]
        if self._options.latitude_threshold <= 0.0:  # Southern Hemisphere
            return np.where(latitude <= self._options.latitude_threshold)[0]

    def get_latitude_range(self):
        if self._options.latitude_threshold > 0.0:  # Northern Hemisphere
            return [self._options.latitude_threshold, 90.0]
        if self._options.latitude_threshold <= 0.0:  # Southern Hemisphere
            return [-90.0, self._options.latitude_threshold]


class LatitudeLongitudeBox(ROIBase):
    """ Box defined by latitude and longitude range """
    def __init__(self):
        super(LatitudeLongitudeBox, self).__init__()
        self.default = {"lat_range": [-91, 90], "lon_range": [-181, 181]}

    def get_roi_list(self, longitude, latitude):
        in_lon = np.logical_and(longitude >= self._options.lon_range[0],
                                longitude <= self._options.lon_range[1])
        in_lat = np.logical_and(latitude >= self._options.lat_range[0],
                                latitude <= self._options.lat_range[1])
        in_roi = np.where(np.logical_and(in_lon, in_lat))[0]
        return in_roi


def get_roi_class(name):
    return globals()[name]()

