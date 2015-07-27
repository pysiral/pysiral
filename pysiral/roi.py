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
        self._options = TreeDict.fromdict(opt_dict, expand_nested=True)


class LowerLatLimit(ROIBase):
    """ Lower Latitude Treshold """
    def __init__(self):
        super(LowerLatLimit, self).__init__()

    def get_roi_list(self, longitude, latitude):
        if self._options.latitude_threshold > 0.0:  # Northern Hemisphere
            return np.where(latitude >= self._options.latitude_threshold)[0]
        if self._options.latitude_threshold <= 0.0:  # Southern Hemisphere
            return np.where(latitude <= self._options.latitude_threshold)[0]
