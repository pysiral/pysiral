# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 17:15:39 2016

@author: shendric
"""

from pysiral.config import options_from_dictionary


class L2FreeboardAlgorithmBaseClass(object):

    def __init__(self):
        pass

    def set_options(self, **opt_dict):
        self._options = options_from_dictionary(**opt_dict)

    def get_freeboard(self, l2):
        freeboard, msg = self._get_freeboard(l2)
        return freeboard, msg


class SnowGeometricCorrection(L2FreeboardAlgorithmBaseClass):
    # TODO: Add functionality to use snow density
    """ Applies geometric corrections for snow wave progagation """

    def __init__(self):
        super(SnowGeometricCorrection, self).__init__()

    def _get_freeboard(self, l2):
        correction_factor = self._options.vacuum_light_speed_reduction
        geometric_correction = correction_factor * l2.snow_depth
        sea_ice_only = l2.surface_type.sea_ice.indices
        freeboard = l2.afrb
        freeboard[sea_ice_only] += geometric_correction[sea_ice_only]
        return freeboard, ""


def get_frb_algorithm(name):
    return globals()[name]()
