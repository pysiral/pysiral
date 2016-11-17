# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 17:15:39 2016

@author: shendric
"""

from pysiral.config import options_from_dictionary
from pysiral.errorhandler import ErrorStatus

import numpy as np


class L2FreeboardAlgorithmBaseClass(object):

    def __init__(self):
        self.error = ErrorStatus()

    def set_options(self, **opt_dict):
        self._options = options_from_dictionary(**opt_dict)

    def get_freeboard(self, l2):
        freeboard, freeboard_uncertainty = self._get_freeboard(l2)
        return freeboard, freeboard_uncertainty


class SnowGeometricCorrection(L2FreeboardAlgorithmBaseClass):
    # TODO: Add functionality to use snow density
    """ Applies geometric corrections for snow wave progagation """

    def __init__(self):
        super(SnowGeometricCorrection, self).__init__()

    def _get_freeboard(self, l2):
        """ Compute the freeboard and its uncertainty """

        # Init parameter arrays
        shape = l2.afrb.shape
        freeboard = np.full(shape, np.nan, dtype=np.float32)
        freeboard_uncertainty = np.full(shape, np.nan, dtype=np.float32)

        # Compute the snow wave speed correction factor
        correction_factor = self._options.vacuum_light_speed_reduction
        geometric_correction = correction_factor * l2.snow_depth

        # Compute freeboard only for sea ice waveforms
        is_ice = l2.surface_type.sea_ice.indices
        freeboard[is_ice] = l2.afrb[is_ice] + geometric_correction[is_ice]
        freeboard_uncertainty[is_ice] = 0.0

        return freeboard, freeboard_uncertainty


def get_frb_algorithm(name):
    return globals()[name]()
