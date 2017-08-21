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

    def get_freeboard(self, l1b, l2):
        freeboard, freeboard_uncertainty = self._get_freeboard(l1b, l2)
        return freeboard, freeboard_uncertainty


class SnowGeometricCorrection(L2FreeboardAlgorithmBaseClass):
    # TODO: Add functionality to use snow density
    """ Applies geometric corrections for snow wave progagation """

    def __init__(self):
        super(SnowGeometricCorrection, self).__init__()

    def _get_freeboard(self, l1b, l2):
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

        # XXX: This is not correct yet
        deriv_snow = correction_factor
        uncertainty = np.sqrt((deriv_snow*l2.snow_depth.uncertainty)**2. +
                              l2.afrb.uncertainty**2.)
        freeboard_uncertainty[is_ice] = uncertainty[is_ice]

        return freeboard, freeboard_uncertainty


class SnowFreeboardAssumption(L2FreeboardAlgorithmBaseClass):
    # TODO: Add functionality to use snow density
    """ Applies geometric corrections for snow wave progagation """

    def __init__(self):
        super(SnowFreeboardAssumption, self).__init__()

    def _get_freeboard(self, l1b, l2):
        """ Compute the freeboard and its uncertainty """

        # Init parameter arrays
        shape = l2.afrb.shape
        freeboard = np.full(shape, np.nan, dtype=np.float32)
        freeboard_uncertainty = np.full(shape, np.nan, dtype=np.float32)

        # Transfer freeboard only for sea ice waveforms
        is_ice = l2.surface_type.sea_ice.indices
        freeboard[is_ice] = l2.afrb[is_ice]

        # Also just transfer "radar" freeboard uncertainty
        freeboard_uncertainty[is_ice] = l2.afrb.uncertainty[is_ice]

        return freeboard, freeboard_uncertainty


class L2RadarFreeboardAlgorithmBaseClass(object):

    def __init__(self):
        self.error = ErrorStatus()

    def set_options(self, **opt_dict):
        self._options = options_from_dictionary(**opt_dict)

    def get_radar_freeboard(self, l1b, l2):
        rfrb, rfrb_unc = self._get_radar_freeboard(l1b, l2)
        return rfrb, rfrb_unc


class RadarFreeboardDefault(L2RadarFreeboardAlgorithmBaseClass):
    """
    Default Class for computing Radar Freeboard
    (simple computation based on elevation, mean sea surface and
     sea surface anomaly)
    """

    def __init__(self):
        super(RadarFreeboardDefault, self).__init__()

    def _get_radar_freeboard(self, l1b, l2):
        """ Compute the radar freeboard and its uncertainty """

        # radar freeboard is simple
        rfrb = l2.elev - l2.mss - l2.ssa

        # radar freeboard uncertainty is not
        rfrb_unc = self._get_radar_freeboard_uncertainty(l2)

        return rfrb, rfrb_unc

    def _get_radar_freeboard_uncertainty(self, l2):
        """
        Get the radar freeboard uncertainty based on  error propagation
        from assumed uncorrelated uncertainties of
        - elevation (range)
        - sea surface anomaly
        mss uncertainties is ignored, respectively part of ssa uncertainty
        """
        rfrb_unc = np.sqrt((l2.elev.uncertainty)**2. +
                           (l2.ssa.uncertainty)**2.)
        return rfrb_unc


def get_frb_algorithm(name):
    return globals()[name]()
