# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 17:15:39 2016

@author: shendric
"""

import numbers
import numpy as np
from attrdict import AttrDict
from pysiral.errorhandler import ErrorStatus


class L2FreeboardAlgorithmBaseClass(object):

    def __init__(self):
        self.error = ErrorStatus()

    def set_options(self, **opt_dict):
        self._options = AttrDict(**opt_dict)

    def get_freeboard(self, l1b, l2):
        freeboard, freeboard_uncertainty = self._get_freeboard(l1b, l2)
        return freeboard, freeboard_uncertainty


class SnowGeometricCorrection(L2FreeboardAlgorithmBaseClass):
    """
    Computes the freeboard from radar freeboard by application of geometric corrections for snow wave progagation
    only. No form of penetration correction is applied based on the assumption that radar freeboard is the
    sea-ice freeboard with missing correction for the slower wave propagation speed inside the snow layer.

    The correction can be be computed in two cases:

        1. [Static Correction] A correction factor 0 < c < 1 need to be supplied that describes
           the fraction of slower wave propagation speed in the snow layer compared to the
           vacuum speed. In this case the geometric correction is c * snow_depth which will
           be added to the radar freeboard.

        2. [Density Dependend Correction] In this case, the string "mallett2020" indicates that
           the geometric correction should be computed as function of density according to
           Mallett et al., 2020.


    Reference:

        Mallett, R. D. C., Lawrence, I. R., Stroeve, J. C., Landy, J. C., and Tsamados, M.:
        Brief communication: Conventional assumptions involving the speed of radar waves in snow
        introduce systematic underestimates to sea ice thickness and seasonal growth rate estimates,
        The Cryosphere, 14, 251–260, https://doi.org/10.5194/tc-14-251-2020, 2020.


    Updates:

        - [July 2020] added the option to use density-dependent formulation from Mallett et al., 2020.


    Configuration Example (Level-2 Processor Definition):

        frb:
            pyclass: SnowGeometricCorrection
            options:
                vacuum_light_speed_reduction: [0.22|mallett2020]

    """

    def __init__(self):
        """
        Initialize the class.
        """
        super(SnowGeometricCorrection, self).__init__()

    def _get_freeboard(self, l1b, l2):
        """
        Compute and apply the geometric correction to the radar freeboard to compute sea-ice freeboard.
        Only l2 data is needed for this class.
        :param l1b: The Level-1b data class
        :param l2: The Level-2 data class
        :return: freeboard and freeboard uncertainty as arrays with expected dimensions
        """

        # Init parameter arrays
        shape = l2.afrb.shape
        freeboard = np.full(shape, np.nan, dtype=np.float32)
        freeboard_uncertainty = np.full(shape, np.nan, dtype=np.float32)

        # Compute the snow wave speed correction factor
        correction_factor = self.get_correction_factor(l2)
        geometric_correction = correction_factor * l2.sd

        # Compute freeboard only for sea ice waveforms
        is_ice = l2.surface_type.sea_ice.indices
        freeboard[is_ice] = l2.afrb[is_ice] + geometric_correction[is_ice]

        # Compute the uncertainty
        deriv_snow = correction_factor
        uncertainty = np.sqrt((deriv_snow*l2.sd.uncertainty)**2. + l2.afrb.uncertainty**2.)
        freeboard_uncertainty[is_ice] = uncertainty[is_ice]

        # All done, return values
        return freeboard, freeboard_uncertainty

    def get_correction_factor(self, l2):
        """
        Returns the correction factor that is multiplied by snow depth. The difference between the
        possible options (see class docstring) is implemented in this method
        :param l2: Level-2 data class
        :return: Correction factor c (correction = c * snow_depth). Scalar or array
        """

        # Get the option value from the Level-2 processor definition file
        config_value = self._options.vacuum_light_speed_reduction

        # Case 1: Number is float value
        # -> return as is
        if isinstance(config_value, numbers.Real):
            return config_value

        # Case 2: Value a string and the correct option for the Mallett et al., 2020 correction.
        # Note
        c = 299792458.   # vacuum light speed
        if config_value == "mallett2020":
            c_s = c * (1. + 0.51 * l2.sdens[:]/1000.)**-1.5
            correction_factor = c/c_s - 1.
            return correction_factor

        # Invalid option: Raise a ValueError for the moment.
        else:
            raise ValueError("Invalid option `vacuum_light_speed_reduction` [{}]".format(str(config_value)))


        correction_factor = self._options.vacuum_light_speed_reduction


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
        self._options = AttrDict(**opt_dict)

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
