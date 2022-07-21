# -*- coding: utf-8 -*-
"""
@author: Stefan Hendricks

A python module dedicated to freeboard estimation.

"""

import numbers
import numpy as np

from pysiral.l2data import Level2Data
from pysiral.l1bdata import Level1bData
from pysiral.l2proc.procsteps import Level2ProcessorStep


class SnowGeometricCorrection(Level2ProcessorStep):
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

    def __init__(self, *args, **kwargs):
        """
        Initialize the class.
        """
        super(SnowGeometricCorrection, self).__init__(*args, **kwargs)

    def execute_procstep(self, l1b, l2):
        """
        Compute and apply the geometric correction to the radar freeboard to compute sea-ice freeboard.
        Only l2 data is needed for this class and the object will be modified in-place
        :param l1b: The Level-1b data class
        :param l2: The Level-2 data class
        :return: error status flag
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

        # Set the values in the Level-2 data object
        l2.frb.set_value(freeboard)
        l2.frb.set_uncertainty(freeboard_uncertainty)

        # Register the snow geometric correction as a auxiliary variable
        is_not_ice = np.logical_not(l2.surface_type.sea_ice.flag)
        geometric_correction[is_not_ice] = np.nan
        l2.set_auxiliary_parameter("sgcor", "snow_geometric_correction", geometric_correction)

        # Return the error flag (where frb is NaN)
        return np.isnan(l2.frb[:])

    def get_correction_factor(self, l2):
        """
        Returns the correction factor that is multiplied by snow depth. The difference between the
        possible options (see class docstring) is implemented in this method
        :param l2: Level-2 data class
        :return: Correction factor c (correction = c * snow_depth). Scalar or array
        """

        # Get the option value from the Level-2 processor definition file
        config_value = self.cfg.options.vacuum_light_speed_reduction

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

    @property
    def l2_input_vars(self):
        return ["afrb", "sd", "sdens", "surface_type"]

    @property
    def l2_output_vars(self):
        return ["frb", "sgcorr"]

    @property
    def error_bit(self):
        return self.error_flag_bit_dict["frb"]


class SnowFreeboardAssumption(Level2ProcessorStep):
    """
    Assumes the altimeter freeboard is the snow freeboard
    """

    def __init__(self, *args, **kwargs):
        super(SnowFreeboardAssumption, self).__init__(*args, **kwargs)

    def execute_procstep(self, l1b, l2):
        """
        Compute the freeboard and its uncertainty
        :param l1b:
        :param l2:
        :return:
        """

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

    @property
    def l2_input_vars(self):
        return ["afrb", "surface_type"]

    @property
    def l2_output_vars(self):
        return ["frb"]


class RadarFreeboardDefault(Level2ProcessorStep):
    """
    Default Class for computing radar freeboard based on elevation, mean sea surface and
    sea level anomaly
    """

    def __init__(self, *args, **kwargs):
        super(RadarFreeboardDefault, self).__init__(*args, **kwargs)

    def execute_procstep(self, l1b, l2):
        """
        Compute the radar freeboard and its uncertainty and modifies the Level-2 data object in-place
        :param l1b:
        :param l2:
        :return:
        """

        shape = l2.elev.shape
        radar_freeboard = np.full(shape, np.nan, dtype=np.float32)
        radar_freeboard_uncertainty = np.full(shape, np.nan, dtype=np.float32)

        # Compute freeboard only for sea ice waveforms
        is_ice = l2.surface_type.sea_ice.indices

        # Compute radar freeboard as the simple difference of
        # surface elevation - ssh (= mss + sla)
        radar_freeboard[is_ice] = l2.elev[is_ice] - l2.mss[is_ice] - l2.sla[is_ice]

        # Compute radar freeboard uncertainty assuming that the uncertainty
        # component of the MSS is negligible
        radar_freeboard_uncertainty[is_ice] = np.sqrt(l2.elev.uncertainty[is_ice] ** 2. + l2.sla.uncertainty[is_ice] ** 2.)

        # Add parameters to Level-2 object
        # NOTE: Conventions of the Level-2 Processor for radar freeboard variable id
        #       are `afrb` (altimeter freeboard)
        l2.afrb.set_value(radar_freeboard)
        l2.afrb.set_uncertainty(radar_freeboard_uncertainty)

        # Return the error flag (where frb is NaN)
        return np.isnan(l2.afrb[:])

    @property
    def l2_input_vars(self):
        return ["elev", "mss", "sla"]

    @property
    def l2_output_vars(self):
        return ["afrb"]

    @property
    def error_bit(self):
        return self.error_flag_bit_dict["frb"]


class LaforgeTFMR50PPCorrection(Level2ProcessorStep):
    """
    Class implementing the pulse-peakiness based range correction for the
    TFMRA50 described in

    > Laforge, A., Fleury, S., Dinardo, S., Garnier, F., Remy, F., Benveniste, J., Bouffard, J., & Verley, J. (2021).
     Toward improved sea ice freeboard observation with SAR altimetry using the physical retracker SAMOSA+.
     In Advances in Space Research (Vol. 68, Issue 2, pp. 732–745). Elsevier BV.
     https://doi.org/10.1016/j.asr.2020.02.001

    NOTE: This range correction is not applicable for CryoSat-2 data and the use of a TFMRA50
    """

    def __init__(self, *args, **kwargs):
        super(LaforgeTFMR50PPCorrection, self).__init__(*args, **kwargs)

    def execute_procstep(self, l1b: Level1bData, l2: Level2Data) -> None:
        """
        Correct range values based on pulse peakiness according the Laforge et al. 2010.

        :param l1b:
        :param l2:
        :return:
        """

        pp = l2.pulse_peakiness
        range_correction = -1390. * pp**3. + 339 * pp**2. - 28.4 * pp + 0.994
        pp_threshold_idx = np.where(pp > 0.3)[0]
        range_correction[pp_threshold_idx] = 0.35
        l2.elev[:] += range_correction
        l2.set_auxiliary_parameter("pp_rc", "pp_range_correction", range_correction)

        breakpoint()

    @property
    def l2_input_vars(self):
        return ["pulse_peakiness", "elev"]

    @property
    def l2_output_vars(self):
        return ["pp_range_correction"]

    @property
    def error_bit(self):
        return self.error_flag_bit_dict["range_correction"]
