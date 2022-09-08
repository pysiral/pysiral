# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 15:48:58 2015

@author: Stefan Hendricks

TODO: move this to it own module and separate retrackers into different files
"""

# Utility methods for retracker:
import os.path

from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import bottleneck as bn

import sys
import time
import numpy as np
from loguru import logger
from attrdict import AttrDict
import logging
from pysiral import InterceptHandler

# cythonized bottleneck functions for cTFMRA
try:
    from .bnfunc.cytfmra import (cytfmra_findpeaks, cytfmra_interpolate,
                                 cytfmra_wfm_noise_level, cytfmra_normalize_wfm)
    CYTFMRA_OK = True
except ImportError:
    logger.error("Cannot import cytfmra")
    CYTFMRA_OK = False

try:
    from samosa.sampy import SAMOSA as initialize_SAMOSAlib
    from samosa.sampy import initialize_epoch, compute_ThNEcho
    from samosa.help_functions import calc_sigma0, func_wind_speed, calc_ssb_j2Table

    logger.info("SAMOSA retracker loaded from the environment")
    logging.getLogger('samosa.sampy').addHandler(InterceptHandler())
    logging.getLogger('samosa.sampy').setLevel('INFO')
    SAMOSA_OK = True
except ImportError:
    logger.error("Unable to import the SAMOSA retracker. Has it been installed?")
    SAMOSA_OK = False

from typing import Tuple

from pysiral.core.flags import ANDCondition, FlagContainer
from pysiral.l2proc.procsteps import Level2ProcessorStep


class BaseRetracker(object):
    """
    Main Retracker Class (all retrackers must be of instance BaseRetracker)
    # TODO: API clean-up is sorely needed.
    """

    def __init__(self):
        self._indices = None
        self._classifier = None
        self._l1b = None
        self._l2 = None
        self._range = None
        self._power = None
        self._options = AttrDict()

        # Dictionary containing potential auxiliary output variables of
        # the retracker algorithm that will be transferred to the l2
        # data object.
        # NOTE: This is a bit clunky because retracker algorithm will
        # not be run on all the waveforms, so a bit extra work is needed
        # that the output a) has the correct dimension and b) does not
        # overwrite itself, if the same algorithm is called consecutively.
        self.auxdata_output = []

    def set_options(self, **opt_dict):
        # TODO: Create options object, respectively use __init__
        self._options = AttrDict(opt_dict)

    def set_indices(self, indices):
        # TODO: Validation
        self._indices = indices

    def set_classifier(self, classifier):
        self._classifier = classifier

    def init(self, n_records):
        # TODO: Move to __init__
        self._create_default_properties(n_records)

    def register_auxdata_output(self, var_id, var_name, value, uncertainty=None):
        """
        Add an auxiliary parameter, that will be transferred to the l2data object after retracking
        """
        self.auxdata_output.append([var_id, var_name, value, uncertainty])

    def retrack(self, l1b, l2):

        # Store the pointer to l1b and l2 data objects
        self._l1b = l1b
        self._l2 = l2

        # Initialize the retracked range with an NaN array
        # -> only waveforms of type "surface_type" will retracked
        self._create_default_properties(l1b.n_records)

        # Give the Retracker the possibility to create additional
        # data arrays (All retracker must have this method)
        self.create_retracker_properties(l1b.n_records)

        # Set the indices for the given surface type
        # Check if there is something to do
        if self._indices is None:
            self._indices = np.arange(l1b.n_records)
        if len(self._indices) == 0:
            return False

        # Loop over each waveform of given surface type
        self.l2_retrack(l1b.waveform.range, l1b.waveform.power, self._indices,
                        l1b.waveform.radar_mode, l1b.waveform.is_valid)
        return True

    def l2_retrack(self, rng, pwr, indices, radar_mode, is_valid):
        """
        Abstract method, not to be called directly but expected to be overwritten
        by the child class
        :return:
        """
        raise NotImplementedError("BaseRetracker.l2_retrack should not be called directly")

    def get_l1b_parameter(self, data_group, parameter_name):
        """ Get any valid level-2 paremeter name """
        try:
            return self._l1b.get_parameter_by_name(data_group, parameter_name)
        except (KeyError, AttributeError):
            return None

    def get_l2_parameter(self, parameter_name):
        """ Get any valid level-2 paremeter name """
        try:
            return self._l2.get_parameter_by_name(parameter_name)
        except (KeyError, AttributeError):
            return None

    def _create_default_properties(self, n_records):
        # XXX: Currently only range and status (False: ok)
        for parameter in ["_range", "_power"]:
            setattr(self, parameter, np.full(n_records, np.nan))
        self._uncertainty = np.full(n_records, 0.0, dtype=np.float32)
        self._flag = np.zeros(shape=n_records, dtype=np.bool)

    def create_retracker_properties(self, n_records):
        # Will have to be overwritten
        pass

    @property
    def range(self):
        return self._range

    @property
    def uncertainty(self):
        return self._uncertainty

    @property
    def power(self):
        return self._power

    @property
    def indices(self):
        return self._indices

    @property
    def error_flag(self):
        return FlagContainer(self._flag)


class Level2RetrackerContainer(Level2ProcessorStep):
    """
    The interface for the Level-2 processor for all retrackers
    """

    def __init__(self, *args, **kwargs):
        """
        Initialize the instance
        :param args:
        :param kwargs:
        """
        super(Level2RetrackerContainer, self).__init__(*args, **kwargs)

    def execute_procstep(self, l1b, l2):
        """
        Mandatory class
        :param l1b:
        :param l2:
        :return:
        """

        # Get the error status
        error_status = self.get_clean_error_status(l2.n_records)

        # Retracker are surface type dependent
        # -> Loop over all requested surface types in the
        #    l2 processor definition file
        for surface_type, retracker_def in list(self.cfg.options.items()):

            # Check first if there are any waveforms of the requested
            # surface type
            surface_type_flag = l2.surface_type.get_by_name(surface_type)
            if surface_type_flag.num == 0:
                logger.info(f"- no waveforms of type {surface_type}")
                continue

            # Benchmark retracker performance
            timestamp = time.time()

            # Retrieve the retracker associated with surface type from the l2 settings
            retracker = get_retracker_class(retracker_def["pyclass"])

            # Set options (if any)
            if retracker_def["options"] is not None:
                retracker.set_options(**retracker_def["options"])

            # set subset of waveforms
            retracker.set_indices(surface_type_flag.indices)

            # Add classifier data (some retracker need that)
            retracker.set_classifier(l1b.classifier)

            # Start the retracker procedure
            retracker.retrack(l1b, l2)

            # Update range value of the Level-2 data object
            l2.update_retracked_range(retracker)

            # XXX: Let the retracker return other parameters?
            l2.radar_mode = l1b.waveform.radar_mode

            # retrieve potential error status and update surface type flag
            if retracker.error_flag.num > 0:
                l2.surface_type.add_flag(retracker.error_flag.flag, "invalid")
            logger.info("- Retrack class %s with %s in %.3f seconds" % (
                surface_type, retracker_def["pyclass"],
                time.time()-timestamp))

        return error_status

    @property
    def l2_input_vars(self):
        return []

    @property
    def l2_output_vars(self):
        return ["radar_mode", "range", "elev"]

    @property
    def error_bit(self):
        return self.error_flag_bit_dict["retracker"]


class SICCI2TfmraEnvisat(BaseRetracker):
    """ Default Retracker from AWI CryoSat-2 production system """

    DOCSTR = r"Threshold first maximum retracker (TFMRA)"

    def __init__(self):
        super(SICCI2TfmraEnvisat, self).__init__()

    def set_default_options(self):
        self.set_options(**self.default_options_dict)

    @property
    def default_options_dict(self):
        return {
            "threshold": dict(type="fixed", value=0.5),
            "offset": 0.0,
            "wfm_oversampling_factor": 10,
            "wfm_oversampling_method": "linear",
            "wfm_smoothing_window_size": [11, 11, 51],
            "first_maximum_normalized_threshold": [0.15, 0.15, 0.45],
            "first_maximum_local_order": 1
        }

    def create_retracker_properties(self, n_records):
        # None so far
        pass

    def l2_retrack(self, rng, wfm, indices, radar_mode, is_valid):
        """ API Calling method """

        # Get sigma0 information
        sigma0 = self.get_l1b_parameter("classifier", "sigma0")
        lew1 = self.get_l1b_parameter("classifier", "leading_edge_width_first_half")
        lew2 = self.get_l1b_parameter("classifier", "leading_edge_width_second_half")
        lew = lew1+lew2
        sitype = self._l2.sitype
        tfmra_threshold = self.get_tfmra_threshold(sigma0, lew, sitype, indices)

        for i in indices:

            # Remove artifact bins
            wfm_skipped = wfm[i, :]
            wfm_skipped[0:5] = 0

            # Get the filtered waveform, index of first maximum & norm
            filt_rng, filt_wfm, fmi, norm = self.get_filtered_wfm(rng[i, :], wfm_skipped, radar_mode[i])

            # first maximum finder might have failed
            if fmi == -1:
                self._range[i] = np.nan
                self._power[i] = np.nan
                return

            # Get track point and its power
            tfmra_range, tfmra_power = self.get_threshold_range(filt_rng, filt_wfm, fmi, tfmra_threshold[i])

            # Mandatory return function
            self._range[i] = tfmra_range + self._options.offset
            self._power[i] = tfmra_power * norm

        if "uncertainty" in self._options and self._options.uncertainty.type == "fixed":
            self._uncertainty[:] = self._options.uncertainty.value

    def get_tfmra_threshold(self, sigma0, lew, sitype, indices):

        # short link to options
        option = self._options.threshold

        threshold = np.full(sigma0.shape, np.nan)

        # legacy option where threshold is float in settings file
        if type(option) is float:
            threshold[indices] = option
            return threshold

        # fixed threshold
        if option.type == "fixed":
            threshold[indices] = option.value

        # depended on sigma0
        elif option.type == "sigma_func":
            value = np.zeros(sigma0.shape)
            for i, coef in enumerate(option.coef):
                value += coef * sigma0**i
            threshold[indices] = value[indices]

        # dependent in sea ice type and sigma0
        elif option.type == "sitype_sigma_func":
            value = np.zeros(sigma0.shape)
            for i, coef in enumerate(option.coef_fyi):
                value += coef * sigma0**i
            value_myi = np.zeros(sigma0.shape)
            for i, coef in enumerate(option.coef_myi):
                value_myi += coef * sigma0**i
            myi_list = np.where(sitype > 0.5)[0]
            value[myi_list] = value_myi[myi_list]
            threshold[indices] = value[indices]

        # dependent on sigma0 and leading-edge width 3rd order polynomial fit
        elif option.type == "poly_plane_fit":
            value = np.zeros(sigma0.shape)
            value += option.intercept
            for i, coef in enumerate(option.coef_lew):
                value += coef * lew**(i+1)
            for i, coef in enumerate(option.coef_sig0):
                value += coef * sigma0**(i+1)
            threshold[indices] = value[indices]

        else:
            msg = "treshold type not recognized: %s" % str(option.type)
            raise ValueError(msg)

        return threshold

    def get_preprocessed_wfm(self, rng, wfm, radar_mode, is_valid):
        """
        Returns the intermediate product (oversampled range bins,
        oversampled and filtered waveforms, indices of first maxima
        and peak power norm for custom applications
        """

        oversample_factor = self._options.wfm_oversampling_factor
        wfm_shape = (wfm.shape[0], wfm.shape[1]*oversample_factor)

        filt_rng = np.full(wfm_shape, np.nan)
        filt_wfm = np.full(wfm_shape, np.nan)
        fmi = np.full(wfm_shape[0], -1, dtype=np.int32)
        norm = np.full(wfm_shape[0], np.nan)

        for i in np.arange(wfm.shape[0]):
            if not is_valid[i]:
                continue
            result = self.get_filtered_wfm(rng[i, :], wfm[i, :], radar_mode[i])
            filt_rng[i, :] = result[0]
            filt_wfm[i, :] = result[1]
            fmi[i] = result[2]
            norm[i] = result[3]

        return filt_rng, filt_wfm, fmi, norm

    def get_thresholds_distance(self, rng, wfm, fmi, t0, t1):
        """
        Return the distance between two thresholds t0 < t1
        """
        width = np.full(rng.shape[0], np.nan, dtype=np.float32)
        for i in np.arange(rng.shape[0]):
            if fmi[i] is None:
                continue
            r0 = self.get_threshold_range(rng[i, :], wfm[i, :], fmi[i], t0)
            r1 = self.get_threshold_range(rng[i, :], wfm[i, :], fmi[i], t1)
            width[i] = r1[0] - r0[0]

        # some irregular waveforms might produce negative width values
        is_negative = np.where(width < 0.)[0]
        width[is_negative] = np.nan
        return width

    def get_filtered_wfm(self, rng, wfm, radar_mode):

        # Echo oversampling & filtering
        filt_rng, filt_wfm = self.filter_waveform(rng, wfm, radar_mode)

        # Normalize filtered waveform
        filt_wfm, norm = self.normalize_wfm(filt_wfm)

        # Get noise level in normalized units
        oversampling = self._options.wfm_oversampling_factor
        noise_level = wfm_get_noise_level(filt_wfm, oversampling)

        # Find first maxima
        # (needs to be above radar mode dependent noise threshold)
        fmnt = self._options.first_maximum_normalized_threshold[radar_mode]
        peak_minimum_power = fmnt + noise_level
        fmi = self.get_first_maximum_index(filt_wfm, peak_minimum_power)

        return filt_rng, filt_wfm, fmi, norm

    def filter_waveform(self, rng, wfm, radar_mode):
        """
        Return a filtered waveform: block filter smoothing of
        oversampled original waveform
        """

        # Parameter from options dictionary if omitted
        opt = self._options
        oversampling = opt.wfm_oversampling_factor
        interp_type = opt.wfm_oversampling_method
        window_size = opt.wfm_smoothing_window_size[radar_mode]

        # Waveform Oversampling
        n = len(rng)
        range_os = np.linspace(np.nanmin(rng), np.nanmax(rng), n*oversampling)
        interpolator = interp1d(rng, wfm, kind=interp_type)
        wfm_os = interpolator(range_os)

        # Smoothing
        wfm_os = smooth(wfm_os, window_size)

        return range_os, wfm_os

    @staticmethod
    def normalize_wfm(y):
        norm = np.nanmax(y)
        return y/norm, norm

    @staticmethod
    def get_first_maximum_index(wfm, peak_minimum_power):
        """
        Return the index of the first peak (first maximum) on
        the leading edge before the absolute power maximum.
        The first peak is only valid if its power exceeds a certain threshold
        """

        # Get the main maximum first
        absolute_maximum_index = np.argmax(wfm)

        # Find relative maxima before the absolute maximum
        try:
            peaks = findpeaks(wfm[0:absolute_maximum_index])
        except:
            return -1

        # Check if relative maximum are above the required threshold
        leading_maxima = np.where(wfm[peaks] >= peak_minimum_power)[0]

        # Identify the first maximum
        first_maximum_index = absolute_maximum_index
        if len(leading_maxima) > 0:
            # first_maximum_index = leading_maxima[0]
            first_maximum_index = peaks[leading_maxima[0]]

        return first_maximum_index

    @staticmethod
    def get_threshold_range(rng, wfm, first_maximum_index, threshold):
        """
        Return the range value and the power of the retrack point at
        a given threshold of the firsts maximum power
        """

        # get first index greater as threshold power
        first_maximum_power = wfm[first_maximum_index]

        # Get power of retracked point
        tfmra_power = threshold*first_maximum_power

        # Use linear interpolation to get exact range value
        potential_points = np.where(wfm[:first_maximum_index] > tfmra_power)[0]
        diff_points = np.diff(potential_points)
        diff_points_gap = np.where(diff_points > 1)[0]

        if len(potential_points) == 0:
            return np.nan, np.nan
        else:
            if len(diff_points_gap) > 0:
                retrack_point = potential_points[diff_points_gap[-1]+1]
            else:
                retrack_point = potential_points[0]

        i0, i1 = retrack_point-1, retrack_point
        gradient = (wfm[i1]-wfm[i0])/(rng[i1]-rng[i0])
        tfmra_range = (tfmra_power - wfm[i0]) / gradient + rng[i0]

        return tfmra_range, tfmra_power


class TFMRA(BaseRetracker):
    """
    Default Retracker from AWI CryoSat-2 production system (pure python implementation)
    FIXME: This implementation is deprecated and should be removed in favor of cTFMRA
    """

    DOCSTR = r"Threshold first maximum retracker (TFMRA)"

    def __init__(self):
        super(TFMRA, self).__init__()

    def set_default_options(self):
        self.set_options(**self.default_options_dict)

    @property
    def default_options_dict(self):
        default_options_dict = {
            "threshold": 0.5,
            "offset": 0.0,
            "wfm_oversampling_factor": 10,
            "wfm_oversampling_method": "linear",
            "wfm_smoothing_window_size": [11, 11, 51],
            "first_maximum_normalized_threshold": [0.15, 0.15, 0.45],
            "first_maximum_local_order": 1}
        return default_options_dict

    def create_retracker_properties(self, n_records):
        # None so far
        pass

    def l2_retrack(self, rng, wfm, indices, radar_mode, is_valid):
        """ API Calling method """

        for i in indices:

            # Get the filtered waveform, index of first maximum & norm
            filt_rng, filt_wfm, fmi, norm = self.get_filtered_wfm(
                rng[i, :], wfm[i, :], radar_mode[i])

            # first maximum finder might have failed
            if fmi == -1:
                self._range[i] = np.nan
                self._power[i] = np.nan
                return

            # Get track point and its power
            tfmra_threshold = self._options.threshold
            tfmra_range, tfmra_power = self.get_threshold_range(
                filt_rng, filt_wfm, fmi, tfmra_threshold)

            # Mandatory return function
            self._range[i] = tfmra_range + self._options.offset
            self._power[i] = tfmra_power * norm

        # Apply a radar mode dependent range bias if option is in
        # level-2 settings file
        if "range_bias" in self._options:
            for radar_mode_index in np.arange(3):
                indices = np.where(radar_mode == radar_mode_index)[0]
                if len(indices) == 0:
                    continue
                range_bias = self._options.range_bias[radar_mode_index]
                self._range[indices] -= range_bias

    def get_preprocessed_wfm(self, rng, wfm, radar_mode, is_valid):
        """
        Returns the intermediate product (oversampled range bins,
        oversampled and filtered waveforms, indices of first maxima
        and peak power norm for custom applications
        """

        oversample_factor = self._options.wfm_oversampling_factor
        wfm_shape = (wfm.shape[0], wfm.shape[1]*oversample_factor)

        filt_rng = np.full(wfm_shape, np.nan)
        filt_wfm = np.full(wfm_shape, np.nan)
        fmi = np.full(wfm_shape[0], -1, dtype=np.int32)
        norm = np.full(wfm_shape[0], np.nan)

        for i in np.arange(wfm.shape[0]):
            if not is_valid[i]:
                continue
            result = self.get_filtered_wfm(rng[i, :], wfm[i, :], radar_mode[i])
            filt_rng[i, :] = result[0]
            filt_wfm[i, :] = result[1]
            fmi[i] = result[2]
            norm[i] = result[3]

        return filt_rng, filt_wfm, fmi, norm

    def get_thresholds_distance(self, rng, wfm, fmi, t0, t1):
        """
        Return the distance between two thresholds t0 < t1
        """
        width = np.full(rng.shape[0], np.nan, dtype=np.float32)
        for i in np.arange(rng.shape[0]):
            if fmi[i] is None:
                continue
            r0 = self.get_threshold_range(rng[i, :], wfm[i, :], fmi[i], t0)
            r1 = self.get_threshold_range(rng[i, :], wfm[i, :], fmi[i], t1)
            width[i] = r1[0] - r0[0]

        # some irregular waveforms might produce negative width values
        is_negative = np.where(width < 0.)[0]
        width[is_negative] = np.nan
        return width

    def get_filtered_wfm(self, rng, wfm, radar_mode):

        # Echo oversampling & filtering
        filt_rng, filt_wfm = self.filter_waveform(rng, wfm, radar_mode)

        # Normalize filtered waveform
        filt_wfm, norm = self.normalize_wfm(filt_wfm)

        # Get noise level in normalized units
        oversampling = self._options.wfm_oversampling_factor
        noise_level = wfm_get_noise_level(filt_wfm, oversampling)

        # Find first maxima
        # (needs to be above radar mode dependent noise threshold)
        fmnt = self._options.first_maximum_normalized_threshold[radar_mode]
        peak_minimum_power = fmnt + noise_level
        fmi = self.get_first_maximum_index(filt_wfm, peak_minimum_power)

        return filt_rng, filt_wfm, fmi, norm

    def filter_waveform(self, rng, wfm, radar_mode):
        """
        Return a filtered waveform: block filter smoothing of
        oversampled original waveform
        """

        # Parameter from options dictionary if omitted
        opt = self._options
        oversampling = opt.wfm_oversampling_factor
        interp_type = opt.wfm_oversampling_method
        window_size = opt.wfm_smoothing_window_size[radar_mode]

        # Waveform Oversampling
        n = len(rng)
        range_os = np.linspace(np.nanmin(rng), np.nanmax(rng), n*oversampling)
        interpolator = interp1d(rng, wfm, kind=interp_type)
        wfm_os = interpolator(range_os)

        # Smoothing
        wfm_os = smooth(wfm_os, window_size)

        return range_os, wfm_os

    @staticmethod
    def normalize_wfm(y):
        norm = np.nanmax(y)
        return y/norm, norm

    @staticmethod
    def get_first_maximum_index(wfm, peak_minimum_power):
        """
        Return the index of the first peak (first maximum) on
        the leading edge before the absolute power maximum.
        The first peak is only valid if its power exceeds a certain threshold
        """

        # Get the main maximum first
        absolute_maximum_index = np.argmax(wfm)

        # Find relative maxima before the absolute maximum
        try:
            peaks = findpeaks(wfm[0:absolute_maximum_index])
        except:
            return -1

        # Check if relative maximum are above the required threshold
        leading_maxima = np.where(wfm[peaks] >= peak_minimum_power)[0]

        # Identify the first maximum
        first_maximum_index = absolute_maximum_index
        if len(leading_maxima) > 0:
            # first_maximum_index = leading_maxima[0]
            first_maximum_index = peaks[leading_maxima[0]]

        return first_maximum_index

    @staticmethod
    def get_threshold_range(rng, wfm, first_maximum_index, threshold):
        """
        Return the range value and the power of the retrack point at
        a given threshold of the firsts maximum power
        """

        # get first index greater as threshold power
        first_maximum_power = wfm[first_maximum_index]

        # Get power of retracked point
        tfmra_power = threshold*first_maximum_power

        # Use linear interpolation to get exact range value
        points = np.where(wfm[:first_maximum_index] > tfmra_power)[0]

        if len(points) == 0:
            return np.nan, np.nan

        i0, i1 = points[0]-1, points[0]
        gradient = (wfm[i1]-wfm[i0])/(rng[i1]-rng[i0])
        tfmra_range = (tfmra_power - wfm[i0]) / gradient + rng[i0]

        return tfmra_range, tfmra_power


class cTFMRA(BaseRetracker):
    """
    The default TFMRA retracker implementation using cythonized
    functions for numerical perfomrance
    """

    DOCSTR = r"Threshold first maximum retracker (TFMRA)"

    def __init__(self):
        super(cTFMRA, self).__init__()

    def set_default_options(self, specified_options=None):
        """
        Set the default options for the cTFMRA retracker. This can be modified
        by the specified_options keyword that will overwrite the values in the
        default options dictionary
        :param specified_options: (dict or dict-like) specific options
        :return:
        """

        # Create a copy of the default options
        options = dict(**self.default_options_dict)

        # Update options (if applicable)
        specified_options = specified_options if specified_options is not None else {}
        options.update(specified_options)

        # Set class options
        self.set_options(**options)

    @property
    def default_options_dict(self):
        default_options_dict = {
            # The power threshold for retracking point (normalized first maximum power)
            "threshold": 0.5,
            # Offset applied to the retracked range
            "offset": 0.0,
            # Indices of ranges bins for the noise power computation
            "noise_level_range_bin_idx": [0, 5],
            # Filter setting: Waveform oversampling factor
            "wfm_oversampling_factor": 10,
            # Filter setting: Waveform oversampling method
            "wfm_oversampling_method": "linear",
            # Filter setting: Oversampled waveform boxfilter size [lrm, sar, sarin]
            "wfm_smoothing_window_size": [11, 11, 51],
            # First maximum detection: Minimum power of first maximum (normalized maximum power)
            "first_maximum_normalized_threshold": [0.15, 0.15, 0.45],
            # First maximum detection: Peak detection setting
            "first_maximum_local_order": 1,
            # First maximum detection: First valid range bin index for first maximum
            "first_maximum_ignore_leading_bins": 0}
        return default_options_dict

    def create_retracker_properties(self, n_records):
        """
        Mandatory method, but unused
        :param n_records:
        :return:
        """
        pass

    def l2_retrack(self, rng, wfm, indices, radar_mode, is_valid):
        """
         API Calling method for retrackers.
        :param rng: (np.array, dim:(n_records, n_bins)
        :param wfm:
        :param indices:
        :param radar_mode:
        :param is_valid:
        :return:
        """

        # --- Additional input/output arrays ---

        # Get the threshold for each waveform and store the data
        tfmra_threshold = self.get_tfmra_threshold(indices)

        # Auxiliary output: The noise power determined by the TFMRA
        tfmra_noise_power = np.full(tfmra_threshold.shape, np.nan)

        # Auxiliary output: The index of the first maximum
        tfmra_first_maximum_index = np.full(tfmra_threshold.shape, -1, dtype=int)

        # --- TFMRA configuration options ---

        # Apply a fixed range offset, e.g. in the case of a known and constant retracker bias
        fixed_range_offset = self._options.offset

        # A factor by how much points the waveform should be oversampled
        # before smoothing
        oversampling_factor = self._options.wfm_oversampling_factor

        # The bin range of the waveform where the noise level should be detected
        # NOTE: In the config file the bins values refer to the actual range bins
        #       and not to the oversampled waveforms. This settings depends on
        #       the specific sensor and a default value should not be used.
        #       Thus a ValueError is raised if the option is missing
        noise_level_range_idx = self._options.get("noise_level_range_bin_idx", None)
        if noise_level_range_idx is None:
            raise ValueError("Missing ctfmra option: noise_level_range_bin_idx")

        # Waveforms from some altimeter have artefacts at the beginning of the
        # waveform that may be misclassified as a first maximum.
        # NOTE: This setting also depends on the sensor, but a default can be used.
        fmi_first_valid_idx = self._options.get("first_maximum_ignore_leading_bins", 0)
        # The value above with the oversampling factor
        fmi_first_valid_idx_filt = fmi_first_valid_idx * oversampling_factor

        # The window size for the box filter (radar mode dependant list)
        wfm_smoothing_window_size = self._options.wfm_smoothing_window_size

        # The power threshold for the first maximum (radar mode dependant list)
        first_maximum_normalized_threshold = self._options.first_maximum_normalized_threshold

        # --- Waveform Retracking ---
        # NOTE: This is a loop over all waveforms in the trajectory data but only those
        #       listed in the list of indices will be processed. The reason for this
        #       approach is that different retracker settings are needed for different
        #       surface types and radar modes and the initial thought was to avoid
        #       a copying of data.
        # TODO: This is a candidate for multi-processing

        for i in indices:

            # Do not retrack waveforms that are marked as invalid
            if not is_valid[i]:
                continue

            # Get the filtered waveform
            window_size = wfm_smoothing_window_size[radar_mode[i]]
            filt_rng, filt_wfm, norm = self.get_filtered_wfm(rng[i, :], wfm[i, :], oversampling_factor, window_size)

            # Get noise level in normalized units
            i0, i1 = [idx * oversampling_factor for idx in noise_level_range_idx]
            noise_level_normed = cytfmra_wfm_noise_level(filt_wfm, i0, i1)
            tfmra_noise_power[i] = noise_level_normed * norm

            # Find first maxima
            # (needs to be above radar mode dependent noise threshold)
            fmnt = first_maximum_normalized_threshold[radar_mode[i]]
            peak_minimum_power = fmnt + noise_level_normed
            fmi = self.get_first_maximum_index(filt_wfm, peak_minimum_power, fmi_first_valid_idx_filt)
            tfmra_first_maximum_index[i] = fmi

            # first maximum finder might have failed
            if fmi == -1:
                self._range[i] = np.nan
                self._power[i] = np.nan
                continue

            # Get track point and its power
            tfmra_range, tfmra_power, _ = self.get_threshold_range(filt_rng,
                                                                   filt_wfm,
                                                                   fmi,
                                                                   tfmra_threshold[i],
                                                                   fmi_first_valid_idx_filt)

            # Set the values
            self._range[i] = tfmra_range + fixed_range_offset
            self._power[i] = tfmra_power * norm

        # Add auxiliary variables to the l2 data
        self.register_auxdata_output("tfmrathr", "tfmra_threshold", tfmra_threshold)
        self.register_auxdata_output("tfmrafmi", "tfmra_first_maximum_index", tfmra_first_maximum_index)
        self.register_auxdata_output("tfmranp", "tfmra_noise_power", tfmra_noise_power)

        # Apply a radar mode dependent range bias if option is in
        # level-2 settings file
        if "range_bias" in self._options:
            for radar_mode_index in np.arange(3):
                indices = np.where(radar_mode == radar_mode_index)[0]
                if len(indices) == 0:
                    continue
                range_bias = self._options.range_bias[radar_mode_index]
                self._range[indices] -= range_bias

        if "uncertainty" in self._options:
            if self._options.uncertainty.type == "fixed":
                self._uncertainty[:] = self._options.uncertainty.value

    def get_tfmra_threshold(self, indices):
        """
        Compute the TFMRA threshold for each waveform with several options.
        :param indices: A list of array indices for which to compute thresholds (e.g. for sea ice waveforms)
        :return: An array of thresholds with the dimensions of the l2 data object, but only values for <indices>
        """

        # short link to options
        option = self._options.threshold

        # Init an empty threshold array
        threshold = np.full(self._l2.longitude.shape, np.nan)

        # [Legacy] Threshold is float in settings file
        if type(option) is float:
            threshold[indices] = option
            return threshold

        # Option 1: fixed threshold for all waveforms
        if option.type == "fixed":
            threshold[indices] = option.value

        # Option 2 (deprecated): Threshold as a function of sigma0
        elif option.type == "sigma_func":
            sigma0 = self.get_l1b_parameter("classifier", "sigma0")
            value = np.zeros(sigma0.shape)
            for i, coef in enumerate(option.coef):
                value += coef * sigma0**i
            threshold[indices] = value[indices]

        # Option 3 (deprecated): Threshold as a function of sigma0 and sea ice type
        elif option.type == "sitype_sigma_func":
            sigma0 = self.get_l1b_parameter("classifier", "sigma0")
            sitype = self._l2.sitype
            value = np.zeros(sigma0.shape)
            for i, coef in enumerate(option.coef_fyi):
                value += coef * sigma0**i
            value_myi = np.zeros(sigma0.shape)
            for i, coef in enumerate(option.coef_myi):
                value_myi += coef * sigma0**i
            myi_list = np.where(sitype > 0.5)[0]
            value[myi_list] = value_myi[myi_list]
            threshold[indices] = value[indices]

        # Option 4 (deprecated): Threshold as a function of sigma0 and leading-edge
        #   width 3rd order polynomial fit. Algorithm for CCI/C3S Envisat CDR
        elif option.type == "poly_plane_fit":

            # Get the required classifier
            sigma0 = self.get_l1b_parameter("classifier", "sigma0")
            lew = self.get_l1b_parameter("classifier", "leading_edge_width")
            value = np.zeros(sigma0.shape)
            value += option.intercept
            for i, coef in enumerate(option.coef_lew):
                value += coef * lew**(i+1)
            for i, coef in enumerate(option.coef_sig0):
                value += coef * sigma0**(i+1)
            threshold[indices] = value[indices]

        # TODO: remove dependency of TFMRA threshold to l2 object
        elif option.type == "l2_variable":
            variable_name = option.get("variable_name", None)
            if variable_name is None:
                msg = "Missing option `variable_name` (options.threshold.variable_name) for threshold type `variable`"
                raise ValueError(msg)
            preset_threshold = self._l2.get_parameter_by_name(variable_name)
            threshold[indices] = preset_threshold[indices]

        # Catching invalid processor settings
        else:
            msg = "treshold type not recognized: %s" % str(option.type)
            raise ValueError(msg)

        # All done, return the threshold
        return threshold

    def get_preprocessed_wfm(self, rng, wfm, radar_mode, is_valid):
        """
        Returns the intermediate product (oversampled range bins,
        oversampled and filtered waveforms, indices of first maxima
        and peak power norm for custom applications
        """

        # ---  Get TFMRA options ---
        # The bin range of the waveform where the noise level should be detected
        # NOTE: In the config file the bins values refer to the actual range bins
        #       and not to the oversampled waveforms. This settings depends on
        #       the specific sensor and a default value should not be used.
        #       Thus a ValueError is raised if the option is missing
        noise_level_range_idx = self._options.get("noise_level_range_bin_idx", None)
        if noise_level_range_idx is None:
            raise ValueError("Missing ctfmra option: noise_level_range_bin_idx")

        # Waveforms from some altimeter have artefacts at the beginning of the
        # waveform that may be misclassified as a first maximum.
        # NOTE: This setting also depends on the sensor, but a default can be used.
        fmi_first_valid_idx = self._options.get("first_maximum_ignore_leading_bins", 0)

        # A factor by how much points the waveform should be oversampled
        # before smoothing
        oversampling_factor = self._options.wfm_oversampling_factor

        # The window size for the box filter (radar mode dependant list)
        wfm_smoothing_window_size = self._options.wfm_smoothing_window_size

        # The power threshold for the first maximum (radar mode dependant list)
        first_maximum_normalized_threshold = self._options.first_maximum_normalized_threshold

        # ---  Prepare Output ---

        # Dimensions of the waveforms
        wfm_shape = (wfm.shape[0], wfm.shape[1]*oversampling_factor)

        # Window delay (range)
        filt_rng = np.full(wfm_shape, np.nan)
        # Waveform Power
        filt_wfm = np.full(wfm_shape, np.nan)
        # First Maximum Index
        fmi = np.full(wfm_shape[0], -1, dtype=np.int32)
        # Waveform norm (max power)
        norm = np.full(wfm_shape[0], np.nan)

        # --- Loop over all Waveforms ---
        for i in np.arange(wfm.shape[0]):

            if not is_valid[i]:
                continue

            # Get the filtered waveform
            window_size = wfm_smoothing_window_size[radar_mode[i]]
            result = self.get_filtered_wfm(rng[i, :], wfm[i, :], oversampling_factor, window_size)
            filt_rng[i, :], filt_wfm[i, :], norm[i] = result

            # Get noise level in normalized units
            i0, i1 = [idx * oversampling_factor for idx in noise_level_range_idx]
            noise_level_normed = cytfmra_wfm_noise_level(filt_wfm[i, :], i0, i1)

            # Find first maxima
            # (needs to be above radar mode dependent noise threshold)
            fmnt = first_maximum_normalized_threshold[radar_mode[i]]
            peak_minimum_power = fmnt + noise_level_normed
            fmi_first_valid_idx_filt = fmi_first_valid_idx * oversampling_factor
            fmi[i] = self.get_first_maximum_index(filt_wfm[i, :], peak_minimum_power, fmi_first_valid_idx_filt)

        return filt_rng, filt_wfm, fmi, norm

    def get_thresholds_distance(self, rng, wfm, fmi, t0, t1):
        """
        Return the distance between two thresholds t0 < t1
        """

        width = np.full(rng.shape[0], np.nan, dtype=np.float32)
        noise_level_range_bin_idx = self._options.noise_level_range_bin_idx
        oversampling_factor = self._options.wfm_oversampling_factor
        first_valid_idx = noise_level_range_bin_idx[1] * oversampling_factor
        for i in np.arange(rng.shape[0]):
            if fmi[i] is None:
                continue

            r0, p0, i0 = self.get_threshold_range(rng[i, :], wfm[i, :], fmi[i], t0, first_valid_idx=first_valid_idx)
            r1, p1, i1 = self.get_threshold_range(rng[i, :], wfm[i, :], fmi[i], t1, first_valid_idx=first_valid_idx)
            width[i] = r1 - r0

        # some irregular waveforms might produce negative width values
        is_negative = np.where(width < 0.)[0]
        width[is_negative] = np.nan
        return width

    @staticmethod
    def get_filtered_wfm(rng, wfm, oversampling_factor, window_size):
        """
        Return a filtered version of the the waveform. This process inclused
        oversampling, smoothing and normalization to the first maximum power.
        :param rng: (np.array, dim:n_records) window delay for each range bin
        :param wfm: (np.array, dim:n_records) the power for each range bin with
            sensor dependent units
        :param oversampling_factor: (int) The waveform oversamling factor
        :param window_size: (int) The filter size of the box filter
        :return:
        """

        # Use cython implementation of waveform oversampling
        filt_rng, wfm_os = cytfmra_interpolate(rng.astype(np.float64), wfm.astype(np.float64), oversampling_factor)

        # Smooth the waveform using a box smoother
        filt_wfm = bnsmooth(wfm_os, window_size)

        # Normalize filtered waveform
        filt_wfm, norm = cytfmra_normalize_wfm(filt_wfm)

        # All done, return
        return filt_rng, filt_wfm, norm

    @staticmethod
    def get_first_maximum_index(wfm, peak_minimum_power, first_valid_idx=0):
        """
        Return the index of the first peak (first maximum) on the leading edge
        before the absolute power maximum. The first peak is only valid if
        its power exceeds a certain threshold
        :param wfm: (np.array, dim=(n_range_bins)): Normalied waveform power
        :param peak_minimum_power: (float) threshold for normalized power that
            a peak must be surpass to be regarded as a first maximum candidate
        :param first_valid_idx: (int):
        :return:
        """

        # Get the main maximum first
        try:
            absolute_maximum_index = bn.nanargmax(wfm)
        except ValueError:
            return -1

        # Find relative maxima before the absolute maximum
        # NOTE: The search function uses a subset before the absolute maximum
        #       with an option to ignore the first range bins, which contain
        #       FFT artefacts for some platforms (e.g. ERS-1/2 & Envisat)
        #       The first valid index for the first maximum needs to be
        #       specified in the config file and defaults to 0.
        peaks = cytfmra_findpeaks(wfm[first_valid_idx:absolute_maximum_index+1])+first_valid_idx

        # Check if relative maximum are above the required threshold
        leading_maxima = np.where(wfm[peaks] >= peak_minimum_power)[0]

        # Identify the first maximum
        first_maximum_index = int(absolute_maximum_index)
        if len(leading_maxima) > 0:
            first_maximum_index = peaks[leading_maxima[0]]

        return first_maximum_index

    @staticmethod
    def get_threshold_range(rng: np.ndarray,
                            wfm: np.ndarray,
                            first_maximum_index: int,
                            threshold: float,
                            first_valid_idx: int = 0) -> Tuple[float, float, int]:
        """
        Return the range value and the power of the retrack point at
        a given threshold of the firsts maximum power
        :param rng: (np.array, dim=(n_range_bins) Window delay in meters
        :param wfm: (np.array, dim=(n_range_bins) Waveform power in normalized units
        :param first_maximum_index: (int) Index of first maximum
        :param threshold: (float) Power threshold
        :param first_valid_idx: (int) First valid index for first maximum / leading edge
        :return: tfmra range (float), tfmra power (float)
        """

        # get first index greater as threshold power
        first_maximum_power = wfm[first_maximum_index]

        # Get power of retracked point
        tfmra_power = threshold*first_maximum_power

        # Use linear interpolation to get exact range value
        points = np.where(wfm[first_valid_idx:first_maximum_index] > tfmra_power)[0]+first_valid_idx

        # Check if something went wrong with the first maximum
        if len(points) == 0:
            return np.nan, np.nan, np.nan

        i0, i1 = points[0]-1, points[0]
        gradient = (wfm[i1]-wfm[i0])/(rng[i1]-rng[i0])
        tfmra_range = (tfmra_power - wfm[i0]) / gradient + rng[i0]

        return tfmra_range, tfmra_power, i0


class SICCILead(BaseRetracker):

    def __init__(self):
        super(SICCILead, self).__init__()

    def create_retracker_properties(self, n_records):
        parameter = [
            "retracked_bin",
            "maximum_power_bin",
            "sigma",
            "k",
            "alpha",
            "power_in_echo_tail",
            "rms_echo_and_model"]
        for parameter_name in parameter:
            setattr(self, parameter_name, np.ndarray(shape=n_records, dtype=np.float32) * np.nan)

    def l2_retrack(self, rng, wfm, indices, radar_mode, is_valid):
        # Run the retracker
        self._sicci_lead_retracker(rng, wfm, indices)
        # Filter the results
        if self._options.filter.use_filter:
            self._filter_results()

    def _sicci_lead_retracker(self, range, wfm, indices):

        # retracker options (see l2 settings file)
        skip = self._options.skip_first_bins
        initial_guess = list(self._options.initial_guess)
        maxfev = self._options.maxfev
        time = np.arange(wfm.shape[1]-skip).astype(float)
        x = np.arange(wfm.shape[1])

        # Loop over lead indices
        for index in indices:
            wave = wfm[index, skip:]
            initial_guess[0] = np.argmax(wave)
            initial_guess[3] = np.max(wave)
            popt, cov = curve_fit(P_lead, time, wave.astype(float),
                                  p0=initial_guess, maxfev=maxfev)

            # import matplotlib.pyplot as plt
            #
            # plt.figure(dpi=150)
            # plt.scatter(time, wave, s=1, color="black")
            # time_oversampled = np.linspace(time[0], time[-1], 1000)
            # plt.plot(time_oversampled, P_lead(time_oversampled, *popt), color="red", alpha=0.5)
            # plt.xlim(popt[0]-20, popt[0]+30)
            # plt.show()


            # try:
            #     popt, cov = curve_fit(P_lead, time, wave.astype(float),
            #                           p0=initial_guess, maxfev=maxfev)
            # except:
            #     continue
                # popt = [np.nan, np.nan, np.nan, np.nan]

            # Store retracker parameter for filtering
            # tracking point in units of range bins
            # tracking point in units of range bins
            self.retracked_bin[index] = skip + popt[0]
            self.k[index] = popt[1]
            self.sigma[index] = popt[2]
            self.alpha[index] = popt[3]
            self.maximum_power_bin[index] = np.argmax(wave)

            # Get derived parameter
            self.power_in_echo_tail[index] = power_in_echo_tail(
                wfm[index, :], self.retracked_bin[index], self.alpha[index])
            self.rms_echo_and_model[index] = rms_echo_and_model(
                wfm[index, :], self.retracked_bin[index],
                self.k[index], self.sigma[index], self.alpha[index])

            # Get the range by interpolation of range bin location
            try:
                self._range[index] = interp1d(
                    x, range[index, :], kind='linear', copy=False)(
                        self.retracked_bin[index])
            except ValueError:
                self._range[index] = np.nan
                self._range[index] = np.nan

    def _filter_results(self):
        """ Filter the lead results based on threshold defined in SICCI """

        thrs = self._options.filter
        clf = self._classifier

        valid = ANDCondition()
        valid.add(self.sigma < thrs.maximum_std_of_gaussion_rise)
        valid.add(self.maximum_power_bin > thrs.minimum_bin_count_maxpower)

        # sea ice backscatter not available for ERS?
        if thrs.minimum_echo_backscatter is not None:
            valid.add(clf.sigma0 > thrs.minimum_echo_backscatter)

        bin_seperation = np.abs(self.retracked_bin - self.maximum_power_bin)
        valid.add(bin_seperation < thrs.maximum_retracker_maxpower_binsep)

        valid.add(self.power_in_echo_tail < thrs.maximum_power_in_echo_tail)
        valid.add(self.rms_echo_and_model < thrs.maximum_rms_echo_model_diff)
        valid.add(self.retracked_bin > thrs.sensible_lead_retracked_bin[0])
        valid.add(self.retracked_bin < thrs.sensible_lead_retracked_bin[1])

        # Error flag is also computed for other surface types, do not
        # overide those
        error_flag = self._flag
        error_flag[self.indices] = np.logical_not(valid.flag[self.indices])
        self._flag = error_flag

        # import matplotlib.pyplot as plt
        #
        # f, ax = plt.subplots(6, sharex=True, facecolor="white",  figsize=(10, 16))
        # ax[0].plot(self.retracked_bin[self.indices], lw=0.5, color="#00ace5")
        # ax[0].set_title("retracked_bin")
        # ax[0].axhline(thrs.sensible_lead_retracked_bin[0], color="green")
        # ax[0].axhline(thrs.sensible_lead_retracked_bin[1], color="red")
        #
        # ax[1].plot(self.maximum_power_bin[self.indices],
        #           lw=0.5, color="#00ace5")
        # ax[1].set_title("maximum_power_bin")
        # ax[1].axhline(thrs.minimum_bin_count_maxpower, color="green")
        #
        # ax[2].plot(self.sigma[self.indices], lw=0.5, color="#00ace5")
        # ax[2].set_title("sigma")
        # ax[2].axhline(thrs.maximum_std_of_gaussion_rise, color="red")
        #
        # ax[3].plot(np.abs(self.retracked_bin[self.indices]- self.maximum_power_bin[self.indices]),
        #           lw=0.5, color="#00ace5")
        # ax[3].set_title("retracked bin - max power bin")
        # ax[3].axhline(thrs.minimum_echo_backscatter, color="green")
        #
        # ax[4].plot(self.power_in_echo_tail[self.indices],
        #           lw=0.5, color="#00ace5")
        # ax[4].set_title("power_in_echo_tail")
        # ax[4].axhline(thrs.maximum_power_in_echo_tail, color="red")
        # ax[4].set_ylim(0, 1)
        #
        # ax[5].plot(self.rms_echo_and_model[self.indices],
        #           lw=0.5, color="#00ace5")
        # ax[5].set_title("rms_echo_and_model")
        # ax[5].axhline(thrs.maximum_rms_echo_model_diff, color="red")
        # # ax[5].set_ylim(0, 30)
        #
        # for i in np.arange(6):
        #    ax[i].yaxis.grid(True, which='minor')
        #    ax[i].yaxis.set_tick_params(direction='out')
        #    ax[i].yaxis.set_ticks_position('left')
        #    ax[i].xaxis.set_ticks([])
        #    spines_to_remove = ["top", "right", "bottom"]
        #    for spine in spines_to_remove:
        #        ax[i].spines[spine].set_visible(False)
        #
        # plt.tight_layout()

#
#        plt.figure()
#        plt.plot(valid.flag[self.indices])
#        plt.ylim(-0.1, 1.1)
#
#        plt.show(block=True)
#        # stop


class SICCIOcog(BaseRetracker):

    def __init__(self):
        super(SICCIOcog, self).__init__()
        self.retracked_bin = None
        self.leading_edge_width = None
        self.tail_shape = None

    def create_retracker_properties(self, n_records):
        parameter = ["retracked_bin", "leading_edge_width", "tail_shape"]
        for parameter_name in parameter:
            setattr(self, parameter_name,
                    np.ndarray(shape=(n_records), dtype=np.float32) * np.nan)

    def l2_retrack(self, range, wfm, indices, radar_mode, is_valid):
        # Run the retracker
        self._sicci_ice_retracker(range, wfm, indices)
        # Filter the results
        if self._options.filter.use_filter:
            self._filter_results()

    def _sicci_ice_retracker(self, range, wfm, indices):
        # All retracker options need to be defined in the l2 settings file
        # skip the first bins due to fft noise
        skip = self._options.skip_first_bins

        # percentage of the OCOG retracking points
        percentage = self._options.percentage

        # percentage of earlier retracking point to estimate leading edge width
        lew_percentage = self._options.leading_edge_width_percentage

        # dummy array for interpolation of actual retracker range window
        x = np.arange(wfm.shape[1])

        # Loop over waveform indices marked as leads
        for index in indices:
            wave = np.array(wfm[index, skip:]).astype("float64")
            range_bin = ocog_func(wave, percentage, skip)
            range_bin_lew = ocog_func(wave, lew_percentage, skip)
            try:
                self._range[index] = interp1d(
                    x, range[index, :], kind='linear', copy=False)(range_bin)
            except ValueError:
                self.retracked_bin[index] = np.nan
                self.leading_edge_width[index] = np.nan
                self.tail_shape[index] = np.nan
                continue

            # Store additional retracker parameter
            self.retracked_bin[index] = range_bin
            self.leading_edge_width[index] = range_bin - range_bin_lew
            try:
                self.tail_shape[index] = ocog_tail_shape(
                    wfm[index, :], range_bin)
            except (ValueError, TypeError):
                self.tail_shape[index] = np.nan

    def _filter_results(self):
        """ These threshold are based on the SICCI code"""
        thrs = self._options.filter

        valid = ANDCondition()
        valid.add(self.leading_edge_width < thrs.maximum_leading_edge_width)
        valid.add(self.tail_shape < thrs.maximum_echo_tail_line_deviation)
        valid.add(self.retracked_bin > thrs.sensible_seaice_retracked_bin[0])
        valid.add(self.retracked_bin < thrs.sensible_seaice_retracked_bin[1])

        # Error flag is also computed for other surface types, do not
        # overide those
        error_flag = self._flag
        error_flag[self.indices] = np.logical_not(valid.flag[self.indices])
        self._flag = error_flag

#        import matplotlib.pyplot as plt
#
#        f, ax = plt.subplots(3, sharex=True, facecolor="white",
#                             figsize=(10, 16))
#        ax[0].plot(self.retracked_bin[self.indices], lw=0.5, color="#00ace5")
#        ax[0].set_title("retracked_bin")
#        ax[0].axhline(thrs.sensible_seaice_retracked_bin[0], color="green")
#        ax[0].axhline(thrs.sensible_seaice_retracked_bin[1], color="red")
#
#        ax[1].plot(self.tail_shape[self.indices],
#                   lw=0.5, color="#00ace5")
#        ax[1].set_title("tail_shape")
#        ax[1].axhline(thrs.maximum_echo_tail_line_deviation, color="red")
#
#        ax[2].plot(self.leading_edge_width[self.indices],
#                   lw=0.5, color="#00ace5")
#        ax[2].set_title("leading_edge_width")
#        # ax[2].axhline(thrs.maximum_leading_edge_width, color="red")
#
#        for i in np.arange(3):
#            ax[i].yaxis.grid(True, which='minor')
#            ax[i].yaxis.set_tick_params(direction='out')
#            ax[i].yaxis.set_ticks_position('left')
#            ax[i].xaxis.set_ticks([])
#            spines_to_remove = ["top", "right", "bottom"]
#            for spine in spines_to_remove:
#                ax[i].spines[spine].set_visible(False)
#
#        plt.tight_layout()
#
#        plt.figure()
#        plt.plot(valid.flag[self.indices])
#        plt.ylim(-0.1, 1.1)
#
#        plt.show(block=True)
#        # stop


class ICESatGLAH13Elevation(BaseRetracker):
    """ Pseudo-retracker for ICESat GLAH13 data """

    DOCSTR = r"Pseudo-retracker for ICESat GLAH13 data"

    def __init__(self):
        super(ICESatGLAH13Elevation, self).__init__()

    def create_retracker_properties(self, n_records):
        pass

    def l2_retrack(self, rng, wfm, indices, radar_mode, is_valid):
        """ API Calling method """

        # For l1bdata files from glah13 data, the range that leads to
        # the sea ice surface elevation is encoded in the first entry
        # for record in the waveform range array
        for i in indices:
            # Mandatory return function
            self._range[i] = rng[i, 0] if is_valid[i] else np.nan
            self._power[i] = 1.0

        if "uncertainty" in self._options:
            if self._options.uncertainty.type == "fixed":
                self._uncertainty[:] = self._options.uncertainty.value


class ERSPulseDeblurring(Level2ProcessorStep):
    """
    A processing step applying the pulse deblurring correction of
    Peacock 1998 to retracked ranges.

    NOTE: Should only be used for ERS-1/2 data
    """

    def __init__(self, *args, **kwargs):
        super(ERSPulseDeblurring, self).__init__(*args, **kwargs)

    def execute_procstep(self, l1b, l2):
        """
        Compute the pulse deblurring correction

            Hcor = h + eps / 5. (eps < 0)

        based on the classifier data transferred from the l1p (epss: epsilon in seconds).
        :param l1b:
        :param l2:
        :return: error_status_flag
        """

        # Get a clean error status
        error_status = self.get_clean_error_status(l2.n_records)

        # Compute epsilon in meter (eps_m = eps_sec * c / 2.)
        eps = l2.epss * 0.5 * 299792458.

        # Compute and apply pulse deblurring correction
        pulse_deblurring_correction = np.array(eps < 0.).astype(float) * eps / 5.0
        for target_variable in self.target_variables:
            var = l2.get_parameter_by_name(target_variable)
            var[:] = var[:] + pulse_deblurring_correction
            l2.set_parameter(target_variable, var[:], var.uncertainty[:])

        # Add pulse deblurring correction to level-2 auxiliary data
        l2.set_auxiliary_parameter("pdbc", "pulse_deblurring_correction", pulse_deblurring_correction)

        # Return clean error status (for now)
        return error_status

    @property
    def l2_input_vars(self):
        return ["elev", "epss"]

    @property
    def l2_output_vars(self):
        return ["pdbc"]

    @property
    def target_variables(self):
        return self.cfg.options.get("target_variables", ["elevation"])

    @property
    def error_bit(self):
        return self.error_flag_bit_dict["retracker"]


class SGDRMultipleElevations(Level2ProcessorStep):
    """
    A processing step that computes elevation from a set of
    already retracked ranges and sets these as auxiliary
    variables

    NOTE: This step is only useful if the retracked ranges
          can be taken from a L2/SGDR file

          Also, if any corrections to the elevations are
          to be maded, the elevation parameter names
          have to be added as target_variables in the
          corresponding L2 processor steps
    """

    def __init__(self, *args, **kwargs):
        super(SGDRMultipleElevations, self).__init__(*args, **kwargs)

    def execute_procstep(self, l1b, l2):
        """
        Computes elevations from pre-defined retrackers
        :param l1b:
        :param l2:
        :return: error_status_flag
        """

        # Get a clean error status
        error_status = self.get_clean_error_status(l2.n_records)

        # Get options (and break if incorrect)
        classifier_name_fmt = self.cfg.options.get("classifier_name_fmt", None)
        output_name_fmt = self.cfg.options.get("output_name_fmt", None)
        if None in [classifier_name_fmt, output_name_fmt]:
            logger.error("- Invalid options to {}, skipping".format(self.__class__.__name__))
            error_status[:] = True
            return error_status

        # Get initial elevations (altitude - range)
        # NOTE: It is assumed here that instrument and center of gravity correction
        #       have already been applied to the ranges in the l1 data
        for target_retracker in self.target_retrackers:

            # Get the retracked range from the l1 classifier data group
            classifer_var_name = classifier_name_fmt.format(target_retracker)
            retracker_range = l1b.classifier.get_parameter(classifer_var_name)

            # Compute elevation and add to l2
            elevation_value = l2.altitude[:] - retracker_range
            aux_id = self.auxid_fmt.format(target_retracker)
            aux_name = output_name_fmt.format(target_retracker)
            l2.set_auxiliary_parameter(aux_id, aux_name, elevation_value)

        # Return clean error status (for now)
        return error_status

    @property
    def target_retrackers(self):
        return self.cfg.options.get("predefined_retrackers", [])

    @property
    def l2_input_vars(self):
        return []

    @property
    def l2_output_vars(self):
        return [self.auxid_fmt.format(retracker_id for retracker_id in self.target_retrackers)]

    @property
    def auxid_fmt(self):
        return "elev{}"

    @property
    def error_bit(self):
        return self.error_flag_bit_dict["other"]


class TFMRAMultiThresholdFreeboards(Level2ProcessorStep):
    """
    Level-2 processor step to compute elevations from a range of TFMRA thresholds.
    NOTE: Computational expensive, should be handled with care
    """

    def __init__(self, *args, **kwargs):
        super(TFMRAMultiThresholdFreeboards, self).__init__(*args, **kwargs)

    def execute_procstep(self, l1b, l2):
        """
        Computes the elevation of a range of defined of TFMRA retrackers thresholds to all
        lead or ice elevations.
        :param l1b:
        :param l2:
        :return: error_status_flag
        """

        # Get a clean error status
        error_status = self.get_clean_error_status(l2.n_records)

        # Get the retracker and retrack only sea-ice surfaces
        tfmra = cTFMRA()
        tfmra.set_default_options()
        filt_rng, filt_wfm, fmi, norm = tfmra.get_preprocessed_wfm(
            l1b.waveform.range,
            l1b.waveform.power,
            l1b.waveform.radar_mode,
            l2.surface_type.sea_ice.flag)

        # Get the target thresholds
        threshold_min, threshold_max = self.cfg.options.threshold_range
        threshold_res = self.cfg.options.threshold_res
        thresholds = np.arange(threshold_min, threshold_max+0.1*threshold_res, threshold_res)

        # Construct the output array
        ranges = np.full((l2.n_records, thresholds.size), np.nan)

        for i in np.where(l2.surface_type.sea_ice.flag)[0]:
            for j, threshold in enumerate(thresholds):
                retracked_range, _, _ = tfmra.get_threshold_range(
                    filt_rng[i, :],
                    filt_wfm[i, :],
                    fmi[i],
                    threshold)
                ranges[i, j] = retracked_range

        # Convert ranges to freeboards using the already pre-computed steps.
        # NOTE: This specifically assumes that the sea surface height is constant
        freeboards = np.full(ranges.shape, np.nan)
        for j in np.arange(thresholds.size):
            freeboards[:, j] = l2.altitude - ranges[:, j] - l2.rctotal - l2.mss - l2.sla + l2.sgcor

        # Register results as auxiliary data variable
        dim_dict = {"new_dims": (("tfmra_thresholds", thresholds.size), ),
                    "dimensions": ("time", "tfmra_thresholds"),
                    "add_dims": (("tfmra_thresholds", thresholds), )}
        l2.set_multidim_auxiliary_parameter("thfrbs", "threshold_freeboards", freeboards, dim_dict)

        # Return clean error status (for now)
        return error_status

    @property
    def l2_input_vars(self):
        return ["range"]

    @property
    def l2_output_vars(self):
        return ["tfmra_elevation_range"]

    @property
    def auxid_fmt(self):
        return "elev{}"

    @property
    def error_bit(self):
        return self.error_flag_bit_dict["other"]

# %% Function for CryoSat-2 based retracker


def wfm_get_noise_level(wfm, oversample_factor):
    """ According to CS2AWI TFMRA implementation """
    return bn.nanmean(wfm[:5*oversample_factor])


def smooth(x, window):
    """ Numpy implementation of the IDL SMOOTH function """
    return np.convolve(x, np.ones(window)/window, mode='same')


def bnsmooth(x, window):
    """ Bottleneck implementation of the IDL SMOOTH function """
    pad = int((window-1)/2)
    n = len(x)
    xpad = np.ndarray(shape=(n+window))
    xpad[0:pad] = 0.0
    xpad[pad:n+pad] = x
    xpad[n+pad:] = 0.0
    return bn.move_mean(xpad, window=window, axis=0)[window-1:(window+n-1)]


def peakdet(v, delta, x=None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html

    Returns two arrays

    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %
    %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %        in MAXTAB and MINTAB are replaced with the corresponding
    %        X-values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.

    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.

    """
    maxtab = []
    mintab = []

    if x is None:
        x = np.arange(len(v))

    v = np.asarray(v)

    mn, mx = np.Inf, -np.Inf
    mnpos, mxpos = np.NaN, np.NaN

    lookformax = True

    for i in np.arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]

        if lookformax:
            if this < mx-delta:
                maxtab.append(mxpos)
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append(mnpos)
                mx = this
                mxpos = x[i]
                lookformax = True

    return np.array(maxtab), np.array(mintab)


def findpeaks(data, spacing=1, limit=None):
    """Finds peaks in `data` which are of `spacing` width and >=`limit`.
    :param data: values
    :param spacing: minimum spacing to the next peak (should be 1 or more)
    :param limit: peaks should have value greater or equal
    :return:
    """
    len = data.size
    x = np.zeros(len+2*spacing)
    x[:spacing] = data[0]-1.e-6
    x[-spacing:] = data[-1]-1.e-6
    x[spacing:spacing+len] = data
    peak_candidate = np.zeros(len)
    peak_candidate[:] = True
    for s in range(spacing):
        start = spacing - s - 1
        h_b = x[start:start + len]  # before
        start = spacing
        h_c = x[start:start + len]  # central
        start = spacing + s + 1
        h_a = x[start:start + len]  # after
        peak_candidate = np.logical_and(
            peak_candidate, np.logical_and(h_c > h_b, h_c > h_a))

    ind = np.argwhere(peak_candidate)
    ind = ind.reshape(ind.size)
    if limit is not None:
        ind = ind[data[ind] > limit]
    return ind


# %% Functions for SICCI retracker

def ocog_func(wave, percentage, skip):
    waveform = wave*wave
    sq_sum = np.sum(waveform)
    waveform = waveform*waveform
    qa_sum = np.sum(waveform)
    # Calculate retracking threshold (2.6.2 in ATDBv0)
    threshold = percentage * np.sqrt(qa_sum / sq_sum)
    ind_first_over = np.min(np.nonzero(wave > threshold))
    decimal = (wave[ind_first_over-1] - threshold) / \
        (wave[ind_first_over-1] - wave[ind_first_over])
    x_range_bin = skip + ind_first_over - 1 + decimal
    return x_range_bin


def P_lead(t, t_0, k, sigma, a):
    """ Lead waveform model (for SICCILead curve fitting) """
    # Time for F to be F_L
    t_b = k*sigma**2
    # Helper coefficient
    sq_ktb = np.sqrt(k*t_b)
    # Polynomial coefficients for F_L
    aa = ((5*k*sigma) - (4*sq_ktb)) / (2*sigma*t_b*sq_ktb)
    aaa = ((2*sq_ktb)-(3*k*sigma)) / (2*sigma*t_b*t_b*sq_ktb)
    # We are mostly interested in time in reference to t_0
    t_diff = (t - t_0)
    # F is piecewise. These are (Truth where applicable) * (Function)
    F_1 = (t <= t_0) * (t_diff/sigma)
    F_2 = (t >= (t_b+t_0)) * (np.sqrt(np.abs(t_diff)*k))
    F_L = np.logical_and(t > t_0, t < (t_b+t_0)) * (aaa*t_diff**3 + aa*t_diff**2+t_diff/sigma)
    # Compile F
    F = F_1 + F_L + F_2
    return a*np.exp(-F*F)  # Return e^-f^2(t)


def power_in_echo_tail(wfm, retracked_bin, alpha, pad=3):
    """
    The tail power is computed by summing the bin count in all bins from
    2 bins beyond the tracking point to the end of the range gate, then
    normalised by dividing by the value of alpha returned from the lead
    retracking
    source: SICCI
    """
    try:
        tracking_point = int(retracked_bin)
        return sum(wfm[tracking_point+pad:])/alpha
    except:
        return np.nan


def ocog_tail_shape(wfm, tracking_point, tail_pad=3):
    """ From SICCI module """
    tail = wfm[tracking_point+tail_pad:]
    tail = tail / np.mean(tail)
    P = np.polyfit(np.arange(len(tail)), tail, 1)
    residual = (tail-np.polyval(P, np.arange(len(tail))))
    return np.sqrt(np.sum(residual*residual) / len(residual))


def rms_echo_and_model(wfm, retracked_bin, k, sigma, alpha):
    """
    The root sum squared difference between the echo and the fitted function
    in the lead retracking is computed. The 5 bins before the tracking point
    are used as the echo rise.
    """
    try:
        tracking_point = int(retracked_bin)
        time = np.arange(len(wfm)).astype(float)
        modelled_wave = P_lead(time, retracked_bin, k, sigma, alpha)
        diff = wfm[tracking_point-4:tracking_point+1] - \
            modelled_wave[tracking_point-4:tracking_point+1]
        return np.sqrt(np.sum(diff*diff)/5)/alpha
    except:
        return np.nan

class SAMOSAPlus(BaseRetracker):
    """
    Interface to the SAMOSA+ retracker by CLS.
    Retracker must be installed as a package into the environment for this class to be used.

    """

    def __init__(self):
        super(SAMOSAPlus, self).__init__()

    def create_retracker_properties(self, n_records):
        # False branches here and below were used for debugging
        if True:
            parameter = ["misfit", "swh", "wind_speed", "oceanlike_flag", "epoch", "guess", "Pu", "rval", "kval", "pval", "cval"]
        else:
            parameter = ["misfit", "swh", "wind_speed", "oceanlike_flag"]
        for parameter_name in parameter:
            setattr(self, parameter_name,
                    np.ndarray(shape=(n_records), dtype=np.float32) * np.nan)

    def l2_retrack(self, range, wfm, indices, radar_mode, is_valid):
        # Run the retracker
        self._samosa_plus_retracker(range, wfm, indices, radar_mode)
        # Filter the results
        # Needs a filter option in the config file
        #if self._options.filter.use_filter:
        #    self._filter_results()

    def _samosa_plus_retracker(self, range, wfm, indices, radar_mode):
        # Range contains the range to each bin in each waveform

        # All retracker options need to be defined in the l2 settings file
        # How to get an option
        # opt_val = self._options.name_in_file
        ### CST is a structure collecting universal constants

        CST = type('', (), {})()

        CST.c0 = 299792458.  ## speed of light in m/sec
        CST.R_e = 6378137.  ## Reference Ellipsoid Earh Radius in m
        CST.f_e = 1 / 298.257223563  ## Reference Ellipsoid Earth Flatness
        CST.gamma_3_4 = 1.2254167024651779  ## Gamma Function Value at 3/4

        ### OPT is a structure collecting parameters relative to the minimization scheme settings

        OPT = type('', (), {})()

        OPT.method = 'trf'  ## acronym of the minimization solver, see scipy.optimize.least_squares for details
        OPT.ftol = 1e-2  ## exit tolerance on f
        OPT.gtol = 1e-2  ## exit tolerance on gradient norm of f
        OPT.xtol = 2 * 1e-3  ## exit tolerance on x
        OPT.diff_step = None  ## relative step size for the finite difference approximation of the Jacobian
        OPT.max_nfev = None  ## maximum number of function evaluations
        OPT.loss = 'linear'  ## loss function , see scipy.optimize.least_squares for details

        ### RDB is a structure collecting parameters relative to the sensor radar database

        RDB = type('', (), {})()

        RDB.Np_burst = 64  # number of pulses per burst
        RDB.Npulse = 128  # number of the range gates per pulse (without zero-padding)
        RDB.PRF_SAR = 18181.8181818181  # Pulse Repetition Frequency in SAR mode , given in Hz
        RDB.BRI = 0.0117929625  # Burst Repetition Interval, given in sec
        RDB.f_0 = 13.575e9  # Carrier Frequency in Hz
        RDB.Bs = 320e6  # Sampled Bandwidth in Hz
        RDB.theta_3x = np.deg2rad(1.10)  # (rad) Antenna 3 dB beamwidth (along-track)
        RDB.theta_3y = np.deg2rad(1.22)  # (rad) Antenna 3 dB beamwidth (cross-track)
        RDB.G_0 = 10. ** (42.6 / 10)  # Boresight One-Way Antenna Power Radiation Gain (natural units)
        RDB.bias_sigma0 = -3.04  # static bias in sigma0 (dB)

        ### LUT is a structure collecting filenames of the all SAMOSA LUT
        ### All the LUT files must be in a folder named auxi and located in the same folder as sampy.py

        # FIXME: Need to add a configuration parameter that gives a path to these files. Per surface type. Possibly
        # per mode as well?
        LUT = type('', (), {})()

        LUT.F0 = 'LUT_F0.txt'  ## filename of the F0 LUT
        LUT.F1 = 'LUT_F1.txt'  ## filename of the F1 LUT
        LUT.alphap_noweight = 'alphap_table_DX3000_ZP20_SWH20_10_Sept_2019(CS2_NOHAMMING).txt'  ## filename of the alphap LUT ( case no weighting)
        LUT.alphap_weight = 'alphap_table_DX3000_ZP20_SWH20_10_Sept_2019(CS2_HAMMING).txt'  ## filename of the alphap LUT ( case weighting)
        LUT.alphapower_noweight = 'alphaPower_table_CONSTANT_SWH20_10_Feb_2020(CS2_NOHAMMING).txt'  ## filename of the alpha power LUT ( case no weighting)
        LUT.alphapower_weight = 'alphaPower_table_CONSTANT_SWH20_10_Feb_2020(CS2_NOHAMMING).txt'  ## filename of the alpha power LUT ( case weighting)

        ### time array tau : it gives the relative time of each range gate of the radar waveform with respect a time zero
        ### time zero corresponds at the time of the reference gate

        wf_zp = np.shape(wfm)[1] / RDB.Npulse  #### zero-padding factor of the waveform
        logger.info('Waveform zero padding factor is {:f}'.format(wf_zp))
        Nstart = RDB.Npulse * wf_zp
        Nend = RDB.Npulse * wf_zp
        dt = 1. / (RDB.Bs * wf_zp)  #### time sampling step for the array tau, it includes the zero-padding factor

        # print(np.arange(-(4 / 2), ((4 - 1) / 2)))
        # [-2. -1.  0.  1.]
        # Zero bin as at len()//2. So use the range at this location as the base to adjust from
        tau = np.arange(-(Nstart / 2) * dt, ((Nend - 1) / 2) * dt, dt)

        NstartNoise = 2  ## noise range gate counting from 1, no oversampling
        NendNoise = 6  ## noise range gate counting from 1, no oversampling

        #window_del_20_hr_ku_deuso = window_del_20_hr_ku * (uso_cor_20_hr_ku + 1)
        window_del_20_hr_ku_deuso = self._l1b.classifier.window_delay
        #Raw_Elevation = alt_20_hr_ku - CST.c0 / 2 * window_del_20_hr_ku_deuso
        raw_range = CST.c0 * window_del_20_hr_ku_deuso * 0.5
        Raw_Elevation = self._l1b.time_orbit.altitude - raw_range

        sin_index = np.where(radar_mode == 2)[0]
        if len(sin_index > 0):
            Raw_Elevation[sin_index] = self._l1b.time_orbit.altitude[sin_index] - range[sin_index,np.shape(wfm)[1]//2]
            raw_range[sin_index] = range[sin_index,np.shape(wfm)[1]//2]
            logger.info('Using L1b range array for {:d} SARIn mode records'.format(len(sin_index)))

        ThNEcho = compute_ThNEcho(wfm.T, NstartNoise * wf_zp,
                                  NendNoise * wf_zp)  ### computing Thermal Noise from the waveform matric

        # initialize_epoch relies on the waveform being in counts, not watts, so revert
        wf_norm = np.zeros_like(wfm)
        for rec in np.arange(np.shape(wfm)[0]):
            wf_norm[rec,:] = (65535.0 * wfm[rec,:] / np.nanmax(wfm[rec,:])).round().astype(np.uint16)

        epoch0 = initialize_epoch(wf_norm.T, tau, Raw_Elevation, CST,
                                  size_half_block=10)  ### initializing the epoch (first-guess epoch) from the waveform matrix

        samlib = initialize_SAMOSAlib(CST, RDB, OPT,
                                      LUT)  #### initializing the SAMOSA library sampy, it's a mandatory step

        n = np.shape(wfm)[0]

        GEO = type('', (), {})()

        CONF = type('', (), {})()

        CONF.flag_slope = False                    ### flag True commands to include in the model the slope of orbit and surface (this effect usually is included in LookAngles Array)
        CONF.beamsamp_factor = 1                   ### 1 means only one beam per resolution cell is generated in the DDM, the other ones are decimated
        CONF.wf_weighted = True                   ### flag True if the waveform under iteration is weighted
        CONF.N_Look_min = -90                      ### number of the first Look to generate in the DDM (only used if LookAngles array is not passed in input: i.e. set to  None)
        CONF.N_Look_max = 90                       ### number of the last Look to generate in the DDM (only used if LookAngles array is not passed in input: i.e. set to  None)
        CONF.guess_swh = 2                         ### first-guess SWH in meter
        CONF.guess_pu = 1                          ### first-guess Pu
        CONF.guess_nu = 2                          ### first-guess nu (only used in second step of SAMOSA+)
        CONF.lb_epoch = None                       ### lower bound on epoch in sec. If set to None, lower bound will be set to the first time in input array tau
        CONF.lb_swh = -0.5                         ### lower bound on SWH in m
        CONF.lb_pu = 0.2                           ### lower bound on Pu
        CONF.lb_nu = 0                             ### lower bound on nu (only used in second step of SAMOSA+)
        CONF.ub_epoch = None                       ### upper bound on epoch in sec. If set to None, upper bound will be set to the last time in input array tau
        CONF.ub_swh = 30                           ### upper bound on SWH in m
        CONF.ub_pu = 1.5                           ### upper bound on Pu
        CONF.ub_nu = 1e9                           ### upper bound on nu (only used in second step of SAMOSA+)
        CONF.rtk_type = 'samosa+'                  ### choose between 'samosa' or 'samosa+'
        CONF.wght_factor= 1.4705                   ### widening factor of PTR main lobe after Weighting Window Application
        CONF.Lz = CST.c0 / (2. * RDB.Bs)

        # dummy array for interpolation of actual retracker range window
        x = np.arange(wfm.shape[1])

        # Make an altitude rate. Currently zero as L1b rate is nan. Not used anyway with flag_slope false.
        hrate = np.zeros_like(self._l1b.time_orbit.altitude_rate)

        vel = np.sqrt(self.get_l1b_parameter("classifier", "satellite_velocity_x")**2
                      + self.get_l1b_parameter("classifier", "satellite_velocity_y")**2
                      + self.get_l1b_parameter("classifier", "satellite_velocity_z")**2)

        # Loop over waveform indices marked as surface type
        for index in indices:

            LookAngles = 90 - np.linspace(np.rad2deg(self._l1b.classifier.look_angle_start[index]),
                                          np.rad2deg(self._l1b.classifier.look_angle_stop[index]),
                                          num=int(self._l1b.classifier.stack_beams[index]), endpoint=True)
            MaskRanges = None
            GEO.LAT=self._l1b.time_orbit.latitude[index]                              ### latitude in degree for the waveform under iteration
            GEO.LON=self._l1b.time_orbit.longitude[index]                              ### longitude in degree between -180, 180 for the waveform under iteration
            GEO.Height=self._l1b.time_orbit.altitude[index]                            ### Orbit Height in meter for the waveform under iteration
            GEO.Vs=np.squeeze(vel)[index]                       ### Satellite Velocity in m/s for the waveform under iteration
            GEO.Hrate=hrate[index]                   ### Orbit Height rate in m/s for the waveform under iteration
            GEO.Pitch=np.radians(self._l1b.time_orbit.antenna_pitch[index])     ### Altimeter Reference Frame Pitch in radiant
            GEO.Roll=np.radians(self._l1b.time_orbit.antenna_roll[index])       ### Altimeter Reference Frame Roll in radiant
            GEO.nu=0                                                         ### Inverse of the mean square slope
            GEO.track_sign=0                                                 ### if Track Ascending => -1, if Track Descending => +1, set it to zero if flag_slope=False in CONF
            GEO.ThN=np.squeeze(ThNEcho)[index]                                   ### Thermal Noise

            wf = np.array(wfm[index, :]).astype("float64")

            CONF.guess_epoch = epoch0[index]

            # Do retrack in units of Watts as we don't have the scaling factors available for calc_sigma0
            epoch_sec,swh,Pu,misfit,oceanlike_flag=samlib.Retrack_Samosa(tau,wf,LookAngles,MaskRanges,GEO,CONF)

            # SAMOSA returns a dR based upon the retracker chosen bin sampled from tau
            self._range[index] = raw_range[index] + epoch_sec * CST.c0 * 0.5
            sigma0,pval,cval,rval,kval = calc_sigma0(None, Pu, CST, RDB, GEO, epoch_sec, window_del_20_hr_ku_deuso[index],
                                 GEO.LAT, GEO.Height, GEO.Vs,
                                 self._l1b.classifier.transmit_power[index])
            wind_speed = func_wind_speed([sigma0])

            self._power[index] = sigma0

            # Store additional retracker parameters
            self.swh[index] = swh
            self.misfit[index] = misfit
            self.wind_speed[index] = wind_speed
            self.oceanlike_flag[index] = oceanlike_flag
            if True:
                self.epoch[index] = epoch_sec
                self.guess[index] = epoch0[index]
                self.Pu[index] = 65535.0 * Pu/np.max(wf)
                self.pval[index] = pval
                self.cval[index] = cval
                self.rval[index] = rval
                self.kval[index] = kval

        self.register_auxdata_output("samswh", "samosa_swh", self.swh)
        self.register_auxdata_output("samwsp", "samosa_wind_speed", self.wind_speed)

        if "range_bias" in self._options:
            for radar_mode_index in np.arange(3):
                indices = np.where(radar_mode == radar_mode_index)[0]
                if len(indices) == 0:
                    continue
                range_bias = self._options.range_bias[radar_mode_index]
                self._range[indices] -= range_bias

        if "uncertainty" in self._options:
            if self._options.uncertainty.type == "fixed":
                self._uncertainty[:] = self._options.uncertainty.value

        if True:
            import xarray as xr
            outds = xr.Dataset({'SSHunc': (['time_20_ku'], self._l1b.time_orbit.altitude-self._range),
                                'raw_elev' : (['time_20_ku'], Raw_Elevation),
                                'range': (['time_20_ku'], self._range),
                                'wd_range': (['time_20_ku'], raw_range),
                                'epoch': (['time_20_ku'], self.epoch),
                                'guess': (['time_20_ku'], self.guess),
                                'Pu': (['time_20_ku'], self.Pu),
                                'lat': (['time_20_ku'], self._l1b.time_orbit.latitude),
                                'height': (['time_20_ku'], self._l1b.time_orbit.altitude),
                                'vel': (['time_20_ku'], vel),
                                'pval': (['time_20_ku'], self.pval),
                                'cval': (['time_20_ku'], self.cval),
                                'rval': (['time_20_ku'], self.rval),
                                'kval': (['time_20_ku'], self.kval),
                                'wf': (['time_20_ku','bins'], wf_norm),
                                'misfit': (['time_20_ku'], self.misfit),
                                'oceanlike_flag': (['time_20_ku'], self.oceanlike_flag),
                                'SWH': (['time_20_ku'], self.swh),
                                'sigma0': (['Sigma0_20Hz'], self._power),
                                'wind_speed': (['U10_20Hz'], self.wind_speed)},
                               coords={'time_20_ku': self._l1b.time_orbit.timestamp,
                                       'bins': np.arange(256),
                                       'lon_20_ku': (['time_20_ku'], self._l1b.time_orbit.longitude),
                                       'lat_20_ku': (['time_20_ku'], self._l1b.time_orbit.latitude)},
                               attrs={'description': "Parameters from SAMOSA+ retracker"})
            outds.time_20_ku.encoding = {'calendar': 'gregorian',
                                         'units': 'seconds since 2000-01-01 0:0:0'}
            if os.path.isfile('samosa_debug.nc'):
                os.remove('samosa_debug.nc')
            outds.to_netcdf('samosa_debug.nc')


    def _filter_results(self):
        pass
        """ These threshold are based on the SICCI code"""
        #thrs = self._options.filter

        #valid = ANDCondition()
        #valid.add(self.leading_edge_width < thrs.maximum_leading_edge_width)

        # Error flag is also computed for other surface types, do not
        # overide those
        #error_flag = self._flag
        #error_flag[self.indices] = np.logical_not(valid.flag[self.indices])
        #self._flag = error_flag

class SSBCorrectionJason2(Level2ProcessorStep):
    """
    A processing step applying a sea state bias correction.

    NOTE: Designed for use with SAMOSA retracker and taken from the Jason 2 mission

    Sea state bias model by Tran et al. (2012), CLS-CNES [1]

    [1]Tran N., S. Philipps, J.-C. Poisson, S. Urien, E. Bronner, N. Picot, "Impact of GDR_D standards on SSB
    corrections", Presentation OSTST2012 in Venice, http://www.aviso.altimetry.fr/fileadmin/
    documents/OSTST/2012/oral/02_friday_28/01_instr_processing_I/01_IP1_Tran.pdf

    """

    def __init__(self, *args, **kwargs):
        super(SSBCorrectionJason2, self).__init__(*args, **kwargs)

    def execute_procstep(self, l1b, l2):
        """
        Compute the ssb correction.
        :param l1b:
        :param l2:
        :return: error_status_flag
        """

        # Get a clean error status
        error_status = self.get_clean_error_status(l2.n_records)

        # Compute ssb
        try:
            ssb = calc_ssb_j2Table(l2.samswh, l2.samwsp)
            logger.info('Computing sea state bias correction')
        except AttributeError:
            # This happens if no SAMOSA retracking because no ocean or lead records
            ssb = np.zeros(l2.n_records)

        # Change NaN values for correction to zero
        ssb = np.nan_to_num(ssb, copy=False)

        for target_variable in self.target_variables:
            var = l2.get_parameter_by_name(target_variable)
            var[:] = var[:] - ssb
            l2.set_parameter(target_variable, var[:], var.uncertainty[:])

        # Add pulse deblurring correction to level-2 auxiliary data
        l2.set_auxiliary_parameter("ssb", "sea_state_bias", ssb)

        # Return clean error status (for now)
        return error_status

    @property
    def l2_input_vars(self):
        return ["elev", "samswh", "samwsp"]

    @property
    def l2_output_vars(self):
        return ["ssb"]

    @property
    def target_variables(self):
        return self.cfg.options.get("target_variables", ["elevation"])

    @property
    def error_bit(self):
        return self.error_flag_bit_dict["retracker"]

# %% Retracker getter funtion


def get_retracker_class(name):

    if name == "cTFMRA" and not CYTFMRA_OK:
        msg = "pysiral error: cTFMRA selected but pysiral.bnfunc.cytfmra " + \
              "not available locally\n"
        msg += "(See documentation on compilation with cython)"
        sys.exit(msg)

    if name == "SAMOSAPlus" and not SAMOSA_OK:
        msg = "pysiral error: SAMOSAPlus selected but" + \
              "not available locally\n"
        msg += "(Obtain SAMOSA and install via pip first)"
        sys.exit(msg)

    return globals()[name]()
