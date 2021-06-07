# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 15:48:58 2015

@author: Stefan Hendrics
"""

# Utility methods for retracker:
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import bottleneck as bn

import sys
import time
import numpy as np
from loguru import logger
from attrdict import AttrDict

# cythonized bottleneck functions for cTFMRA
try:
    from .bnfunc.cytfmra import (cytfmra_findpeaks, cytfmra_interpolate,
                                 cytfmra_wfm_noise_level, cytfmra_normalize_wfm)
    CYTFMRA_OK = True
except ImportError:
    logger.error("Cannot import cytfmra")
    CYTFMRA_OK = False

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
                logger.info("- no waveforms of type %s" % surface_type)
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
        default_options_dict = {
            "threshold": dict(type="fixed", value=0.5),
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

        if "uncertainty" in self._options:
            if self._options.uncertainty.type == "fixed":
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
            fmi = self.get_first_maximum_index(filt_wfm, peak_minimum_power, fmi_first_valid_idx)
            tfmra_first_maximum_index[i] = fmi

            # first maximum finder might have failed
            if fmi == -1:
                self._range[i] = np.nan
                self._power[i] = np.nan
                continue

            # Get track point and its power
            tfmra_range, tfmra_power = self.get_threshold_range(filt_rng, filt_wfm, fmi, tfmra_threshold[i])

            # waveform_index = "{:06g}".format(i)
            # treshold_str = "{:03g}".format(tfmra_threshold[i]*100)
            # surface_type = {"095": "Lead", "050": "Ice"}
            # title = "{} Waveform (index: {})".format(surface_type[treshold_str], waveform_index)
            # import matplotlib.pyplot as plt
            # from pathlib import Path
            # from matplotlib.ticker import MultipleLocator
            # plt.figure(figsize=(10, 10))
            # plt.title(title)
            # ax = plt.gca()
            # plt.scatter(rng[i, :]-rng[i, 0], wfm[i, :]/np.nanmax(wfm[i, :]), marker="s", color="0.0",
            #             label="raw waveform")
            # plt.plot(filt_rng-filt_rng[0], filt_wfm, lw=2, label="filtered waveform")
            # label = "tracking point (threshold: {})".format(tfmra_threshold[i])
            # retracked_range = tfmra_range - rng[i, 0]
            # plt.axvline(retracked_range, color="red", lw=2, label=label)
            # plt.legend()
            # plt.xlim(retracked_range-2, retracked_range+10)
            # ax.xaxis.set_major_locator(MultipleLocator(1))
            # ax.xaxis.set_minor_locator(MultipleLocator(0.2))
            # plt.grid(which="major", axis="x", color="0.75", lw=0.2)
            # plt.grid(which="minor", axis="x", color="0.25", lw=0.1)
            # plt.xlabel("Range Window (meter)")
            # plt.ylabel("Normalized Waveform Power")
            # output_path = Path(r"D:\temp\ers_tfmra_test\filter_width_11")
            # output_filename = "tfmra_waveform_{}_{}.png".format(treshold_str, waveform_index)
            # plt.savefig(output_path / output_filename)

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
            fmi[i] = self.get_first_maximum_index(filt_wfm[i, :], peak_minimum_power, fmi_first_valid_idx)

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
        filt_rng, wfm_os = cytfmra_interpolate(rng.astype(np.float32), wfm.astype(np.float32), oversampling_factor)

        # Smooth the waveform using a box smoother
        filt_wfm = bnsmooth(wfm_os, window_size)

        # Normalize filtered waveform
        filt_wfm, norm = cytfmra_normalize_wfm(filt_wfm)

        # All done, return
        return filt_rng, filt_wfm, norm

    @staticmethod
    def get_first_maximum_index(wfm, peak_minimum_power, first_valid_idx):
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
        # but exclude the first range bins for Artifacts in ERS/Envisat
        skip = 55
        # points = np.where(wfm[:first_maximum_index] > tfmra_power)[0]
        points = np.where(wfm[skip:first_maximum_index] > tfmra_power)[0]+skip

        if len(points) == 0:
            return np.nan, np.nan

        i0, i1 = points[0]-1, points[0]
        gradient = (wfm[i1]-wfm[i0])/(rng[i1]-rng[i0])
        tfmra_range = (tfmra_power - wfm[i0]) / gradient + rng[i0]

        return tfmra_range, tfmra_power


class SICCILead(BaseRetracker):

    def __init__(self):
        super(SICCILead, self).__init__()

    def create_retracker_properties(self, n_records):
        parameter = ["retracked_bin", "maximum_power_bin", "sigma", "k",
                     "alpha", "power_in_echo_tail", "rms_echo_and_model"]
        for parameter_name in parameter:
            setattr(self, parameter_name, np.ndarray(shape=n_records, dtype=np.float32) * np.nan)

    def l2_retrack(self, range, wfm, indices, radar_mode, is_valid):
        # Run the retracker
        self._sicci_lead_retracker(range, wfm, indices)
        # Filter the results
        if self._options.filter.use_filter:
            self._filter_results()

    def _sicci_lead_retracker(self, range, wfm, indices):

        # retracker options (see l2 settings file)
        skip = self._options.skip_first_bins
        initial_guess = self._options.initial_guess
        maxfev = self._options.maxfev
        time = np.arange(wfm.shape[1]-skip).astype(float)
        x = np.arange(wfm.shape[1])

        # Loop over lead indices
        for index in indices:
            wave = wfm[index, skip:]
            initial_guess[3] = np.max(wave)
            try:
                popt, cov = curve_fit(P_lead, time, wave.astype(float),
                                      p0=initial_guess, maxfev=maxfev)
            except:
                continue
                # popt = [np.nan, np.nan, np.nan, np.nan]

            # Store retracker parameter for filtering
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

#        import matplotlib.pyplot as plt
#
#        f, ax = plt.subplots(6, sharex=True, facecolor="white",
#                             figsize=(10, 16))
#        ax[0].plot(self.retracked_bin[self.indices], lw=0.5, color="#00ace5")
#        ax[0].set_title("retracked_bin")
#        ax[0].axhline(thrs.sensible_lead_retracked_bin[0], color="green")
#        ax[0].axhline(thrs.sensible_lead_retracked_bin[1], color="red")
#
#        ax[1].plot(self.maximum_power_bin[self.indices],
#                   lw=0.5, color="#00ace5")
#        ax[1].set_title("maximum_power_bin")
#        ax[1].axhline(thrs.minimum_bin_count_maxpower, color="green")
#
#        ax[2].plot(self.sigma[self.indices], lw=0.5, color="#00ace5")
#        ax[2].set_title("sigma")
#        ax[2].axhline(thrs.maximum_std_of_gaussion_rise, color="red")
#
#        ax[3].plot(np.abs(self.retracked_bin[self.indices]- self.maximum_power_bin[self.indices]),
#                   lw=0.5, color="#00ace5")
#        ax[3].set_title("retracked bin - max power bin")
#        ax[3].axhline(thrs.minimum_echo_backscatter, color="green")
#
#        ax[4].plot(self.power_in_echo_tail[self.indices],
#                   lw=0.5, color="#00ace5")
#        ax[4].set_title("power_in_echo_tail")
#        ax[4].axhline(thrs.maximum_power_in_echo_tail, color="red")
#        ax[4].set_ylim(0, 1)
#
#        ax[5].plot(self.rms_echo_and_model[self.indices],
#                   lw=0.5, color="#00ace5")
#        ax[5].set_title("rms_echo_and_model")
#        ax[5].axhline(thrs.maximum_rms_echo_model_diff, color="red")
#        # ax[5].set_ylim(0, 30)
#
#        for i in np.arange(6):
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
            except ValueError:
                self.tail_shape[index] = np.nan
            except TypeError:
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
        filt_rng, filt_wfm, fmi, norm = tfmra.get_preprocessed_wfm(l1b.waveform.range,
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
                retracked_range, _ = tfmra.get_threshold_range(filt_rng[i, :], filt_wfm[i, :], fmi[i], threshold)
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
    return bn.nanmean(wfm[0:5*oversample_factor])


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

# %% Retracker getter funtion


def get_retracker_class(name):

    if name == "cTFMRA" and not CYTFMRA_OK:
        msg = "pysiral error: cTFMRA selected but pysiral.bnfunc.cytfmra " + \
              "not available locally\n"
        msg += "(See documentation on compilation with cython)"
        sys.exit(msg)

    return globals()[name]()
