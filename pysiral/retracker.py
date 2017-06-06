# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 15:48:58 2015

@author: Stefan
"""

# cythonized bottleneck functions for cTFMRA
try:
    from pysiral.bnfunc.cytfmra import (cytfmra_findpeaks, cytfmra_interpolate,
                                        cytfmra_wfm_noise_level,
                                        cytfmra_normalize_wfm)
    CYTFMRA_OK = True
except:
    CYTFMRA_OK = False


from pysiral.flag import ANDCondition, FlagContainer

from treedict import TreeDict
import numpy as np
import sys

# Utility methods for retracker:
from scipy.interpolate import interp1d
import bottleneck as bn


class BaseRetracker(object):
    """
    Main Retracker Class (all retrackers must be of instance BaseRetracker)
    """

    def __init__(self):
        self._indices = None
        self._classifier = None
        self._l1b = None
        self._l2 = None

    def set_options(self, **opt_dict):
        # TODO: Create options object
        self._options = TreeDict.fromdict(opt_dict, expand_nested=True)

    def set_indices(self, indices):
        # TODO: Validation
        self._indices = indices

    def set_classifier(self, classifier):
        self._classifier = classifier

    def init(self, n_records):
        self._create_default_properties(n_records)

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

    def get_l1b_parameter(self, data_group, parameter_name):
        """ Get any valid level-2 paremeter name """
        try:
            return self._l1b.get_parameter_by_name(data_group, parameter_name)
        except:
            return None

    def get_l2_parameter(self, parameter_name):
        """ Get any valid level-2 paremeter name """
        try:
            return self._l2.get_parameter_by_name(parameter_name)
        except:
            return None

    def _create_default_properties(self, n_records):
        # XXX: Currently only range and status (False: ok)
        for parameter in ["_range", "_power"]:
            setattr(self, parameter, np.full((n_records), np.nan))
        self._uncertainty = np.full((n_records), 0.0, dtype=np.float32)
        self._flag = np.zeros(shape=(n_records), dtype=np.bool)

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
        sitype = self._l2.sitype
        tfmra_threshold = self.get_tfmra_threshold(sigma0, sitype, indices)

        for i in indices:

            # Remove artifact bins
            wfm_skipped = wfm[i, :]
            wfm_skipped[0:5] = 0

            # Get the filtered waveform, index of first maximum & norm
            filt_rng, filt_wfm, fmi, norm = self.get_filtered_wfm(
                rng[i, :], wfm_skipped, radar_mode[i])

            # first maximum finder might have failed
            if fmi == -1:
                self._range[i] = np.nan
                self._power[i] = np.nan
                return

            # Get track point and its power
            tfmra_range, tfmra_power = self.get_threshold_range(
                filt_rng, filt_wfm, fmi, tfmra_threshold[i])

            # Mandatory return function
            self._range[i] = tfmra_range + self._options.offset
            self._power[i] = tfmra_power * norm

        if "uncertainty" in self._options:
            if self._options.uncertainty.type == "fixed":
                self._uncertainty[:] = self._options.uncertainty.value

    def get_tfmra_threshold(self, sigma0, sitype, indices):

        # short link to options
        option = self._options.threshold

        # scale sigma0 (may be necessary for inter-mission sigma0 biases)
        if "sigma_scaling" in option:
            sigma0 = sigma0 * option.sigma_scaling

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
            threshold[indices] = value

        # dependend in sitype and sigma0
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

    def normalize_wfm(self, y):
        norm = np.nanmax(y)
        return y/norm, norm

    def get_first_maximum_index(self, wfm, peak_minimum_power):
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

    def get_threshold_range(self, rng, wfm, first_maximum_index, threshold):
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
    """ Default Retracker from AWI CryoSat-2 production system """

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

    def normalize_wfm(self, y):
        norm = np.nanmax(y)
        return y/norm, norm

    def get_first_maximum_index(self, wfm, peak_minimum_power):
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

    def get_threshold_range(self, rng, wfm, first_maximum_index, threshold):
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
    """ Default Retracker from AWI CryoSat-2 production system """

    DOCSTR = r"Threshold first maximum retracker (TFMRA)"

    def __init__(self):
        super(cTFMRA, self).__init__()

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
                continue

            # Get track point and its power
            tfmra_threshold = self._options.threshold
            tfmra_range, tfmra_power = self.get_threshold_range(
                filt_rng, filt_wfm, fmi, tfmra_threshold)

            # Set the values
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

        filt_rng, filt_wfm = self.filter_waveform(rng, wfm, radar_mode)

        # Normalize filtered waveform
        filt_wfm, norm = cytfmra_normalize_wfm(filt_wfm)

        # Get noise level in normalized units
        oversampling = self._options.wfm_oversampling_factor
        noise_level = cytfmra_wfm_noise_level(filt_wfm, oversampling)

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
        # interp_type = opt.wfm_oversampling_method
        window_size = opt.wfm_smoothing_window_size[radar_mode]

        # Use cython implementation of waveform oversampling
        range_os, wfm_os = cytfmra_interpolate(rng, wfm, oversampling)

        # Smoothing
        wfm_os = bnsmooth(wfm_os, window_size)

        return range_os, wfm_os

    def get_first_maximum_index(self, wfm, peak_minimum_power):
        """
        Return the index of the first peak (first maximum) on
        the leading edge before the absolute power maximum.
        The first peak is only valid if its power exceeds a certain threshold
        """

        # Get the main maximum first
        absolute_maximum_index = bn.nanargmax(wfm)

        # Find relative maxima before the absolute maximum

        try:
            peaks = cytfmra_findpeaks(wfm[0:absolute_maximum_index])
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

    def get_threshold_range(self, rng, wfm, first_maximum_index, threshold):
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


class NoneRetracker(BaseRetracker):
    """
    A dummy retracker that just returns NaN's but does not flag
    the result as invalid. Should be used if a certain surface type
    should not be used
    """

    def __init__(self):
        super(NoneRetracker, self).__init__()

    def create_retracker_properties(self, n_records):
        pass

    def l2_retrack(self, range, wfm, indices, radar_mode, is_valid):
        self._range[indices] = np.nan
        self._power[indices] = np.nan
        self._flag[indices] = False


class SICCILead(BaseRetracker):

    def __init__(self):
        super(SICCILead, self).__init__()

    def create_retracker_properties(self, n_records):
        parameter = ["retracked_bin", "maximum_power_bin", "sigma", "k",
                     "alpha", "power_in_echo_tail", "rms_echo_and_model"]
        for parameter_name in parameter:
            setattr(self, parameter_name,
                    np.ndarray(shape=(n_records), dtype=np.float32) * np.nan)

    def l2_retrack(self, range, wfm, indices, radar_mode, is_valid):
        # Run the retracker
        self._sicci_lead_retracker(range, wfm, indices)
        # Filter the results
        if self._options.filter.use_filter:
            self._filter_results()

    def _sicci_lead_retracker(self, range, wfm, indices):
        from scipy.optimize import curve_fit
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

# %% Function for CryoSat-2 based retracker

def wfm_get_noise_level(wfm, oversample_factor):
    """ According to CS2AWI TFMRA implementation """
    return np.nanmean(wfm[0:5*oversample_factor])


def smooth(x, window):
    """ Numpy implementation of the IDL SMOOTH function """
    return np.convolve(x, np.ones(window)/window, mode='same')

def bnsmooth(x, window):
    """ Bottleneck implementation of the IDL SMOOTH function """
    pad = (window-1)/2
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
    F_L = np.logical_and(
        t > t_0, t < (t_b+t_0)) * (aaa*t_diff**3 + aa*t_diff**2+t_diff/sigma)
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
