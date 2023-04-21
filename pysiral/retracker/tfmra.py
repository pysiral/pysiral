# -*- coding: utf-8 -*-

"""

"""

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"

from typing import Tuple

import bottleneck as bn
import numpy as np
from loguru import logger

from pysiral.l2proc.procsteps import Level2ProcessorStep
from pysiral.retracker import BaseRetracker

# cythonized bottleneck functions for cTFMRA
try:
    from pysiral.retracker.cytfmra import (cytfmra_findpeaks,
                                           cytfmra_interpolate,
                                           cytfmra_normalize_wfm,
                                           cytfmra_wfm_noise_level)
    CYTFMRA_OK = True
except ImportError:
    logger.error("Cannot import cytfmra")
    CYTFMRA_OK = False


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
        return {
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
            "first_maximum_ignore_leading_bins": 0,
        }

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
            first_maximum_equal_total_maximum = self._options.get("first_maximum_equal_total_maximum", False)
            if first_maximum_equal_total_maximum:
                fmi = np.argmax(filt_wfm)
            else:
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

        if (
            "uncertainty" in self._options
            and self._options.uncertainty.type == "fixed"
        ):
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
            a peak must surpass to be regarded as a first maximum candidate
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
        thresholds = np.arange(threshold_min, threshold_max + 0.1 * threshold_res, threshold_res)

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
        dim_dict = {"new_dims": (("tfmra_thresholds", thresholds.size),),
                    "dimensions": ("time", "tfmra_thresholds"),
                    "add_dims": (("tfmra_thresholds", thresholds),)}
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


def bnsmooth(x, window):
    """ Bottleneck implementation of the IDL SMOOTH function """
    pad = int((window-1)/2)
    n = len(x)
    xpad = np.ndarray(shape=(n+window))
    xpad[:pad] = 0.0
    xpad[pad:n+pad] = x
    xpad[n+pad:] = 0.0
    return bn.move_mean(xpad, window=window, axis=0)[window-1:(window+n-1)]