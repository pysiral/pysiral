# -*- coding: utf-8 -*-

"""
This OCOG implementation comes from the first phase of the ESA Climate Change Initiative on Sea Ice
"""

from typing import Tuple

import numpy as np
from scipy.interpolate import interp1d

from pysiral.core.flags import ANDCondition
from pysiral.retracker import BaseRetracker

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"


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
        """ These thresholds are based on the SICCI code"""
        thrs = self._options.filter

        valid = ANDCondition()
        valid.add(self.leading_edge_width < thrs.maximum_leading_edge_width)
        valid.add(self.tail_shape < thrs.maximum_echo_tail_line_deviation)
        valid.add(self.retracked_bin > thrs.sensible_seaice_retracked_bin[0])
        valid.add(self.retracked_bin < thrs.sensible_seaice_retracked_bin[1])

        # Error flag is also computed for other surface types, do not override those
        error_flag = self._flag
        error_flag[self.indices] = np.logical_not(valid.flag[self.indices])
        self._flag = error_flag


def ocog_tail_shape(wfm, tracking_point, tail_pad=3):
    """ From SICCI module """
    tail = wfm[tracking_point+tail_pad:]
    tail = tail / np.mean(tail)
    P = np.polyfit(np.arange(len(tail)), tail, 1)
    residual = (tail-np.polyval(P, np.arange(len(tail))))
    return np.sqrt(np.sum(residual*residual) / len(residual))


def ocog_func(wave, percentage, skip):
    waveform = wave*wave
    sq_sum = np.sum(waveform)
    waveform = waveform*waveform
    qa_sum = np.sum(waveform)
    # Calculate retracking threshold (2.6.2 in ATDBv0)
    threshold = percentage * np.sqrt(qa_sum / sq_sum)
    ind_first_over = np.min(np.nonzero(wave > threshold))
    decimal = (wave[ind_first_over-1] - threshold) / (wave[ind_first_over-1] - wave[ind_first_over])
    return skip + ind_first_over - 1 + decimal


def ocog_properties(waveform: np.ndarray) -> Tuple[float, float, float]:
    """
    Computes the three OCOG properties center of gravity, amplitude & width
    of a waveform in range gate units.

    :param waveform: The waveform array

    :return: center of gravity (cog), amplitude & width
    """

    idx = np.arange(waveform.size) + 1
    waveform_sum = np.sum(waveform)
    waveform_scaled_sum = np.sum(idx*waveform)
    # waveform_squared_sum = np.sum(np.power(waveform, 2))
    waveform_squared_sum = np.sum(np.power(waveform, 2))

    # self._amplitude[i] = np.sqrt((y2 ** 2.0).sum() / y2.sum())
    #     # self._width[i] = ((y2.sum()) ** 2.0) / (y2 ** 2.0).sum()

    cog = waveform_scaled_sum / waveform_sum
    amplitude = np.sqrt(waveform_squared_sum / waveform_sum)
    width = waveform_sum**2. / waveform_squared_sum

    return cog, amplitude, width
