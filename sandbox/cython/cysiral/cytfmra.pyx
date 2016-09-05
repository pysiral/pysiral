# -*- coding: utf-8 -*-
"""
Created on Mon Sep 05 19:09:08 2016

@author: Stefan
"""


from scipy.interpolate import interp1d

import cython
cimport cython
import numpy as np
cimport numpy as np

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

#class cTFMRACore():
#    """ Default Retracker from AWI CryoSat-2 production system """
#
#    DOCSTR = r"Threshold first maximum retracker (TFMRA)"
#
#    def __init__(self, options):
#        self._options = options
#
#
#    def retrack(self, rng, wfm, indices, radar_mode, is_valid):
#        """ API Calling method """
#
#        self.range = np.ndarray(shape=radar_mode.shape)
#        self.power = np.ndarray(shape=radar_mode.shape)
#
#        for i in indices:
#
#            # Get the filtered waveform, index of first maximum & norm
#            filt_rng, filt_wfm, fmi, norm = self.get_filtered_wfm(
#                rng[i, :], wfm[i, :], radar_mode[i])
#
#            # first maximum finder might have failed
#            if fmi == -1:
#                self.range[i] = np.nan
#                self.power[i] = np.nan
#                return
#
#            # Get track point and its power
#            tfmra_threshold = self._options.threshold
#            tfmra_range, tfmra_power = self.get_threshold_range(
#                filt_rng, filt_wfm, fmi, tfmra_threshold)
#
#            # Mandatory return function
#            self.range[i] = tfmra_range + self._options.offset
#            self.power[i] = tfmra_power * norm
#
#    def get_preprocessed_wfm(self, rng, wfm, radar_mode, is_valid):
#        """
#        Returns the intermediate product (oversampled range bins,
#        oversampled and filtered waveforms, indices of first maxima
#        and peak power norm for custom applications
#        """
#
#        oversample_factor = self._options.wfm_oversampling_factor
#        wfm_shape = (wfm.shape[0], wfm.shape[1]*oversample_factor)
#
#        filt_rng = np.full(wfm_shape, np.nan)
#        filt_wfm = np.full(wfm_shape, np.nan)
#        fmi = np.full(wfm_shape[0], -1, dtype=np.int32)
#        norm = np.full(wfm_shape[0], np.nan)
#
#        for i in np.arange(wfm.shape[0]):
#            if not is_valid[i]:
#                continue
#            result = self.get_filtered_wfm(rng[i, :], wfm[i, :], radar_mode[i])
#            filt_rng[i, :] = result[0]
#            filt_wfm[i, :] = result[1]
#            fmi[i] = result[2]
#            norm[i] = result[3]
#
#        return filt_rng, filt_wfm, fmi, norm
#
#    def get_thresholds_distance(self, rng, wfm, fmi, t0, t1):
#        """
#        Return the distance between two thresholds t0 < t1
#        """
#        width = np.full(rng.shape[0], np.nan, dtype=np.float32)
#        for i in np.arange(rng.shape[0]):
#            if fmi[i] is None:
#                continue
#            r0 = self.get_threshold_range(rng[i, :], wfm[i, :], fmi[i], t0)
#            r1 = self.get_threshold_range(rng[i, :], wfm[i, :], fmi[i], t1)
#            width[i] = r1[0] - r0[0]
#
#        # some irregular waveforms might produce negative width values
#        is_negative = np.where(width < 0.)[0]
#        width[is_negative] = np.nan
#        return width
#
#    def get_filtered_wfm(self, rng, wfm, radar_mode):
#
#        # Echo oversampling & filtering
#        filt_rng, filt_wfm = self.filter_waveform(rng, wfm, radar_mode)
#
#        # Normalize filtered waveform
#        filt_wfm, norm = self.normalize_wfm(filt_wfm)
#
#        # Get noise level in normalized units
#        oversampling = self._options.wfm_oversampling_factor
#        noise_level = wfm_get_noise_level(filt_wfm, oversampling)
#
#        # Find first maxima
#        # (needs to be above radar mode dependent noise threshold)
#        fmnt = self._options.first_maximum_normalized_threshold[radar_mode]
#        peak_minimum_power = fmnt + noise_level
#        fmi = self.get_first_maximum_index(filt_wfm, peak_minimum_power)
#
#        return filt_rng, filt_wfm, fmi, norm
#
#    def filter_waveform(self, rng, wfm, radar_mode):
#        """
#        Return a filtered waveform: block filter smoothing of
#        oversampled original waveform
#        """
#
#        # Parameter from options dictionary if omitted
#        opt = self._options
#        oversampling = opt.wfm_oversampling_factor
#        interp_type = opt.wfm_oversampling_method
#        window_size = opt.wfm_smoothing_window_size[radar_mode]
#
#        # Waveform Oversampling
#        n = len(rng)
#        range_os = np.linspace(np.nanmin(rng), np.nanmax(rng), n*oversampling)
##        interpolator = interp1d(rng, wfm, kind=interp_type)
##        wfm_os = interpolator(range_os)
#        wfm_os = np.interp(range_os, rng, wfm)
#
#        # Smoothing
#        wfm_os = cysmooth(wfm_os, window_size)
#
#        return range_os, wfm_os
#
#    def normalize_wfm(self, y):
#        norm = np.nanmax(y)
#        return y/norm, norm
#
#    def get_first_maximum_index(self, wfm, peak_minimum_power):
#        """
#        Return the index of the first peak (first maximum) on
#        the leading edge before the absolute power maximum.
#        The first peak is only valid if its power exceeds a certain threshold
#        """
#
#        # Get the main maximum first
#        absolute_maximum_index = np.argmax(wfm)
#
#        # Find relative maxima before the absolute maximum
#        try:
#            peaks = cyfindpeaks(wfm[0:absolute_maximum_index])
#        except:
#            return -1
#
#        # Check if relative maximum are above the required threshold
#        leading_maxima = np.where(wfm[peaks] >= peak_minimum_power)[0]
#
#        # Identify the first maximum
#        first_maximum_index = absolute_maximum_index
#        if len(leading_maxima) > 0:
#            # first_maximum_index = leading_maxima[0]
#            first_maximum_index = peaks[leading_maxima[0]]
#
#        return first_maximum_index
#
#    def get_threshold_range(self, rng, wfm, first_maximum_index, threshold):
#        """
#        Return the range value and the power of the retrack point at
#        a given threshold of the firsts maximum power
#        """
#
#        # get first index greater as threshold power
#        first_maximum_power = wfm[first_maximum_index]
#
#        # Get power of retracked point
#        tfmra_power = threshold*first_maximum_power
#
#        # Use linear interpolation to get exact range value
#        points = np.where(wfm[:first_maximum_index] > tfmra_power)[0]
#
#        if len(points) == 0:
#            return np.nan, np.nan
#
#        i0, i1 = points[0]-1, points[0]
#        gradient = (wfm[i1]-wfm[i0])/(rng[i1]-rng[i0])
#        tfmra_range = (tfmra_power - wfm[i0]) / gradient + rng[i0]
#
#        return tfmra_range, tfmra_power


#def wfm_get_noise_level(wfm, oversample_factor):
#    """ According to CS2AWI TFMRA implementation """
#    return np.nanmean(wfm[0:5*oversample_factor])


def cysmooth(np.ndarray[DTYPE_t, ndim=1] x, int window):
    """ Numpy implementation of the IDL SMOOTH function """
    cdef np.ndarray[DTYPE_t, ndim=1] kernel = np.ones(window)/float(window)
    cdef np.ndarray[DTYPE_t, ndim=1] y = np.convolve(x, kernel, mode='same')
    return y


def cyfindpeaks(np.ndarray[DTYPE_t, ndim=1] data):
    """
    Finds peaks in `data` which are of `spacing` width and >=`limit`.
    :param data: values
    :param spacing: minimum spacing to the next peak (should be 1 or more)
    :param limit: peaks should have value greater or equal
    :return:
    """

    cdef int spacing = 1
    cdef int len = data.size
    cdef np.ndarray[DTYPE_t, ndim=1] x = np.zeros(len+2*spacing)

    x[:spacing] = data[0]-1.e-6
    x[-spacing:] = data[-1]-1.e-6
    x[spacing:spacing+len] = data

    cdef np.ndarray[DTYPE_t, ndim=1] peak_candidate = np.zeros(len)
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

    return ind
