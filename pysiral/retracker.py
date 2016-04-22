# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 15:48:58 2015

@author: Stefan
"""

from treedict import TreeDict
import numpy as np


# Utility methods for retracker:
from scipy.interpolate import interp1d


class BaseRetracker(object):
    """
    Main Retracker Class (all retrackers must be of instance BaseRetracker)
    """

    def __init__(self):
        self._surface_type_id = None

    def set_options(self, **opt_dict):
        # TODO: Create options object
        self._options = TreeDict.fromdict(opt_dict, expand_nested=True)

    def set_surface_type(self, surface_type):
        # TODO: Validation
        self._surface_type_id = surface_type

    def retrack(self, l1b, l2):
        # Initialize the retracked range with an NaN array
        # -> only waveforms of type "surface_type" will retracked
        self._create_default_properties(l1b.n_records)
        # Give the Retracker the possibility to create additional
        # data arrays (All retracker must have this method)
        self._create_retracker_properties(l1b.n_records)
        # Set the indices for the given surface type
        self._set_surface_type_indices(l2)
        # Check if there is something to do
        if len(self._index) == 0:
            return False
        # Loop over each waveform of given surface type
        for index in self._index:
            print index
            self._retrack(l1b.waveform.range[index, :],
                          l1b.waveform.power[index, :],
                          index)
        return True

    def _set_surface_type_indices(self, l2):
        surface_type = l2.surface_type.get_by_name(self._surface_type_id)
        self._index = surface_type.indices

    def _create_default_properties(self, n_records):
        # XXX: Currently only range and status (true: ok)
        self._range = np.ndarray(shape=(n_records))*np.nan
        self._power = np.ndarray(shape=(n_records))*np.nan
        self._status = np.zeros(shape=(n_records), dtype=np.bool)

    def _add_retracked_point(self, range, power, index, **keyw):
        # XXX: Very basic
        if np.isfinite(range):
            self._status[index] = True
            self._range[index] = range
            self._power[index] = power

    @property
    def range(self):
        return self._range

    @property
    def power(self):
        return self._power

    @property
    def index(self):
        return self._index


class TFMRA(BaseRetracker):

    def __init__(self):
        super(TFMRA, self).__init__()

    def _create_retracker_properties(self, n_records):
        # None so far
        pass

    def _retrack(self, range, wfm, index):
        # Echo oversampling & filtering
        x, y = self._filter_waveform(range, wfm)
        # Normalize filtered waveform
        y, norm = self._normalize_wfm(y)
        # Get noise level in normalized units
        oversampling = self._options.wfm_oversampling_factor
        noise_level = wfm_get_noise_level(wfm, oversampling)
        # Find first maxima
        first_maximum_index = self._get_first_maxima_index(y, noise_level)
        # Get track point
        tfmra_range, tfmra_power = self._get_retracker_range(
            x, y, first_maximum_index, norm)
        # Mandatory return function
        self._add_retracked_point(tfmra_range, tfmra_power, index)

    def _filter_waveform(self, range, wfm):
        n = len(range)
        # Waveform Oversampling
        oversampling = self._options.wfm_oversampling_factor
        range_os = np.linspace(
            np.nanmin(range), np.nanmax(range), n*oversampling)
        interpolator_type = self._options.wfm_oversampling_method
        interpolator = interp1d(range, wfm, kind=interpolator_type)
        wfm_os = interpolator(range_os)
        # Smoothing
        wfm_os = smooth(wfm_os, self._options.wfm_smoothing_window_size)
        return range_os, wfm_os

    def _normalize_wfm(self, y):
        norm = np.nanmax(y)
        return y/norm, norm

    def _get_first_maxima_index(self, y, noise_level):
        from scipy import signal
        # Get the main maximum first
        absolute_maximum_index = np.nanargmax(y)
        # Find relative maxima before the absolute maximum
        order = self._options.first_maximum_local_order
        peaks = signal.argrelmax(y[0:absolute_maximum_index], order=order)[0]
        # Check if relative maximum are above the required threshold
        threshold = self._options.first_maximum_normalized_threshold
        leading_maxima = np.where(y[peaks] >= threshold + noise_level)[0]
        # Identify the first maximum
        first_maximum_index = absolute_maximum_index
        if len(leading_maxima) > 0:
            first_maximum_index = peaks[leading_maxima[0]]
        return first_maximum_index

    def _get_retracker_range(self, x, y, first_maximum_index, norm):
        # get first index greater as threshold power
        first_maximum_power = y[first_maximum_index]
        tfmra_threshold = self._options.threshold
        tfmra_power_normed = tfmra_threshold*first_maximum_power
        points = np.where(y > tfmra_power_normed)[0]
        # Use linear interpolation to get exact range value
        i0, i1 = points[0]-1, points[0]
        gradient = (y[i1]-y[i0])/(x[i1]-x[i0])
        tfmra_range = (tfmra_power_normed - y[i0]) / gradient + x[i0]
        # Add optional offset
        tfmra_range += self._options.offset
        # Get power of retracked point
        tfmra_power = tfmra_power_normed*norm
        return tfmra_range, tfmra_power


def wfm_get_noise_level(wfm, oversample_factor):
    """ According to CS2AWI TFMRA implementation """
    return np.nanmean(wfm[0:5*oversample_factor])


def smooth(x, window):
    """ Numpy implementation of the IDL SMOOTH function """
    return np.convolve(x, np.ones(window)/window, mode='same')
