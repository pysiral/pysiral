# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 14:06:28 2016

@author: shendric

Purpose: Test a faster implementation of the TFMRA retracker

"""

from pysiral.config import ConfigInfo
from pysiral.l1bdata import L1bConstructor
from pysiral.retracker import TFMRA

from scipy.interpolate import interp1d
import numpy as np

import os
import time

TEST_FILE = r"CS_OFFL_SIR_SAR_1B_20160401T013146_20160401T013933_C001.DBL"


def test_cython_tfmra():

    pysiral_config = ConfigInfo()

    # Get the full path to a test data file
    cs2_filename = os.path.join(pysiral_config.sandbox_path, "testdata",
                                "cryosat2", TEST_FILE)

    # Create l1b data object
    l1b = L1bConstructor(pysiral_config)
    l1b.mission = "cryosat2"
    l1b.filename = cs2_filename
    l1b.construct()

    print "Constructing of l1bdata object complete"

    # Get reference data set
    wfm = l1b.waveform.power
    rng = l1b.waveform.range
    rmode = l1b.waveform.radar_mode
    is_valid = l1b.waveform.is_valid
    indices = l1b.surface_type.get_by_name("ocean").indices

    print "Start old TFMRA retracker implementation"
    t0 = time.clock()
    tfmra = cTFMRA()
    tfmra.set_default_options()
    tfmra.init(l1b.n_records)
    tfmra.l2_retrack(rng, wfm, indices, rmode, is_valid)
    t1 = time.clock()
    print " . done in %.2f seconds" % (t1 - t0)

    stop


#    print "Test cython TFMRA retracker implementation"
#    tfmra = cTFMRA()
#    tfmra.set_default_options()
#    tfmra.init(l1b.n_records)
#    tfmra.l2_retrack(rng, wfm, indices, rmode, is_valid)
#    t2 = time.clock()
#    print " . done in %.2f seconds" % (t2 - t1)


from pysiral.retracker import BaseRetracker
from pysiral.performance.cytfmra import cysmooth, cyfindpeaks


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

            # Mandatory return function
            self._range[i] = tfmra_range + self._options.offset
            self._power[i] = tfmra_power * norm

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
        wfm_os = cysmooth(wfm_os, window_size)

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

        #try:
        peaks = cyfindpeaks(wfm[0:absolute_maximum_index])
        #except:
        #    return -1

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


def wfm_get_noise_level(wfm, oversample_factor):
    """ According to CS2AWI TFMRA implementation """
    return np.nanmean(wfm[0:5*oversample_factor])
#
#
#def smooth(x, window):
#    """ Numpy implementation of the IDL SMOOTH function """
#    return np.convolve(x, np.ones(window)/window, mode='same')
#
#
#def findpeaks(data, spacing=1, limit=None):
#    """Finds peaks in `data` which are of `spacing` width and >=`limit`.
#    :param data: values
#    :param spacing: minimum spacing to the next peak (should be 1 or more)
#    :param limit: peaks should have value greater or equal
#    :return:
#    """
#    len = data.size
#    x = np.zeros(len+2*spacing)
#    x[:spacing] = data[0]-1.e-6
#    x[-spacing:] = data[-1]-1.e-6
#    x[spacing:spacing+len] = data
#    peak_candidate = np.zeros(len)
#    peak_candidate[:] = True
#    for s in range(spacing):
#        start = spacing - s - 1
#        h_b = x[start:start + len]  # before
#        start = spacing
#        h_c = x[start:start + len]  # central
#        start = spacing + s + 1
#        h_a = x[start:start + len]  # after
#        peak_candidate = np.logical_and(
#            peak_candidate, np.logical_and(h_c > h_b, h_c > h_a))
#
#    ind = np.argwhere(peak_candidate)
#    ind = ind.reshape(ind.size)
#    if limit is not None:
#        ind = ind[data[ind] > limit]
#    return ind



## @autojit
#def tfmra_func(rng, wfm, radar_mode, indices, is_valid, opt):
#
#    tfmra_range = np.full(rng.shape[0], np.nan)
#    tfmra_power = np.full(rng.shape[0], np.nan)
#
#    for i in indices:
#
#        # Get the filtered waveform, index of first maximum & norm
#        filt_rng, filt_wfm, fmi, norm = tfmra_get_filtered_wfm(
#            rng[i, :], wfm[i, :], radar_mode[i], opt)
#
#        # first maximum finder might have failed
#        if fmi == -1:
#            tfmra_range[i] = np.nan
#            tfmra_power[i] = np.nan
#            return
#
#        # Get track point and its power
#        tfmra_threshold = opt["threshold"]
#        tfmra_rng, tfmra_pow = tfmra_get_threshold_range(
#            filt_rng, filt_wfm, fmi, tfmra_threshold)
#
#        # Mandatory return function
#        tfmra_range[i] = tfmra_rng
#        tfmra_power[i] = tfmra_pow * norm
#
#    return tfmra_range, tfmra_power
#
#
## @autojit
#def tfmra_get_filtered_wfm(rng, wfm, radar_mode, opt):
#
#    # Echo oversampling & filtering
#    # t0 = time.clock()
#    filt_rng, filt_wfm = tfmra_filter_waveform(rng, wfm, radar_mode, opt)
#
#    # t1 = time.clock()
#    # Normalize filtered waveform
#    norm = np.nanmax(filt_wfm)
#    filt_wfm = filt_wfm/norm
#
#    # t2 = time.clock()
#    # Get noise level in normalized units
#    oversampling = opt["wfm_oversampling_factor"]
#    noise_level = wfm_get_noise_level(filt_wfm, oversampling)
#
#    # Find first maxima
#    # (needs to be above radar mode dependent noise threshold)
#    fmnt = opt["first_maximum_normalized_threshold"][radar_mode]
#    peak_minimum_power = fmnt + noise_level
##    for i in range(10):
##        t3 = time.clock()
##        fmi = tfmra_get_first_maximum_index(filt_wfm, peak_minimum_power)
##        t4 = time.clock()
##        print "first maximum finder: %.5f seconds" % (t4 - t3)
#
##    for i in range(10):
##        t4 = time.clock()
#    fmi_c = tfmra_get_first_maximum_index_cython(filt_wfm, peak_minimum_power)
##        t5 = time.clock()
##        print "first maximum finder (cython): %.5f seconds" % (t5 - t4)
#
#    # print fmi, fmi_c
#
##    print "filter: %.5f seconds" % (t1 - t0)
##    print "norm: %.5f seconds" % (t2 - t1)
##    print "noise level: %.5f seconds" % (t3 - t2)
##    print "first maximum finder: %.5f seconds" % (t4 - t3)
##    print "first maximum finder (cython): %.5f seconds" % (t5 - t4)
#
#    # stop
#    return filt_rng, filt_wfm, fmi_c, norm
#
#
## @nb.jit(cache=True)
#def tfmra_filter_waveform(rng, wfm, radar_mode, opt):
#    """
#    Return a filtered waveform: block filter smoothing of
#    oversampled original waveform
#    """
#
#    # Parameter from options dictionary if omitted
#    oversampling = opt["wfm_oversampling_factor"]
#    interp_type = opt["wfm_oversampling_method"]
#    window_size = opt["wfm_smoothing_window_size"][radar_mode]
#
#    # Waveform Oversampling
#    n = len(rng)
#    range_os = np.linspace(np.nanmin(rng), np.nanmax(rng), n*oversampling)
#    interpolator = interp1d(rng, wfm, kind=interp_type)
#    wfm_os = interpolator(range_os)
#
#    # Smoothing
#    wfm_os = smooth(wfm_os, window_size)
#
#    return range_os, wfm_os
#
#
## @autojit
## numba compiler fails here
## @nb.jit(cache=True)
#def tfmra_get_first_maximum_index(wfm, peak_minimum_power):
#    """
#    Return the index of the first peak (first maximum) on
#    the leading edge before the absolute power maximum.
#    The first peak is only valid if its power exceeds a certain threshold
#    """
#
#    # Get the main maximum first
#    absolute_maximum_index = np.argmax(wfm)
#
#    # Find relative maxima before the absolute maximum
#    try:
#        peaks = findpeaks(wfm[0:absolute_maximum_index])
#    except:
#        return -1
#
#    # Check if relative maximum are above the required threshold
#    leading_maxima = np.where(wfm[peaks] >= peak_minimum_power)[0]
#
#    # Identify the first maximum
#    first_maximum_index = absolute_maximum_index
#    if len(leading_maxima) > 0:
#        # first_maximum_index = leading_maxima[0]
#        first_maximum_index = peaks[leading_maxima[0]]
#
#    return first_maximum_index
#
#from cyfindpeaks import cyfindpeaks
#
#def tfmra_get_first_maximum_index_cython(wfm, peak_minimum_power):
#    """
#    Return the index of the first peak (first maximum) on
#    the leading edge before the absolute power maximum.
#    The first peak is only valid if its power exceeds a certain threshold
#    """
#
#    # Get the main maximum first
##    t0 = time.clock()
#    absolute_maximum_index = np.argmax(wfm)
#
#    # Find relative maxima before the absolute maximum
#    try:
#        peaks = cyfindpeaks(wfm[0:absolute_maximum_index])
#    except:
#        return -1
#
#    # Check if relative maximum are above the required threshold
#    leading_maxima = np.where(wfm[peaks] >= peak_minimum_power)[0]
#
#    # Identify the first maximum
#    first_maximum_index = absolute_maximum_index
#    if len(leading_maxima) > 0:
#        # first_maximum_index = leading_maxima[0]
#        first_maximum_index = peaks[leading_maxima[0]]
#
#    return first_maximum_index
#
## @autojit
#def tfmra_get_threshold_range(rng, wfm, first_maximum_index, threshold):
#    """
#    Return the range value and the power of the retrack point at
#    a given threshold of the firsts maximum power
#    """
#
#    # get first index greater as threshold power
#    first_maximum_power = wfm[first_maximum_index]
#
#    # Get power of retracked point
#    tfmra_power = threshold*first_maximum_power
#
#    # Use linear interpolation to get exact range value
#    points = np.where(wfm[:first_maximum_index] > tfmra_power)[0]
#
#    if len(points) == 0:
#        return np.nan, np.nan
#
#    i0, i1 = points[0]-1, points[0]
#    gradient = (wfm[i1]-wfm[i0])/(rng[i1]-rng[i0])
#    tfmra_range = (tfmra_power - wfm[i0]) / gradient + rng[i0]
#
#    return tfmra_range, tfmra_power
#
#
## @autojit
#def wfm_get_noise_level(wfm, oversample_factor):
#    """ According to CS2AWI TFMRA implementation """
#    return np.nanmean(wfm[0:5*oversample_factor])
#
#
#@nb.jit(nb.i8[:](nb.f8[:]))
#def findpeaks(data):
#    """Finds peaks in `data` which are of `spacing` width and >=`limit`.
#    :param data: values
#    :param spacing: minimum spacing to the next peak (should be 1 or more)
#    :param limit: peaks should have value greater or equal
#    :return:
#    """
#    spacing = 1
#    # limit = None
#    len = data.size
#    x = np.zeros(len+2*spacing)
#    x[:spacing] = data[0]-1.e-6
#    x[-spacing:] = data[-1]-1.e-6
#    x[spacing:spacing+len] = data
#    peak_candidate = np.zeros(len)
#    peak_candidate[:] = True
#    for s in range(spacing):
#        start = spacing - s - 1
#        h_b = x[start:start + len]  # before
#        start = spacing
#        h_c = x[start:start + len]  # central
#        start = spacing + s + 1
#        h_a = x[start:start + len]  # after
#        condition = numba_and(h_c > h_b, h_c > h_a)
#        peak_candidate = np.logical_and(peak_candidate, condition)
#
#    ind = np.argwhere(peak_candidate)
#    ind = ind.reshape(ind.size)
##    if limit is not None:
##        ind = ind[data[ind] > limit]
#    return ind
#
#
#@nb.jit(nb.b1[:](nb.b1[:], nb.b1[:]))
#def numba_and(a, b):
#    return np.logical_and(a, b)
#
#
## @jit(n)
#def smooth(x, window):
#    """ Numpy implementation of the IDL SMOOTH function """
#    return np.convolve(x, np.ones(window)/window, mode='same')


if __name__ == "__main__":
    # simple_numba_test()
    test_cython_tfmra()
#    import profile
#    profile.run('test_cython_tfmra()')

