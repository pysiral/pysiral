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
ctypedef np.float32_t DTYPE_tf


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def cysmooth(double[:] x, int window):
    """ Numpy implementation of the IDL SMOOTH function """
    cdef np.ndarray[DTYPE_t, ndim=1] kernel
    kernel = np.ones(shape=(window))/float(window)
    y = np.convolve(x, kernel, mode='same')
    return y


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def cyfindpeaks(np.ndarray[DTYPE_t, ndim=1] data):
    """
    Finds peaks in `data` which are of `spacing` width and >=`limit`.
    :param data: values
    :param spacing: minimum spacing to the next peak (should be 1 or more)
    :param limit: peaks should have value greater or equal
    :return:
    """

    cdef int spacing, s, start
    # cdef double[:] h_b, h_c, h_a
    cdef int len = data.size
    cdef np.ndarray[DTYPE_t, ndim=1] x, h_b, h_c, h_a

    spacing = 1
    x = np.zeros(shape=(len+2*spacing))

    x[:spacing] = data[0]-1.e-6
    x[-spacing:] = data[-1]-1.e-6
    x[spacing:spacing+len] = data

    start = 0
    h_b = x[start:start + len]  # before
    start = spacing
    h_c = x[start:start + len]  # central
    start = spacing + 1
    h_a = x[start:start + len]  # after
    peak_candidate = np.logical_and(h_c > h_b, h_c > h_a)

    ind = np.argwhere(peak_candidate)
    ind = ind.reshape(ind.size)

    return ind


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def cyinterpolate(np.ndarray[DTYPE_tf, ndim=1] rng,
                  np.ndarray[DTYPE_tf, ndim=1] wfm,
                  int oversampling):
    cdef int n, n_os
    cdef float minval, maxval
    cdef np.ndarray[DTYPE_t, ndim=1] range_os
    n = len(rng)
    n_os = n*oversampling
    minval = np.amin(rng)
    maxval = np.amax(rng)
    range_os = np.linspace(minval, maxval, n_os)
    wfm_os = np.interp(range_os, rng, wfm)
    return range_os, wfm_os


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def cy_wfm_get_noise_level(double[:] wfm, int oversample_factor):
    """ According to CS2AWI TFMRA implementation """
    cpdef double[:] early_wfm = wfm[0:5*oversample_factor]
    cdef double noise_level
    noise_level = np.sum(early_wfm)/float(len(early_wfm))
    return noise_level


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def cy_normalize_wfm(np.ndarray[DTYPE_t, ndim=1] y):
    cdef int i
    cdef double norm = np.amax(y)
    cdef np.ndarray[DTYPE_t, ndim=1] normed_y = y/norm
    return normed_y, norm

#def cy_filter_waveform(np.ndarray[DTYPE_tf, ndim=1] rng,
#                       np.ndarray[DTYPE_tf, ndim=1] wfm,
#                       int oversampling,
#                       str interp_type,
#                       int window_size):
#    """
#    Return a filtered waveform: block filter smoothing of
#    oversampled original waveform
#    """
#
#    # Waveform Oversampling
#    cdef int n = len(rng)
#    cdef np.ndarray[DTYPE_t, ndim=1] range_os = np.linspace(
#        np.nanmin(rng), np.nanmax(rng), n*oversampling)
#    interpolator = interp1d(rng, wfm, kind=interp_type)
#    cdef np.ndarray[DTYPE_t, ndim=1] wfm_os = interpolator(range_os)
#
#    # Smoothing
#    wfm_os = cysmooth(wfm_os, window_size)
#
#    return range_os, wfm_os