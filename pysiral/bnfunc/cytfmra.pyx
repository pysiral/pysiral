# -*- coding: utf-8 -*-
"""
Created on Mon Sep 05 19:09:08 2016

@author: Stefan
"""

from scipy.interpolate import interp1d

import cython
import numpy as np
import bottleneck as bn

cimport cython
cimport numpy as np
# cimport bottleneck as bn

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t
ctypedef np.float32_t DTYPE_tf


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def cytfmra_findpeaks(np.ndarray[DTYPE_t, ndim=1] data):
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
    cdef int xshape = len+2*spacing
    x = np.ndarray(shape=(xshape))

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

    ind = np.where(peak_candidate)[0]

    return ind


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def cytfmra_interpolate(np.ndarray[DTYPE_tf, ndim=1] rng,
                        np.ndarray[DTYPE_tf, ndim=1] wfm,
                       int oversampling):
    cdef int n, n_os
    cdef float minval, maxval
    cdef np.ndarray[DTYPE_t, ndim=1] range_os
    n = len(rng)
    n_os = n*oversampling

    minval = rng[0]
    maxval = rng[n-1]

    # np.arange is faster than linspace
    cdef float step = float((maxval-minval))/float(n_os-1)
    range_os = np.arange(n_os)*step+minval

    wfm_os = np.interp(range_os, rng, wfm)
    return range_os, wfm_os


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def cytfmra_wfm_noise_level(double[:] wfm, int oversample_factor):
    """ According to CS2AWI TFMRA implementation """
    cpdef double[:] early_wfm = wfm[0:5*oversample_factor]
    cdef double noise_level
    noise_level = np.sum(early_wfm)/float(len(early_wfm))
    return noise_level


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def cytfmra_normalize_wfm(np.ndarray[DTYPE_t, ndim=1] y):
    cdef double norm = bn.nanmax(y)
    cdef np.ndarray[DTYPE_t, ndim=1] normed_y = y/norm
    return normed_y, norm
