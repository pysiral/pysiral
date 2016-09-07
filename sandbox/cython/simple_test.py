# -*- coding: utf-8 -*-
"""
Created on Mon Sep 05 21:18:28 2016

@author: Stefan
"""

import numpy as np
import time

from pysiral.performance.cytfmra import (cyinterpolate, cysmooth,
                                         cy_wfm_get_noise_level,
                                         cy_normalize_wfm)

def test_smooth():

    n_tests = 10000
    a = np.arange(10000, dtype=np.float64)
    print "Test: smoothing"

    t0 = time.clock()
    for i in np.arange(n_tests):
        a = smooth(a, 10)
    t1 = time.clock()
    print "- python: %.4f seconds" % (t1-t0)

    t0 = time.clock()
    for i in np.arange(n_tests):
        a = cysmooth(a[:], 10)
    t1 = time.clock()
    print "- cython: %.4f seconds" % (t1-t0)

def smooth(x, window):
    """ Numpy implementation of the IDL SMOOTH function """
    return np.convolve(x, np.ones(window)/window, mode='same')


def test_interpolate():

    n_tests = 500
    a = np.arange(10000, dtype=np.float32)
    b = np.arange(10000, dtype=np.float32)
    print "\nTest: interpolate"

    t0 = time.clock()
    for i in np.arange(n_tests):
        c = interpolate(a, b, 10)
    t1 = time.clock()
    print "- python: %.4f seconds" % (t1-t0)

    t0 = time.clock()
    for i in np.arange(n_tests):
        c = cyinterpolate(a, b, 10)
    t1 = time.clock()
    print "- cython: %.4f seconds" % (t1-t0)


def test_noise_level():

    n_tests = 2000
    a = np.arange(10000, dtype=np.float64)

    print "\nTest: wfm_get_noise_level"
    t0 = time.clock()
    for i in np.arange(n_tests):
        c = wfm_get_noise_level(a, 10)
    t1 = time.clock()
    print "- python: %.4f seconds" % (t1-t0)

    t0 = time.clock()
    for i in np.arange(n_tests):
        c = cy_wfm_get_noise_level(a[:], 10)
    t1 = time.clock()

    print "- cython: %.4f seconds" % (t1-t0)

def test_normalize():

    n_tests = 2000
    a = np.arange(10000, dtype=np.float64)

    print "\nTest: normalize_wfm"
    t0 = time.clock()
    for i in np.arange(n_tests):
        c = normalize_wfm(a)
    t1 = time.clock()
    print "- python: %.4f seconds" % (t1-t0)

    t0 = time.clock()
    for i in np.arange(n_tests):
        c = cy_normalize_wfm(a[:])
    t1 = time.clock()

    print "- cython: %.4f seconds" % (t1-t0)

from scipy.interpolate import interp1d

def interpolate(rng, wfm, oversampling):
    """ Numpy implementation of the IDL SMOOTH function """
    n = len(rng)
    range_os = np.linspace(np.nanmin(rng), np.nanmax(rng), n*oversampling)
    interpolator = interp1d(rng, wfm, kind="linear")
    wfm_os = interpolator(range_os)
    return wfm_os

def wfm_get_noise_level(wfm, oversample_factor):
    """ According to CS2AWI TFMRA implementation """
    return np.nanmean(wfm[0:5*oversample_factor])


def normalize_wfm(y):
    norm = np.nanmax(y)
    return y/norm, norm

if __name__ == "__main__":
    test_smooth()
    test_interpolate()
    test_noise_level()
    test_normalize()
