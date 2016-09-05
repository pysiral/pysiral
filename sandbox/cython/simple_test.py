# -*- coding: utf-8 -*-
"""
Created on Mon Sep 05 21:18:28 2016

@author: Stefan
"""

import numpy as np
import time


def test_smooth():

    t0 = time.clock()
    for i in np.arange(10000):
        a = smooth(np.arange(10000, dtype=np.float64), 10)
    t1 = time.clock()
    print t1-t0

    from cytfmra import cysmooth
    t0 = time.clock()
    for i in np.arange(10000):
        a = cysmooth(np.arange(10000, dtype=np.float64), 10)
    t1 = time.clock()
    print t1-t0

def smooth(x, window):
    """ Numpy implementation of the IDL SMOOTH function """
    return np.convolve(x, np.ones(window)/window, mode='same')

if __name__ == "__main__":
    test_smooth()
