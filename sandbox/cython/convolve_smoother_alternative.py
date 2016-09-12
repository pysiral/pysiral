# -*- coding: utf-8 -*-
"""
Created on Wed Sep 07 16:51:15 2016

@author: Stefan
"""

import time
import numpy as np
import matplotlib.pyplot as plt

def test_concolve_smoother_alternative():

    window = 11
    signal = np.ones(100)
    signal[::20] = 0

    t0 = time.clock()
    smoothed_signal_1 = smooth(signal, window)
    t1 = time.clock()
    print t1-t0

    t0 = time.clock()
    smoothed_signal_2 = smooth_alternative(signal, window)
    t1 = time.clock()
    print t1-t0

#    t0 = time.clock()
#    smoothed_signal_3 = smooth_alternative2(signal, window)
#    t1 = time.clock()
#    print t1-t0

    plt.figure("signals")
    plt.plot(signal, 'o')
    plt.plot(smoothed_signal_1, color="black", lw=2)
    plt.plot(smoothed_signal_2, color="red")
#    plt.plot(smoothed_signal_3, color="blue")

    plt.figure("smoother difference")
    plt.plot(smoothed_signal_1-smoothed_signal_2)
    plt.show()



def smooth(x, window):
    """ Numpy implementation of the IDL SMOOTH function """
    return np.convolve(x, np.ones(window)/window, mode='same')

import bottleneck as bn

def smooth_alternative(x, window):
    """ Numpy implementation of the IDL SMOOTH function """
    pad = (window-1)/2
    n = len(x)
    xpad = np.ndarray(shape=(n+window))
    xpad[0:pad] = 0.0
    xpad[pad:n+pad] = x
    xpad[n+pad:] = 0.0
    return bn.move_mean(xpad, window=window, axis=0)[window-1:(window+n-1)]

#def smooth_alternative2(x, window):
#    """ Numpy implementation of the IDL SMOOTH function """
#    pad = (window-1)/2
#    n = len(x)
#    xpad = np.ndarray(shape=(n+window))
#    xpad[0:pad] = 0.0
#    xpad[pad:n+pad] = x
#    xpad[n+pad:] = 0.0
#    kernel = np.ones(window)/window
#    return np.tensordot(xpad, kernel)[window-1:(window+n-1)]

#def smooth_alternative(x, window):
#    kernel = np.ones(window)/window
#    pad = (window-1)/2
#    n = len(x)
#    xpad = np.ndarray(shape=(n+window))
#    y = np.ndarray(shape=(n))
#    xpad[0:pad] = 0.0
#    xpad[pad:n+pad] = x
#    xpad[n+pad:] = 0.0
#    indices = np.arange(n)
#    for i in indices:
#        y[i] = sum(xpad[i:(i+window)]*kernel)
#    return y

if __name__ == "__main__":
    test_concolve_smoother_alternative()