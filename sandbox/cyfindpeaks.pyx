# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 17:10:27 2016

@author: shendric
"""
import numpy as np

def cyfindpeaks(data):
    """Finds peaks in `data` which are of `spacing` width and >=`limit`.
    :param data: values
    :param spacing: minimum spacing to the next peak (should be 1 or more)
    :param limit: peaks should have value greater or equal
    :return:
    """
    spacing = 1
    # limit = None
    len = data.size
    x = np.zeros(len+2*spacing)
    x[:spacing] = data[0]-1.e-6
    x[-spacing:] = data[-1]-1.e-6
    x[spacing:spacing+len] = data
    peak_candidate = np.zeros(len)
    peak_candidate[:] = True
    for s in range(spacing):
        start = spacing - s - 1
        h_b = x[start:start + len]  # before
        start = spacing
        h_c = x[start:start + len]  # central
        start = spacing + s + 1
        h_a = x[start:start + len]  # after
        condition = np.logical_and(h_c > h_b, h_c > h_a)
        peak_candidate = np.logical_and(peak_candidate, condition)

    ind = np.argwhere(peak_candidate)
    ind = ind.reshape(ind.size)
#    if limit is not None:
#        ind = ind[data[ind] > limit]
    return ind