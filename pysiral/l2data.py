# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 16:30:24 2015

@author: Stefan
"""
import numpy as np


class Level2Data(object):

    def __init__(self, l1b):
        self.info = np.copy(l1b.info)

    def set_surface_type(self, surface_type):
        self.surface_type = surface_type

    def update_retracked_range(self, retracked_range):
        pass

class L2ElevationArray(np.ndarray):
    """
    Recipe from:
    http://docs.scipy.org/doc/numpy/user/basics.subclassing.html
    XXX: not yet full slicing capability! -> __getitem__ trouble
         always use cls[list] and cls.uncertainty[list]
         cls[list].uncertainty will fail
    """

    def __new__(subtype, shape, dtype=float, buffer=None, offset=0,
                strides=None, order=None, info=None):
        obj = np.ndarray.__new__(
            subtype, shape, dtype, buffer, offset, strides, order)*np.nan
        obj.uncertainty = np.zeros(shape=shape, dtype=float)
        obj.bias = np.zeros(shape=shape, dtype=float)
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.uncertainty = getattr(obj, 'uncertainty', None)
        self.bias = getattr(obj, 'bias', None)

    def _setattr(self, obj, index):
        return obj

    def __getslice__(self, i, j):
        r = np.ndarray.__getslice__(self, i, j)
        r.uncertainty = r.uncertainty[i:j]
        r.bias = r.bias[i:j]
        return r
