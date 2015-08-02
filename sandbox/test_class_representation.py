# -*- coding: utf-8 -*-
"""
Created on Sat Aug 01 13:26:07 2015

@author: Stefan
"""

import numpy as np


class ClassA():

    def __new__(self):

        return np.array(shape=(10))

    def __init__(self):

        self.uncertainty = np.zeros(shape=(10))
        self.bias = np.zeros(shape=(10))


class L2ElevationArray(np.ndarray):
    """
    Recipe from:
    http://docs.scipy.org/doc/numpy/user/basics.subclassing.html
    XXX: not yet full slicing capability!
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

    def __getslice__(self, i, j):
        r = np.ndarray.__getslice__(self, i, j)
        r.uncertainty = r.uncertainty[i:j]
        r.bias = r.bias[i:j]
        return r

    def set_value(self, value):
        uncertainty = self.uncertainty
        bias = self.bias
        self = L2ElevationArray(value)
        setattr(self, "uncertainty", uncertainty)
        setattr(self, "bias", bias)

if __name__ == "__main__":
    a = L2ElevationArray(shape=(10))
    print a.uncertainty
    a.uncertainty = np.arange(10)

    a.set_value(np.arange(10))
    print np.shape(a)
    print a.uncertainty
#    print b.bias
    # print a[np.arange(1, 6, 2)].uncertainty

#    print a
#    print a.uncertainty
#    print a.bias
#    a.slice(np.arange(4))
#    print a
#    print a.uncertainty
#    print a.bias

