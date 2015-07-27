# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 18:09:30 2015

@author: Stefan
"""

import numpy as np


class BaseClassifier(object):

    def __init__(self):
        pass


class OCOGParameter(BaseClassifier):
    """
    Calculate OCOG Parameters (Amplitude, Width)
    Algorithm Source: retrack_ocog.pro from CS2AWI lib
    """

    def __init__(self, wfm):
        super(OCOGParameter, self).__init__()
        self._n = np.shape(wfm)[0]
        self._amplitude = np.ndarray(shape=(self._n))
        self._width = np.ndarray(shape=(self._n))
        self._calc_parameters(wfm)

    def _calc_parameters(self, wfm):
        for i in np.arange(self._n):
            y = wfm[i, :].flatten()
            noise = np.nanmean(y[0:11])
            y -= noise
            y[np.where(y < 0.0)[0]] = 0.0
            y2 = y**2.0
            self._amplitude[i] = np.sqrt((y2**2.0).sum() / y2.sum())
            self._width[i] = ((y2.sum())**2.0) / (y2**2.0).sum()

    @property
    def amplitude(self):
        return self._amplitude

    @property
    def width(self):
        return self._width
