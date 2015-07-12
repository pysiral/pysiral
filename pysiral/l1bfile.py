# -*- coding: utf-8 -*-
"""
Created on Tue Jul 07 14:10:34 2015

@author: Stefan
"""


class L1bData(object):
    """
    Unified L1b Data Class
    """
    def __init__(self):

        self._waveform_group = None
        self._orbit_group = None
        self._geocorrection_group = None
        self._classifier_group = None

