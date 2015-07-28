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
