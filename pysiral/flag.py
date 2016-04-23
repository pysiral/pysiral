# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 15:44:22 2016

@author: Stefan
"""

import numpy as np


class FlagContainer(object):

    def __init__(self, flag):
        self._flag = flag

    @property
    def indices(self):
        return np.where(self._flag)[0]

    @property
    def flag(self):
        return self._flag

    @property
    def num(self):
        return len(self.indices)


class ANDCondition(object):

    def __init__(self):
        self.flag = None

    def add(self, flag):
        if self.flag is None:
            self.flag = flag
        else:
            self.flag = np.logical_and(self.flag, flag)


class ORCondition(object):

    def __init__(self):
        self.flag = None

    def add(self, flag):
        if self.flag is None:
            self.flag = flag
        else:
            self.flag = np.logical_or(self.flag, flag)