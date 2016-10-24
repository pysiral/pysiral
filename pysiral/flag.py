# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 15:44:22 2016

@author: Stefan
"""

import numpy as np


class FlagContainer(object):

    def __init__(self, flag):
        self._flag = flag

    def set_flag(self, flag):
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


class ANDCondition(FlagContainer):

    def __init__(self):
        super(ANDCondition, self).__init__(None)

    def add(self, flag):
        if self._flag is None:
            self.set_flag(flag)
        else:
            self.set_flag(np.logical_and(self.flag, flag))


class ORCondition(FlagContainer):

    def __init__(self):
        super(ORCondition, self).__init__(None)

    def add(self, flag):
        if self._flag is None:
            self.set_flag(flag)
        else:
            self.set_flag(np.logical_or(self.flag, flag))
