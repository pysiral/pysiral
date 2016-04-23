# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 15:30:53 2016

@author: Stefan
"""
from pysiral.config import options_from_dictionary
from pysiral.flag import FlagContainer, ORCondition


class FilterBaseClass(object):

    def __init__(self):
        self._flag = None

    def set_options(self, **opt_dict):
        self._options = options_from_dictionary(**opt_dict)

    def apply_filter(self, *args, **kwargs):
        self._apply_filter(*args, **kwargs)

    @property
    def flag(self):
        return self._flag


# %% Filter for Level2Processor

class FreeboardValidRange(FilterBaseClass):
    """
    Filters freeboard outliers by simple min/max thresholding
    Requires l2 data container and target (either: "afrb", "rfrb")
    """

    def __init__(self):
        super(FreeboardValidRange, self).__init__()

    def _apply_filter(self, l2, target):
        freeboard = getattr(l2, target)
        invalid = ORCondition()
        invalid.add(freeboard < self._options.valid_minimum_point_value)
        invalid.add(freeboard > self._options.valid_maximum_point_value)
        self._flag = FlagContainer(invalid.flag)


def get_filter(name):
    return globals()[name]()
