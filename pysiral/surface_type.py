# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 11:25:04 2015

@author: Stefan
"""
import numpy as np
from treedict import TreeDict


class ANDCondition(object):

    def __init__(self):
        self.flag = None

    def add(self, flag):
        if self.flag is None:
            self.flag = flag
        else:
            self.flag = np.logical_and(self.flag, flag)


class TypeContainer(object):

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


class SurfaceType(object):
    """
    Container for surface type information.

    Possible classifications (Adapted from CryoSat-2 conventions)
        - open ocean
        - closed sea/lakes
        - lead
        - large lead/polynya
        - continental ice
        - land
    """
    _VALID_SURFACE_TYPES = [
        "ocean", "closed_sea", "land_ice",
        "sea_ice", "lead", "polynya", "land"]

    def __init__(self):
        self._surface_type_flags = []
        self._n_records = None

    def add_flag(self, flag, type_str):
        """ Add a surface type flag """
        if type_str not in self._VALID_SURFACE_TYPES:
            # TODO: Error Handling
            return
        if self._invalid_n_records(len(flag)):
            # TODO: Error Handling
            return
        setattr(self, type_str, TypeContainer(flag))
        self._surface_type_flags.append(type_str)

    def has_flag(self, type_str):
        return type_str in self._surface_type_flags

    def _invalid_n_records(self, n):
        """ Check if flag array has the correct length """
        if self._n_records is None:  # New flag, ok
            return False
        elif self._n_records == n:   # New flag has correct length
            return False
        else:                        # New flag has wrong length
            return True

    def __getattr__(self, name):
        """
        Return empty lists for surface type flags that have not been
        set yet
        """
        if name in self._VALID_SURFACE_TYPES:
            return TypeContainer(
                np.zeros(shape=(self._n_records), dtype=np.bool))
        raise AttributeError("Unkown surface type: %s" % name)


class IceType(object):
    """
    Container for ice type information

    Possible classifications
        - young thin ice
        - first year ice
        - multi year ice
    """
    _ICE_TYPE_DICT = {
        "thin_ice": 0,
        "fyi": 1,
        "myi": 2}

    def __init__(self):
        self._ice_type_flag = None
