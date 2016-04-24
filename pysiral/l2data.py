# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 16:30:24 2015

@author: Stefan
"""
import numpy as np


class Level2Data(object):

    _L2_DATA_ITEMS = ["mss", "ssa", "elev", "afrb", "rfrb", "range", "sic",
                      "sitype", "snow_depth", "snow_dens", "ice_dens",
                      "sit"]

    _HEMISPHERE_CODES = {"north": "nh", "south": "sh"}

    def __init__(self, l1b):
        # Copy necessary fields form l1b
        self._n_records = l1b.n_records
        self.info = l1b.info
        self.track = l1b.time_orbit
        # Create Level2 Data Groups
        self._create_l2_data_items()

    def set_surface_type(self, surface_type):
        self.surface_type = surface_type

    def update_retracked_range(self, retracker):
        # Update only for indices (surface type) supplied by retracker class
        # XXX: should get an overhaul
        ii = retracker.indices
        self.range[ii] = retracker.range[ii]
        self.elev[ii] = self.track.altitude[ii] - retracker.range[ii]

    def _create_l2_data_items(self):
        for item in self._L2_DATA_ITEMS:
            setattr(self, item, L2ElevationArray(shape=(self._n_records)))

    @property
    def n_records(self):
        return self._n_records

    @property
    def hemisphere(self):
        return self.info.subset_region_name

    @property
    def hemisphere_code(self):
        return self._HEMISPHERE_CODES[self.hemisphere]



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

    def __getslice__(self, i, j):
        r = np.ndarray.__getslice__(self, i, j)
        r.uncertainty = r.uncertainty[i:j]
        r.bias = r.bias[i:j]
        return r

    def set_value(self, value):
#        uncertainty = self.uncertainty
#        bias = self.bias
        self[:] = value[:]
#        setattr(self, "uncertainty", uncertainty)
#        setattr(self, "bias", bias)

    def set_uncertainty(self, uncertainty):
        self.uncertainty = uncertainty
