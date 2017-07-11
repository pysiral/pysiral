# -*- coding: utf-8 -*-
"""
Created on Thu Apr 07 16:19:41 2016

@author: shendric
"""

import numpy as np


class GridMapParameterBase(object):
    """
    Contains data, pcolor grid calculation capability, colormap definition
    and standardized parameter naming
    """

    def __init__(self):
        from pysiral.config import get_parameter_definitions
        self._parameter_definitions = get_parameter_definitions()
        self._projection = None
        self.latitude = None
        self.longitude = None
        self.grid = None
        self.pgrid = None
        self.pardef = None

    def set_grid(self, longitude, latitude):
        self.longitude = longitude
        self.latitude = latitude

    def set_projection(self, **projection):
        from pysiral.maptools import GeoPcolorGrid
        self.pgrid = GeoPcolorGrid(self.longitude, self.latitude)
        self.pgrid.calc_from_proj(**projection)

    def set_nan_mask(self, nan_mask):
        mask_indices = np.where(np.isnan(nan_mask))
        self.grid.mask[mask_indices] = True

    def set_cmap(self, label=None, unit=None, short_name=None, cmap=None,
                 cmap_diff=None):
        from pysiral.config import options_from_dictionary
        if label is not None:
            self.pardef.label = label
        if unit is not None:
            self.pardef.unit = unit
        if short_name is not None:
            self.pardef.short_name = short_name
        if cmap is not None:
            self.pardef.cmap = options_from_dictionary(**cmap)
        if cmap_diff is not None:
            self.pardef.cmap_diff = options_from_dictionary(**cmap_diff)


class GridMapParameter(GridMapParameterBase):

    def __init__(self):
        super(GridMapParameter, self).__init__()

    def get_cmap(self):
        return self.pardef.cmap

    def get_label(self):
        return self.pardef.label+" ("+self.pardef.unit+")"

    def set_parameter(self, grid, parameter_name):
        self.grid = grid
        self.pardef = self._parameter_definitions[parameter_name]

    @property
    def short_name(self):
        return self.pardef.short_name


class GridMapDiffParameter(GridMapParameterBase):
    """
    Contains data, pcolor grid calculation capability, colormap definition
    and standardized parameter naming
    """

    def __init__(self):
        super(GridMapDiffParameter, self).__init__()

    def get_cmap(self):
        return self.pardef.cmap_diff

    def get_label(self):
        return "$\Delta$ "+self.pardef.label+" ("+self.pardef.unit+")"

    def set_parameter(self, grida, gridb, parameter_name):
        self.grid = grida-gridb
        self.pardef = self._parameter_definitions[parameter_name]

    @property
    def short_name(self):
        return self.pardef.short_name + "_diff"


class GridMapAddParameter(GridMapParameterBase):
    """
    Contains data, pcolor grid calculation capability, colormap definition
    and standardized parameter naming
    """

    def __init__(self):
        super(GridMapAddParameter, self).__init__()

    def get_cmap(self):
        return self.pardef.cmap

    def get_label(self):
        return "$\Delta$ "+self.pardef.label+" ("+self.pardef.unit+")"

    def set_parameter(self, grida, gridb, parameter_name):
        self.grid = gridb+grida
        self.pardef = self._parameter_definitions[parameter_name]

    @property
    def short_name(self):
        return self.pardef.short_name + "_diff"


class GridMapAverageParameter(GridMapParameterBase):
    """
    Contains data, pcolor grid calculation capability, colormap definition
    and standardized parameter naming
    """

    def __init__(self):
        super(GridMapAverageParameter, self).__init__()

    def get_cmap(self):
        return self.pardef.cmap

    def get_label(self):
        return "$\Delta$ "+self.pardef.label+" ("+self.pardef.unit+")"

    def set_parameter(self, grid_stack, parameter_name):
        stop
        self.pardef = self._parameter_definitions[parameter_name]

    @property
    def short_name(self):
        return self.pardef.short_name + "_average"

