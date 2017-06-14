# -*- coding: utf-8 -*-
"""
Created on Sun Jun 11 19:24:04 2017

@author: Stefan
"""

from pysiral.config import get_yaml_config
from pysiral.errorhandler import ErrorStatus
from pysiral.logging import DefaultLoggingClass

from treedict import TreeDict
from pyproj import Proj

import numpy as np


class GridDefinition(DefaultLoggingClass):

    def __init__(self, preset=None):
        super(GridDefinition, self).__init__(self.__class__.__name__)
        self.error = ErrorStatus(caller_id=self.__class__.__name__)
        self._preset = preset
        self._metadata = {"grid_id": "n/a", "grid_tag": "n/a",
                          "hemisphere": "n/a", "resolution_tag": "n/a"}
        self._proj = None
        self._proj_dict = {}
        self._extent_dict = {}

    def set_from_griddef_file(self, filename):
        config = get_yaml_config(filename)
        for key in self._metadata.keys():
            self._metadata[key] = config[key]
        self.set_projection(**config.projection)
        self.set_extent(**config.extent)

    def set_projection(self, **kwargs):
        self._proj_dict = kwargs
        self._set_proj()

    def proj(self, longitude, latitude, **kwargs):
        projx, projy = self._proj(longitude, latitude, **kwargs)
        return projx, projy

    def grid_indices(self, longitude, latitude):
        """ Computes the grid indices the given lon/lat pairs would be sorted
        into (no clipping) """
        projx, projy = self.proj(longitude, latitude)
        extent = self.extent
        xi = np.floor((projx + extent.xsize/2.0)/extent.dx)
        yj = np.floor((projy + extent.ysize/2.0)/extent.dy)
        return xi, yj

    def get_grid_coordinates(self, mode="center"):
        """ Returns longitude/latitude points for each grid cell
        Note: mode keyword only for future use. center coordinates are
        returned by default """
        x0, y0 = self.extent.xoff, self.extent.yoff
        xsize, ysize = self.extent.xsize, self.extent.ysize
        numx, numy = self.extent.numx, self.extent.numy
        xmin, xmax = x0-(xsize/2.), x0+(xsize/2.)
        ymin, ymax = y0-ysize/2., y0+ysize/2.
        x = np.linspace(xmin, xmax, num=numx)
        y = np.linspace(ymin, ymax, num=numy)
        xx, yy = np.meshgrid(x, y)
        lon, lat = self.proj(xx, yy, inverse=True)
        return lon, lat

    def set_extent(self, **kwargs):
        self._extent_dict = kwargs

    def _set_proj(self):
        self._proj = Proj(**self._proj_dict)

    @property
    def hemisphere(self):
        return self._metadata["hemisphere"]

    @property
    def grid_id(self):
        return self._metadata["grid_id"]

    @property
    def grid_tag(self):
        return self._metadata["grid_tag"]

    @property
    def resolution_tag(self):
        return self._metadata["resolution_tag"]

    @property
    def extent(self):
        return TreeDict.fromdict(self._extent_dict, expand_nested=True)
