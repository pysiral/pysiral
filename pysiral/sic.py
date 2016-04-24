# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 13:57:56 2016

@author: Stefan
"""
from pysiral.config import options_from_dictionary
from pysiral.iotools import ReadNC

import scipy.ndimage as ndimage
from pyproj import Proj
import numpy as np
import os


class SICBaseClass(object):

    def __init__(self):
        self._options = None
        self._local_repository = None
        self._subfolders = []
        self._msg = ""

    def set_options(self, **opt_dict):
        self._options = options_from_dictionary(**opt_dict)

    def set_local_repository(self, path):
        self._local_repository = path

    def set_filenaming(self, filenaming):
        self._filenaming = filenaming

    def set_subfolders(self, subfolder_list):
        self._subfolders = subfolder_list

    def get_along_track_sic(self, l2):
        sic, msg = self._get_along_track_sic(l2)
        return sic, msg


class OsiSafSIC(SICBaseClass):

    def __init__(self):
        super(OsiSafSIC, self).__init__()
        self._data = None
        self._current_date = [0, 0, 0]
        self._requested_date = [-1, -1, -1]

    @property
    def year(self):
        return "%04g" % self._requested_date[0]

    @property
    def month(self):
        return "%02g" % self._requested_date[1]

    @property
    def day(self):
        return "%02g" % self._requested_date[2]

    def _get_along_track_sic(self, l2):
        self._msg = ""
        self._get_requested_date(l2)
        self._get_data(l2)
        sic = self._get_sic_track(l2)
        return sic, self._msg

    def _get_requested_date(self, l2):
        """ Use first timestamp as reference, date changes are ignored """
        year = l2.track.timestamp[0].year
        month = l2.track.timestamp[0].month
        day = l2.track.timestamp[0].day
        self._requested_date = [year, month, day]

    def _get_data(self, l2):
        """ Loads file from local repository only if needed """
        if self._requested_date == self._current_date:
            # Data already loaded, nothing to do
            return
        path = self._get_local_repository_filename(l2)
        self._data = ReadNC(path)
        self._data.ice_conc = self._data.ice_conc[0, :, :]
        flagged = np.where(self._data.ice_conc < 0)
        self._data.ice_conc[flagged] = 0
        # This step is important for calculation of image coordinates
        self._data.ice_conc = np.flipud(self._data.ice_conc)
        self._msg = "Loaded SIC file: %s" % path
        self._current_date = self._requested_date

    def _get_local_repository_filename(self, l2):
        path = self._local_repository
        for subfolder_tag in self._subfolders:
            subfolder = getattr(self, subfolder_tag)
            path = os.path.join(path, subfolder)
        filename = self._filenaming.format(
            year=self.year, month=self.month, day=self.day,
            hemisphere_code=l2.hemisphere_code)
        path = os.path.join(path, filename)
        return path

    def _get_sic_track(self, l2):
        # Convert grid/track coordinates to grid projection coordinates
        kwargs = self._options[l2.hemisphere].projection
        p = Proj(**kwargs)
        x, y = p(self._data.lon, self._data.lat)
        l2x, l2y = p(l2.track.longitude, l2.track.latitude)
        # Convert track projection coordinates to image coordinates
        # x: 0 < n_lines; y: 0 < n_cols
        dim = self._options[l2.hemisphere].dimension
        x_min = x[dim.n_lines-1, 0]
        y_min = y[dim.n_lines-1, 0]
        ix, iy = (l2x-x_min)/dim.dx, (l2y-y_min)/dim.dy
        # Extract along track data from grid
        sic = ndimage.map_coordinates(self._data.ice_conc, [iy, ix], order=0)
        return sic

        # XXX: Debug stuff below

#        import matplotlib.pyplot as plt
#        from mpl_toolkits.basemap import Basemap
#
#        plt.figure()
#        m = Basemap(projection="ortho", lon_0=-45, lat_0=90)
#        m.drawcoastlines()
#        tx, ty = m(l2.track.longitude, l2.track.latitude)
#        m.plot(tx, ty)
#
#        plt.figure("x")
#        plt.imshow(x)
#
#        plt.figure("y")
#        plt.imshow(y)
#
#        plt.figure()

#        plt.imshow(self._data.ice_conc, cmap=plt.get_cmap("viridis"),
#                   vmin=0, vmax=100, extent=extent, origin="upper",
#                   interpolation="none")
#        plt.plot(ix, iy)

#        plt.figure("sic")
#        plt.plot(sic)
#        plt.show()
#        stop


def get_l2_sic_handler(name):
    return globals()[name]()
