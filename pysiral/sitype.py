# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 13:57:56 2016

@author: Stefan
"""

from pysiral.errorhandler import ErrorStatus
from pysiral.config import options_from_dictionary
from pysiral.iotools import ReadNC

import scipy.ndimage as ndimage
from pyproj import Proj
import numpy as np
import os


class SITypeBaseClass(object):

    def __init__(self):
        self._options = None
        self._local_repository = None
        self._subfolders = []
        self.error = ErrorStatus()

    def set_options(self, **opt_dict):
        self._options = options_from_dictionary(**opt_dict)

    def set_local_repository(self, path):
        self._local_repository = path

    def set_filenaming(self, filenaming):
        self._filenaming = filenaming

    def set_subfolders(self, subfolder_list):
        self._subfolders = subfolder_list

    def get_along_track_sitype(self, l2):
        sitype, msg = self._get_along_track_sitype(l2)
        return sitype, msg


class OsiSafSIType(SITypeBaseClass):

    def __init__(self):
        super(OsiSafSIType, self).__init__()
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

    def _get_along_track_sitype(self, l2):
        self._msg = ""
        self._get_requested_date(l2)
        self._get_data(l2)
        sic = self._get_sitype_track(l2)
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

        # construct filename
        path = self._get_local_repository_filename(l2)

        # Validation
        if not os.path.isfile(path):
            msg = "OsiSafSIType: File not found: %s " % path
            self.error.add_error("auxdata_missing_sitype", msg)
            return

        self._data = ReadNC(path)
        self._data.ice_type = self._data.ice_type[0, :, :]
        # self._data.confidence_level = self._data.confidence_level[0, :, :]

        # This step is important for calculation of image coordinates
        self._data.ice_type = np.flipud(self._data.ice_type)
        # self._data.confidence_level = np.flipud(self._data.confidence_level)
        self._msg = "OsiSafSIType: Loaded SIType file: %s" % path
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

    def _get_sitype_track(self, l2):
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
        sitype = ndimage.map_coordinates(
            self._data.ice_type, [iy, ix], order=0)
        # Convert flags to myi fraction
        translator = np.array([np.nan, np.nan, 0.0, 1.0, 0.5, np.nan])
        fillvalues = np.where(sitype == -1)[0]
        sitype[fillvalues] = 5
        sitype = np.array([translator[value] for value in sitype])
        return sitype


class ICDCNasaTeam(SITypeBaseClass):
    """ MYI Fraction from NASA Team Algorithm (from ICDC UHH) """

    def __init__(self):
        super(ICDCNasaTeam, self).__init__()
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

    def _get_along_track_sitype(self, l2):
        self._msg = ""
        self._get_requested_date(l2)
        self._get_data(l2)
        sic = self._get_sitype_track(l2)
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

        # construct filename
        path = self._get_local_repository_filename(l2)

        # Validation
        if not os.path.isfile(path):
            msg = "ICDCNasaTeam: File not found: %s " % path
            self.error.add_error("auxdata_missing_sitype", msg)
            return

        self._data = ReadNC(path)
        myi_fraction = getattr(self._data, self._options.variable_name)
        self._data.ice_type = myi_fraction[0, :, :]

        # This step is important for calculation of image coordinates
        # self._data.ice_type = np.flipud(self._data.ice_type)
        # self._data.confidence_level = np.flipud(self._data.confidence_level)
        self._msg = "ICDCNasaTeam: Loaded SIType file: %s" % path
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

    def _get_sitype_track(self, l2):

        # Convert grid/track coordinates to grid projection coordinates
        kwargs = self._options[l2.hemisphere].projection
        p = Proj(**kwargs)
        x, y = p(self._data.longitude, self._data.latitude)
        l2x, l2y = p(l2.track.longitude, l2.track.latitude)

        # Convert track projection coordinates to image coordinates
        # x: 0 < n_lines; y: 0 < n_cols
        dim = self._options[l2.hemisphere].dimension
        x_min = x[0, 0]
        y_min = y[0, 0]
        ix, iy = (l2x-x_min)/dim.dx, (l2y-y_min)/dim.dy

        # Extract along track data from grid
        myi_fraction_percent = ndimage.map_coordinates(
            self._data.ice_type, [iy, ix], order=0)

        # Convert percent [0-100] into fraction [0-1]
        sitype = myi_fraction_percent/100.

        # Remove invalid parameter
        invalid = np.where(sitype < 0)[0]
        sitype[invalid] = 0.0

        return sitype


class MYIDefault(SITypeBaseClass):
    """ Returns myi for all ice covered regions """

    def __init__(self):
        super(MYIDefault, self).__init__()

    def _get_along_track_sitype(self, l2):
        """ Every ice is myi (sitype = 1) """
        sitype = np.zeros(shape=l2.sic.shape, dtype=np.float32)
        is_ice = np.where(l2.sic > 0)[0]
        sitype[is_ice] = 1.0
        return sitype, ""


def get_l2_sitype_handler(name):
    return globals()[name]()
