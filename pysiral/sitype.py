# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 13:57:56 2016

@author: Stefan
"""

from pysiral.auxdata import AuxdataBaseClass
from pysiral.iotools import ReadNC

import scipy.ndimage as ndimage
from pyproj import Proj
import numpy as np
import os


class SITypeBaseClass(AuxdataBaseClass):

    def __init__(self):
        super(SITypeBaseClass, self).__init__()
        self._msg = ""

    def get_along_track_sitype(self, l2):
        sitype, sitype_uncertainty, msg = self._get_along_track_sitype(l2)
        return sitype, sitype_uncertainty, msg


class NoneHandler(SITypeBaseClass):
    """ Dummy handler only returning NaN's """
    def __init__(self):
        super(NoneHandler, self).__init__()

    def _get_along_track_sitype(self, l2):
        sitype = np.full((l2.n_records), np.nan)
        uncertainty = np.full((l2.n_records), np.nan)
        return sitype, uncertainty, ""


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
        self._get_requested_date(l2)
        self._get_data(l2)
        if self.error.status:
            return None, None, self.error.message
        sitype, uncertainty, self._msg = self._get_sitype_track(l2)
        return sitype, uncertainty, self._msg

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
        x_min = x[dim.n_lines-1, 0]-(0.5*dim.dx)
        y_min = y[dim.n_lines-1, 0]-(0.5*dim.dy)
        ix, iy = (l2x-x_min)/dim.dx, (l2y-y_min)/dim.dy

        # Extract along track data from grid
        sitype = ndimage.map_coordinates(
            self._data.ice_type, [iy, ix], order=0)
        # Convert flags to myi fraction
        translator = np.array([np.nan, np.nan, 0.0, 1.0, 0.5, np.nan])
        fillvalues = np.where(sitype == -1)[0]
        sitype[fillvalues] = 5
        sitype = np.array([translator[value] for value in sitype])
        sitype_uncertainty = np.full(sitype.shape, 0.0)
        return sitype, sitype_uncertainty, self._msg


class OsiSafSITypeCDR(SITypeBaseClass):
    """ Class for reprocessed OSISAF sea ice type products (e.g. for C3S).
    Needs to be merged into single OsiSafSitype class at some point """

    def __init__(self):
        super(OsiSafSITypeCDR, self).__init__()

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
        self._get_requested_date(l2)
        self._get_data(l2)
        if self.error.status:
            return None, None, self.error.message
        sitype, uncertainty = self._get_sitype_track(l2)
        return sitype, uncertainty, self._msg

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

        # Read and prepare input data
        self._data = ReadNC(path)
        self._data.ice_type = np.flipud(self._data.ice_type[0, :, :])
        self._data.uncertainty = np.flipud(self._data.uncertainty[0, :, :])

        # Logging
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
        """ Extract ice type and ice type uncertainty along the track """

        # Convert grid/track coordinates to grid projection coordinates
        kwargs = self._options[l2.hemisphere].projection
        p = Proj(**kwargs)
        x, y = p(self._data.lon, self._data.lat)
        l2x, l2y = p(l2.track.longitude, l2.track.latitude)

        # Convert track projection coordinates to image coordinates
        # x: 0 < n_lines; y: 0 < n_cols
        dim = self._options[l2.hemisphere].dimension
        x_min = x[dim.n_lines-1, 0]-(0.5*dim.dx)
        y_min = y[dim.n_lines-1, 0]-(0.5*dim.dy)
        ix, iy = (l2x-x_min)/dim.dx, (l2y-y_min)/dim.dy

        # Extract along track data from grid
        sitype = ndimage.map_coordinates(
            self._data.ice_type, [iy, ix], order=0)
        uncertainty = ndimage.map_coordinates(
            self._data.uncertainty, [iy, ix], order=0)

        # Convert flags to myi fraction
        translator = np.array([np.nan, np.nan, 0.0, 1.0, 0.5, np.nan])
        fillvalues = np.where(sitype == -1)[0]
        sitype[fillvalues] = 5
        sitype = np.array([translator[value] for value in sitype])

        # Uncertainty in product is in %
        sitype_uncertainty = uncertainty / 100.

#        import matplotlib.pyplot as plt
#
#        plt.figure(dpi=150)
#        plt.imshow(self._data.ice_type, interpolation="none")
#        plt.plot(ix, iy)
#        plt.scatter(ix[0], iy[0])
#
#        plt.figure("projection coordinates")
#        plt.plot(l2x, l2y)
#
#        plt.figure("sitype along track")
#        plt.plot(sitype)
#        plt.plot(sitype_uncertainty)
#
#        plt.show()
#
#        stop

        return sitype, sitype_uncertainty


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
        self._get_requested_date(l2)
        self._get_data(l2)
        if self.error.status:
            return None, None, self.error.message
        sitype, sitype_uncertainty = self._get_sitype_track(l2)
        return sitype, sitype_uncertainty, self._msg

    def _get_requested_date(self, l2):
        """ Use first timestamp as reference, date changes are ignored """
        year = l2.track.timestamp[0].year
        month = l2.track.timestamp[0].month
        day = l2.track.timestamp[0].day
        self._requested_date = [year, month, day]

    def _get_data(self, l2):
        """ Loads file from local repository only if needed """

        opt = self._options

        # Check if file is already loaded
        if self._requested_date == self._current_date:
            # Data already loaded, nothing to do
            self._msg = "ICDCNasaTeam: Daily grid already present"
            return

        # construct filename
        path = self._get_local_repository_filename(l2)

        # Check if the file exists, add an error if not
        # (error is not raised at this point)
        if not os.path.isfile(path):
            self._msg = "ICDCNasaTeam: File not found: %s " % path
            self.error.add_error("auxdata_missing_sitype", self._msg)
            return

        # Bulk read the netcdf file
        self._data = ReadNC(path)

        # There are multiple myi concentrations fields in the product
        # The one used here is defined in the auxdata definition file
        # in the pysiral config folder (`auxdata_def.yaml`)
        # -> root.sitype.icdc_nasateam.options.variable_name
        myi_fraction = getattr(self._data, opt.variable_name)
        self._data.ice_type = myi_fraction[0, :, :]

        # Same for the uncertainty variable
        # (see description directly above for how to access variable namde
        #  definition)
        myi_fraction_unc = getattr(self._data, opt.uncertainty_variable_name)
        self._data.ice_type_uncertainty = myi_fraction_unc[0, :, :]

        # Report and save current data period
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

        x_min = x[0, 0]-(0.5*dim.dx)
        y_min = y[0, 0]-(0.5*dim.dy)
        ix, iy = (l2x-x_min)/dim.dx, (l2y-y_min)/dim.dy

        # Extract along track data from grid
        myi_concentration_percent = ndimage.map_coordinates(
            self._data.ice_type, [iy, ix], order=0)

        myi_concentration_uncertainty = ndimage.map_coordinates(
            self._data.ice_type_uncertainty, [iy, ix], order=0)

        # Convert percent [0-100] into fraction [0-1]
        sitype = myi_concentration_percent/100.
        sitype_uncertainty = myi_concentration_uncertainty/100.

        # Remove invalid (negative) values
        sitype[np.where(sitype < 0)] = 0.0
        sitype_uncertainty[np.where(sitype_uncertainty < 0)] = 0.0

        return sitype, sitype_uncertainty


class MYIDefault(SITypeBaseClass):
    """ Returns myi for all ice covered regions """

    def __init__(self):
        super(MYIDefault, self).__init__()

    def _get_along_track_sitype(self, l2):
        """ Every ice is myi (sitype = 1) """
        sitype = np.zeros(shape=l2.sic.shape, dtype=np.float32)
        is_ice = np.where(l2.sic > 0)[0]
        sitype[is_ice] = 1.0
        sitype_uncertainty = np.full(sitype.shape, self.uncertainty_default)
        return sitype, sitype_uncertainty, ""

    @property
    def uncertainty_default(self):
        if "uncertainty_default" in self._options:
            return self._options.uncertainty_default
        else:
            return 0.0


class FYIDefault(SITypeBaseClass):
    """ Returns myi for all ice covered regions """

    def __init__(self):
        super(FYIDefault, self).__init__()

    def _get_along_track_sitype(self, l2):
        """ Every ice is fyi (sitype = 0) """
        sitype = np.zeros(shape=l2.sic.shape, dtype=np.float32)
        sitype_uncertainty = np.full(sitype.shape, self.uncertainty_default)
        return sitype, sitype_uncertainty, ""

    @property
    def uncertainty_default(self):
        try:
            if "uncertainty_default" in self._options:
                return self._options.uncertainty_default
        except TypeError:
            pass

        return 0.0


def get_l2_sitype_handler(name):
    pyclass = globals().get(name, None)
    if pyclass is not None:
        return pyclass()
    else:
        return pyclass
