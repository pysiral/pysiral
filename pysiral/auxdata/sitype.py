# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 13:57:56 2016

@author: Stefan

Important Note:

    All sitype data handlers must be subclasses of pysiral.auxdata.AuxdataBaseClass in order to work
    for the Level-2 Processor. If the auxiliary class is based on a static dataset, this should be parsed
    in `__init__`.

    Please review the variables and properties in the parent class, as well as the correspodning config and
    support classes for grid track interpolation in the pysiral.auxdata module for additional guidance.

    The only other hard requirements is the presence of on specific method in order to be a valid subclass of
    AuxdataBaseClass:


        get_l2_track_vars(l2)

            This method will be called during the Level-2 processor. The argument is the Level-2 data object and
            the purpose of the method is to compute the auxilary variable(s) and associated uncertainty. These
            variable need to be registered using the `register_auxvar(id, name, value, uncertainty)` method of
            the base class. All sitype subclasses need to register at minimum the following variable:

                sea ice type (fraction of multi year ice):
                    id: sitype
                    name: sea_ice_type

            e.g., this code line is mandatory for `get_l2_track_vars` (uncertainty can be None):

                # Register Variables
                self.register_auxvar("sitype", "sea_ice_type", value, uncertainty)

"""

from pysiral.auxdata import AuxdataBaseClass, GridTrackInterpol
from pysiral.iotools import ReadNC

import scipy.ndimage as ndimage
from pyproj import Proj
import numpy as np
import os


class OsiSafSIType(AuxdataBaseClass):
    """ This is a class for the OSI-403 product with variables ice_type and confidence_level """

    def __init__(self, *args, **kwargs):
        super(OsiSafSIType, self).__init__(*args, **kwargs)
        self._data = None

    def get_l2_track_vars(self, l2):
        """ Default grid auxiliary data set"""

        # These properties are needed to construct the product path
        self.start_time = l2.info.start_time
        self.hemisphere_code = l2.hemisphere_code

        # Set the requested data
        self.set_requested_date_from_l2(l2)

        # Update the external data
        self.update_external_data()

        # Check if data is available
        if self.error.status:
            self.error.raise_on_error()

        # Get the data
        sitype, uncertainty = self._get_sitype_track(l2)

        # Register the data
        self.register_auxvar("sitype", "sea_ice_type", sitype, uncertainty)

    def load_requested_auxdata(self):
        """ Required subclass method: Load the data file necessary to satisfy condition for requested date"""

        # Retrieve the file path for the requested date from a property of the auxdata parent class
        path = self.requested_filepath

        # Validation
        if not os.path.isfile(path):
            msg = "OsiSafSIType: File not found: %s " % path
            self.add_handler_message(msg)
            self.error.add_error("auxdata_missing_sitype", msg)
            return

        # --- Read the data ---
        self._data = ReadNC(path)

        # Report
        self.add_handler_message("OsiSafSIType: Loaded SIType file: %s" % path)

    def _get_sitype_track(self, l2):

        # Extract from grid
        griddef = self.cfg.options[l2.hemisphere]
        grid_lons, grid_lats = self._data.lon, self._data.lat
        grid2track = GridTrackInterpol(l2.track.longitude, l2.track.latitude, grid_lons, grid_lats, griddef)
        sitype = grid2track.get_from_grid_variable(self._data.ice_type[0, :, :], flipud=True)
        confidence_level = grid2track.get_from_grid_variable(self._data.confidence_level[0, :, :], flipud=True)

        # set fill values to flag 0 -> nan
        fillvalues = np.where(sitype == -1)[0]
        sitype[fillvalues] = 0

        # --- Translate sitype codes into myi fraction ---
        # flag_meanings: -1: fill value, 1: open_water, 2: first_year_ice, 3: multi_year_ice, 4: ambiguous
        translator = np.array([np.nan, np.nan, 0.0, 1.0, 0.5])
        sitype = np.array([translator[value] for value in sitype])

        # Translate confidence level into myi fraction uncertainty
        # flag_meaning: 0: unprocessed, 1: erroneous, 2: unreliable, 3: acceptable, 4: good, 5: excellent
        translator = np.array([np.nan, 1., 0.5, 0.2, 0.1, 0.0])
        sitype_uncertainty = np.array([translator[value] for value in confidence_level])

        return sitype, sitype_uncertainty

    @property
    def requested_filepath(self):
        """ Note: this overwrites the property in the super class due to some
        peculiarities with the filenaming (hemisphere code) """
        path = self.cfg.local_repository
        for subfolder_tag in self.cfg.subfolders:
            subfolder = getattr(self, subfolder_tag)
            path = os.path.join(path, subfolder)
        filename = self.cfg.filenaming.format(
            year=self.year, month=self.month, day=self.day,
            hemisphere_code=self.hemisphere_code)
        path = os.path.join(path, filename)
        return path


class OsiSafSITypeCDR(AuxdataBaseClass):
    """ Class for reprocessed OSISAF sea ice type products (e.g. for C3S).
    Needs to be merged into single OsiSafSitype class at some point """

    def __init__(self, *args, **kwargs):
        super(OsiSafSITypeCDR, self).__init__(*args, **kwargs)
        self._data = None

    def get_l2_track_vars(self, l2):
        """ Mandadory method of AuxdataBaseClass subclass """

        # These properties are needed to construct the product path
        self.start_time = l2.info.start_time
        self.hemisphere_code = l2.hemisphere_code

        # Set the requested data
        self.set_requested_date_from_l2(l2)

        # Update the external data
        self.update_external_data()

        # Check if error with file I/O
        if self.error.status or self._data is None:
            sitype = self.get_empty_array(l2)
            uncertainty = self.get_empty_array(l2)
        else:
            # Get and return the track
            sitype, uncertainty = self._get_sitype_track(l2)

        # Register the data
        self.register_auxvar("sitype", "sea_ice_type", sitype, uncertainty)

    def load_requested_auxdata(self):
        """ Loads file from local repository only if needed """

        # Retrieve the file path for the requested date from a property of the auxdata parent class
        path = self.requested_filepath

        # Validation
        if not os.path.isfile(path):
            msg = "%s: File not found: %s " % (self.__class__.__name__, path)
            self.add_handler_message(msg)
            self.error.add_error("auxdata_missing_sitype", msg)
            return

        # Read and prepare input data
        self._data = ReadNC(path)
        self._data.ice_type = np.flipud(self._data.ice_type[0, :, :])
        self._data.uncertainty = np.flipud(self._data.uncertainty[0, :, :])

    def _get_sitype_track(self, l2):
        """ Extract ice type and ice type uncertainty along the track """

        # Extract from grid
        griddef = self.cfg.options[l2.hemisphere]
        grid_lons, grid_lats = self._data.lon, self._data.lat
        grid2track = GridTrackInterpol(l2.track.longitude, l2.track.latitude, grid_lons, grid_lats, griddef)
        sitype = grid2track.get_from_grid_variable(self._data.ice_type, flipud=False)
        uncertainty = grid2track.get_from_grid_variable(self._data.uncertainty, flipud=False)
        # sic = grid2track.get_from_grid_variable(self._data.ice_conc, flipud=True)
        #
        # # Convert grid/track coordinates to grid projection coordinates
        # kwargs = self.cfg.options[l2.hemisphere].projection
        # p = Proj(**kwargs)
        # x, y = p(self._data.lon, self._data.lat)
        # l2x, l2y = p(l2.track.longitude, l2.track.latitude)
        #
        # # Convert track projection coordinates to image coordinates
        # # x: 0 < n_lines; y: 0 < n_cols
        # dim = self.cfg.options[l2.hemisphere].dimension
        # x_min = x[dim.n_lines-1, 0]-(0.5*dim.dx)
        # y_min = y[dim.n_lines-1, 0]-(0.5*dim.dy)
        # ix, iy = (l2x-x_min)/dim.dx, (l2y-y_min)/dim.dy
        #
        # # Extract along track data from grid
        # sitype = ndimage.map_coordinates(self._data.ice_type, [iy, ix], order=0)
        # uncertainty = ndimage.map_coordinates(self._data.uncertainty, [iy, ix], order=0)

        # Convert flags to myi fraction
        translator = np.array([np.nan, np.nan, 0.0, 1.0, 0.5, np.nan])
        fillvalues = np.where(sitype == -1)[0]
        sitype[fillvalues] = 5
        sitype = np.array([translator[value] for value in sitype])

        # Uncertainty in product is in %
        sitype_uncertainty = uncertainty / 100.

        return sitype, sitype_uncertainty

    @property
    def requested_filepath(self):
        """ Note: this overwrites the property in the super class due to some
        peculiarities with the filenaming (auto product changes etc) """

        # Unique to this class is the possibility to auto merge
        # products. The current implementation supports only two products
        path = self.cfg.local_repository

        # The path needs to be completed if two products shall be used
        if self.cfg.options.has_key("auto_product_change"):
            opt = self.cfg.options.auto_product_change
            product_index = int(self.start_time > opt.date_product_change)
            product_def = opt.osisaf_product_def[product_index]
            path = os.path.join(path, product_def["subfolder"])
            self.cfg.filenaming = product_def["filenaming"]
            self.cfg.long_name = product_def["long_name"]

        for subfolder_tag in self.cfg.subfolders:
            subfolder = getattr(self, subfolder_tag)
            path = os.path.join(path, subfolder)

        filename = self.cfg.filenaming.format(
            year=self.year, month=self.month, day=self.day,
            hemisphere_code=self.hemisphere_code)
        path = os.path.join(path, filename)
        return path


class ICDCNasaTeam(AuxdataBaseClass):
    """ MYI Fraction from NASA Team Algorithm (from ICDC UHH) """

    def __init__(self, *args, **kwargs):
        super(ICDCNasaTeam, self).__init__(*args, **kwargs)
        self._data = None

    def get_l2_track_vars(self, l2):
        self._get_requested_date(l2)
        self._get_data(l2)
        if self.error.status:
            sitype, sitype_uncertainty = self.get_empty_val(l2), self.get_empty_val(l2)
        else:
            sitype, sitype_uncertainty = self._get_sitype_track(l2)
        self.register_auxvar("sitype", "sea_ice_type", sitype, None)

    def _get_requested_date(self, l2):
        """ Use first timestamp as reference, date changes are ignored """
        year = l2.track.timestamp[0].year
        month = l2.track.timestamp[0].month
        day = l2.track.timestamp[0].day
        self._requested_date = [year, month, day]

    def _get_data(self, l2):
        """ Loads file from local repository only if needed """

        opt = self.cfg.options

        # Check if file is already loaded
        if self._requested_date == self._current_date:
            # Data already loaded, nothing to do
            self.add_handler_message("ICDCNasaTeam: Daily grid already present")
            return

        # construct filename
        path = self._get_local_repository_filename(l2)

        # Check if the file exists, add an error if not
        # (error is not raised at this point)
        if not os.path.isfile(path):
            msg = "ICDCNasaTeam: File not found: %s " % path
            self.add_handler_message(msg)
            self.error.add_error("auxdata_missing_sitype", msg)
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
        self.add_handler_message("ICDCNasaTeam: Loaded SIType file: %s" % path)
        self._current_date = self._requested_date

    def _get_local_repository_filename(self, l2):
        path = self.cfg.local_repository
        for subfolder_tag in self.cfg.subfolders:
            subfolder = getattr(self, subfolder_tag)
            path = os.path.join(path, subfolder)
        filename = self.cfg.filenaming.format(
            year=self.year, month=self.month, day=self.day,
            hemisphere_code=l2.hemisphere_code)
        path = os.path.join(path, filename)
        return path

    def _get_sitype_track(self, l2):

        # Convert grid/track coordinates to grid projection coordinates
        kwargs = self.cfg.option[l2.hemisphere].projection
        p = Proj(**kwargs)
        x, y = p(self._data.longitude, self._data.latitude)
        l2x, l2y = p(l2.track.longitude, l2.track.latitude)

        # Convert track projection coordinates to image coordinates
        # x: 0 < n_lines; y: 0 < n_cols
        dim = self.cfg.option[l2.hemisphere].dimension

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


class MYIDefault(AuxdataBaseClass):
    """ Returns myi for all ice covered regions """

    def __init__(self, *args, **kwargs):
        super(MYIDefault, self).__init__(*args, **kwargs)

    def get_l2_track_vars(self, l2):
        """ Every ice is myi (sitype = 1) """
        sitype = np.zeros(shape=l2.sic.shape, dtype=np.float32)
        is_ice = np.where(l2.sic > 0)[0]
        sitype[is_ice] = 1.0
        uncertainty = np.full(sitype.shape, self.uncertainty_default)
        # Register the data
        self.register_auxvar("sitype", "sea_ice_type", sitype, uncertainty)

    @property
    def uncertainty_default(self):
        if "uncertainty_default" in self.cfg.option:
            return self.cfg.option.uncertainty_default
        else:
            return 0.0


class FYIDefault(AuxdataBaseClass):
    """ Returns myi for all ice covered regions """

    def __init__(self, *args, **kwargs):
        super(FYIDefault, self).__init__(*args, **kwargs)

    def get_l2_track_vars(self, l2):
        """ Every ice is fyi (sitype = 0) """
        sitype = np.zeros(shape=l2.sic.shape, dtype=np.float32)
        uncertainty = np.full(sitype.shape, self.uncertainty_default)
        self.register_auxvar("sitype", "sea_ice_type", sitype, uncertainty)

    @property
    def uncertainty_default(self):
        try:
            if "uncertainty_default" in self.cfg.options:
                return self.cfg.options.uncertainty_default
        except TypeError:
            pass
        return 0.0
