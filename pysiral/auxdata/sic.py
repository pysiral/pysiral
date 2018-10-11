# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 13:57:56 2016

@author: Stefan

Important Note:

    All sic data handlers must be subclasses of pysiral.auxdata.AuxdataBaseClass in order to work
    for the Level-2 Processor. By convention, two methods must exist to make the subclass valid.

        subclass_init()

            This method will be called by default at the beginning of the Level-2 processor *after* all
            options are passed to the instance. It can be used for e.g. reading static data sets

        get_l2_track_vars(l2)

            This method will be called during the Level-2 processor. The argument is the Level-2 data object and
            the purpose of the method is to compute the auxilary variable(s) and associated uncertainty. These
            variable need to be registered using the `register_auxvar(id, name, value, uncertainty)` method of
            the base class. All sic subclasses need to register at minimum the following variable:

                sea ice concentration (in percent):
                    id: sic
                    name: sea_ice_concentration

            e.g., this code line is mandatory for `get_l2_track_vars` (uncertainty can be None):

                # Register Variables
                self.register_auxvar("sic", "sea_ice_concentration", value, uncertainty)

"""


from pysiral.auxdata import AuxdataBaseClass, GridTrackInterpol
from pysiral.iotools import ReadNC

import scipy.ndimage as ndimage
from pyproj import Proj
import numpy as np
import os


class OsiSafSIC(AuxdataBaseClass):


    def __init__(self):
        super(OsiSafSIC, self).__init__()
        self._data = None
        self.error.caller_id = self.__class__.__name__

    def subclass_init(self):
        pass

    def get_l2_track_vars(self, l2):
        """ Main entry point of the class """

        # These properties are needed to construct the product path
        self.start_time = l2.info.start_time
        self.hemisphere_code = l2.hemisphere_code

        # Set the requested date
        self.set_requested_date_from_l2(l2)

        # Update the external data
        self.update_external_data()

        # Check if error with file I/O
        if self.error.status:
            sic = self.get_empty_array(l2)
        else:
            # Get and return the track
            sic = self._get_sic_track(l2)

            # Fill pole hole
            if hasattr(self.options, "fill_pole_hole"):
                opt = self.options.fill_pole_hole
                is_near_pole_hole = l2.track.latitude >= opt.pole_hole_lat_threshold
                is_nan = np.isnan(sic)
                indices = np.where(np.logical_and(is_near_pole_hole, is_nan))
                sic[indices] = opt.pole_hole_fill_value

        # All done, register the variable
        self.register_auxvar("sic", "sea_ice_concentration", sic, None)

    def load_requested_auxdata(self):
        """ Required subclass method: Load the data file necessary to satisfy condition for requested date"""

        # Retrieve the file path for the requested date from a property of the auxdata parent class
        path = self.requested_filepath

        #  --- Validation ---
        if not os.path.isfile(path):
            msg = self.pyclass+": File not found: %s " % path
            self.add_handler_message(msg)
            self.error.add_error("auxdata_missing_sic", msg)
            return

        # --- Read the data ---
        self._data = ReadNC(path)

        # --- Pre-process the data ---
        # Remove time dimension
        self._data.ice_conc = self._data.ice_conc[0, :, :]

        # No negative ice concentrations
        flagged = np.where(self._data.ice_conc < 0)
        self._data.ice_conc[flagged] = np.nan

        self.add_handler_message("OsiSafSIC: Loaded SIC file: %s" % path)

    def _get_sic_track(self, l2):
        """ Simple extraction along trajectory"""

        # Extract from grid
        griddef = self._options[l2.hemisphere]
        grid_lons, grid_lats = self._data.lon, self._data.lat
        grid2track = GridTrackInterpol(l2.track.longitude, l2.track.latitude, grid_lons, grid_lats, griddef)
        sic = grid2track.get_from_grid_variable(self._data.ice_conc, flipud=True)

        return sic

    @property
    def requested_filepath(self):
        """ Note: this overwrites the property in the super class due to some
        peculiarities with the filenaming (auto product changes etc) """

        # Unique to this class is the possiblity to auto merge
        # products. The current implementation supports only two products
        path = self._local_repository

        # The path needs to be completed if two products shall be used
        if "auto_product_change" in self.options:
            opt = self.options.auto_product_change
            product_index = int(self.start_time > opt.date_product_change)
            product_def = opt.osisaf_product_def[product_index]
            path = os.path.join(path, product_def["subfolder"])
            self.set_filenaming(product_def["filenaming"])
            self.set_long_name(product_def["long_name"])

        for subfolder_tag in self._subfolders:
            subfolder = getattr(self, subfolder_tag)
            path = os.path.join(path, subfolder)

        filename = self._filenaming.format(
            year=self.year, month=self.month, day=self.day,
            hemisphere_code=self.hemisphere_code)
        path = os.path.join(path, filename)
        return path


class IfremerSIC(AuxdataBaseClass):

    def __init__(self):
        super(IfremerSIC, self).__init__()
        self._data = None
        self._current_date = [0, 0, 0]
        self._requested_date = [-1, -1, -1]
        self.error.caller_id = self.__class__.__name__

    def subclass_init(self):
        """ Read the grid information """
        # XXX: This is a dirty hack, but needed for getting SIC lon/lat grid
        self._grid = {}
        for hemisphere in ["north", "south"]:
            grid_file = os.path.join(
                self._local_repository, "grid_%s_12km.nc" % hemisphere)
            self._grid[hemisphere] = ReadNC(grid_file)

    def get_l2_track_vars(self, l2):
        self._get_requested_date(l2)
        self._get_data(l2)
        sic = self._get_sic_track(l2)
        # All done, register the variable
        self.register_auxvar("sic", "sea_ice_concentration", sic, None)

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

        # Validation
        if not os.path.isfile(path):
            msg ="IfremerSIC: File not found: %s " % path
            self.add_handler_message(msg)
            self.error.add_error("auxdata_missing_sic", msg)
            return

        self._data = ReadNC(path)
        self._data.ice_conc = self._data.concentration[0, :, :]
        flagged = np.where(np.logical_or(self._data.ice_conc < 0, self._data.ice_conc > 100))
        self._data.ice_conc[flagged] = 0

        # This step is important for calculation of image coordinates
        self._data.ice_conc = np.flipud(self._data.ice_conc)
        self.add_handler_message("IfremerSIC: Loaded SIC file: %s" % path)
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
        grid = self._grid[l2.hemisphere]
        x, y = p(grid.longitude, grid.latitude)
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
