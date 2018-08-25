# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 13:57:56 2016

@author: Stefan
"""

from pysiral.auxdata import AuxdataBaseClass, GridTrackInterpol
from pysiral.iotools import ReadNC

import scipy.ndimage as ndimage
from pyproj import Proj
import numpy as np
import os


class SICBaseClass(AuxdataBaseClass):

    def __init__(self):
        super(SICBaseClass, self).__init__()
        self._msg = ""

    def get_along_track_sic(self, l2):
        sic, msg = self._get_along_track_sic(l2)
        return sic, msg


class OsiSafSIC(SICBaseClass):

    def __init__(self):
        super(OsiSafSIC, self).__init__()
        self._data = None
        self.error.caller_id = self.__class__.__name__

    def _initialize(self):
        pass

    def _get_along_track_sic(self, l2):
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
            return None, self._msg

        # Get and return the track
        sic = self._get_sic_track(l2)


        # All done, return
        return sic, self._msg

    def load_requested_auxdata(self):
        """ Required subclass method: Load the data file necessary to satisfy condition for requested date"""

        # Retrieve the file path for the requested date from a property of the auxdata parent class
        path = self.requested_filepath

        #  --- Validation ---
        if not os.path.isfile(path):
            self._msg = self.__class__.__name__+": File not found: %s " % path
            self.error.add_error("auxdata_missing_sic", self._msg)
            return

        # --- Read the data ---
        self._data = ReadNC(path)

        # --- Pre-process the data ---
        # Remove time dimension
        self._data.ice_conc = self._data.ice_conc[0, :, :]

        # No negative ice concentrations
        flagged = np.where(self._data.ice_conc < 0)
        self._data.ice_conc[flagged] = 0

        self._msg = "OsiSafSIC: Loaded SIC file: %s" % path

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

class IfremerSIC(SICBaseClass):

    def __init__(self):
        super(IfremerSIC, self).__init__()
        self._data = None
        self._current_date = [0, 0, 0]
        self._requested_date = [-1, -1, -1]
        self.error.caller_id = self.__class__.__name__

    @property
    def year(self):
        return "%04g" % self._requested_date[0]

    @property
    def month(self):
        return "%02g" % self._requested_date[1]

    @property
    def day(self):
        return "%02g" % self._requested_date[2]

    def _initialize(self):
        """ Read the grid information """
        # XXX: This is a dirty hack, but needed for getting SIC lon/lat grid
        self._grid = {}
        for hemisphere in ["north", "south"]:
            grid_file = os.path.join(
                self._local_repository, "grid_%s_12km.nc" % hemisphere)
            self._grid[hemisphere] = ReadNC(grid_file)

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

        # Validation
        if not os.path.isfile(path):
            self._msg = "IfremerSIC: File not found: %s " % path
            self.error.add_error("auxdata_missing_sic", self._msg)
            return

        self._data = ReadNC(path)
        self._data.ice_conc = self._data.concentration[0, :, :]
        flagged = np.where(
            np.logical_or(self._data.ice_conc < 0, self._data.ice_conc > 100))
        self._data.ice_conc[flagged] = 0

        # This step is important for calculation of image coordinates
        self._data.ice_conc = np.flipud(self._data.ice_conc)
        self._msg = "IfremerSIC: Loaded SIC file: %s" % path
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

#        # XXX: Debug stuff below
#
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
#
#        plt.imshow(self._data.ice_conc, cmap=plt.get_cmap("viridis"),
#                   origin="upper", interpolation="none")
#        plt.plot(ix, iy)
#
#        plt.figure("sic")
#        plt.plot(sic)
#        plt.show()
#        stop


def get_l2_sic_handler(name):
    pyclass = globals().get(name, None)
    if pyclass is not None:
        return pyclass()
    else:
        return pyclass
