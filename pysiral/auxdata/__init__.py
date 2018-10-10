# -*- coding: utf-8 -*-
"""
Created on Fri May 19 14:42:35 2017

@author: Stefan
"""

import os
import re

import numpy as np

import scipy.ndimage as ndimage

from collections import OrderedDict

from pyproj import Proj

from pysiral.config import options_from_dictionary
from pysiral.errorhandler import ErrorStatus


class AuxdataBaseClass(object):
    """
    Base class for all sub-type auxdata base classes (e.g. SICBaseClass).
    This class defines the mandatory set of methods and properties for all
    auxdata classes
    """

    def __init__(self):

        # Error handler
        self.error = ErrorStatus()

        # General messages
        self.msgs = []

        # --- This will be filled by the set_xxx methods ---

        # General options
        self._options = None
        self._local_repository = None
        self._filename = None
        self._filenaming = None
        self._subfolders = []
        self._long_name = ""

        # --- Class internals ---

        # This is for auxiliary data handlers that require to read external product files for
        # a defined period (daily, monthly, ...). The implementation currently keeps only one
        # external product in memory at the time. The period (date list: yyyy, mm, dd) of this
        # currently loaded product is designated as current_date  This date is compared to the
        # requested date and if a new product is loaded upon mismatch of current & requested data
        # NOTE: This will be bypassed by static auxiliary data classes
        # TODO: Load all auxiliary products for processing period in memory (allow parallel processing)
        self._current_date = [0, 0, 0]
        self._requested_date = [-1, -1, -1]

        # A dictionary with the output variables of the auxiliary data set
        self.reset_auxvars()

    def set_long_name(self, docstr):
        """ Set a description of the auxdata source """
        self._long_name = docstr

    def set_options(self, **opt_dict):
        """  Pass a dictionary with options """
        if self._options is None:
            self._options = options_from_dictionary(**opt_dict)
        else:
            self._options.update(options_from_dictionary(**opt_dict))

    def set_local_repository(self, path):
        """ Set the path the local auxdata repository """
        self._local_repository = path

    def set_filename(self, filename):
        """ Set a constant filename (e.g. for mss) """
        self._filename = filename

    def set_filenaming(self, filenaming):
        """ Set the filenaming of the auxdata files """
        self._filenaming = filenaming

    def set_subfolder(self, subfolder_list):
        """ Set a list of folders (year, [month, [day]]) where auxdata files
        can be found """
        self._subfolders = subfolder_list

    def set_requested_date(self, year, month, day):
        """ Use first timestamp as reference, date changes are ignored """
        self._requested_date = [year, month, day]

    def set_requested_date_from_l2(self, l2):
        """ Convenience method, Use first timestamp as reference, date changes are ignored """
        year = l2.track.timestamp[0].year
        month = l2.track.timestamp[0].month
        day = l2.track.timestamp[0].day
        self.set_requested_date(year, month, day)

    def reset_auxvars(self):
        """ Empties the auxiliary data store. To be executed during class initialization and
        before retrieving data (e.g. since the Level-2 processor calls this instance repeatedly) """
        self._auxvars = OrderedDict()

    def add_variables_to_l2(self, l2):
        """ Main Access points for the Level-2 Processor """

        # Call the API get_track class. This is the mandatory method of all auxiliary subclasses (independent
        # of type. Test if this is indeed the case
        if not self.has_mandatory_track_method:
            msg = "Mandatory subclass method `get_l2_track_vars` not implemented for %s " % self.pyclass
            self.error.add_error("not-implemented", msg)
            self.error.raise_on_error()

        # Before calling the get_track_vars of the subclass, we must empty any existing data from a potential
        # previous execution
        self.reset_auxvars()

        # Call the mandatory track extraction method. Each subclass should register its output via the
        # `register_auxvar` method of the parent class
        self.get_l2_track_vars(l2)

        # Update the Level-2 object
        self.update_l2(l2)

    def register_auxvar(self, var_name, var):
        self._auxvars[var_name] = var

    def update_external_data(self):
        """ This method will check if the requested date matches current data
        and call the subclass data loader method if not """
        # Check if data for day is already loaded
        if self._requested_date != self._current_date:
            # NOTE: The implementation of this method needs to be in the subclass
            self.load_requested_auxdata()
            self._current_date = self._requested_date
            self._msg = self.__class__.__name__ + ": Load "+self.requested_filepath
        else:
            self._msg = self.__class__.__name__+": Data already present"

    def initialize(self, *args, **kwargs):
        """ Initialize before Level-2 processing """
        # This executes the _initialize method of each children class
        self._initialize(*args, **kwargs)

    def update_l2(self, l2):
        """ Automatically add all auxiliary variables to a Level-2 data object"""
        # TODO: this has yet no functionality
        pass

    @property
    def pyclass(self):
        return self.__class__.__name__

    @property
    def filename(self):
        return self._filename

    @property
    def requested_filepath(self):
        """ Returns the local file path for the requested date"""

        # Main directory
        path = self._local_repository

        # Add the subfolders
        for subfolder_tag in self._subfolders:
            subfolder = getattr(self, subfolder_tag)
            path = os.path.join(path, subfolder)

        # Get the period dict (will be constructed from filenaming)
        period_dict = {}
        attrs = re.findall("{.*?}", self._filenaming)
        for attr_def in attrs:
            attr_name = attr_def[1:-1]
            period_dict[attr_name] = getattr(self, attr_name)
        filename = self._filenaming.format(**period_dict)
        path = os.path.join(path, filename)
        return path

    @property
    def longname(self):
        return self._long_name

    @property
    def options(self):
        return self._options

    @property
    def year(self):
        return "%04g" % self._requested_date[0]

    @property
    def month(self):
        return "%02g" % self._requested_date[1]

    @property
    def day(self):
        return "%02g" % self._requested_date[2]

    @property
    def auxdata_variable_names(self):
        """ Returns a list of variables that are provided by the auxiliary data class """
        return sorted(self._variables.keys())

    @property
    def has_mandatory_track_method(self):
        """ Test if this object instance has the mandatory method for extracting track data. This method
        is named get_l2_track_vars() and needs to be present in any auxiliary subclass"""
        has_method = False
        get_track_children_method = getattr(self, "get_l2_track_vars", None)
        if callable(get_track_children_method):
            has_method = True
        return has_method

    @property
    def auxvar_names(self):
        return self._auxvars.keys()

class GridTrackInterpol(object):
    """ Implements fast extraction of gridded data along a track using Image Interpolation """

    def __init__(self, lons, lats, grid_lons, grid_lats, griddef):
        """
        lons, lats: ground track
        Example grid definition (dict)
            projection:
                proj: stere
                ellps: WGS84
                lon_0: 0
                lat_0: -90
                lat_ts: -70
                a: 6378273
                b: 6356889.44891
            dimension:
                n_cols: 632
                n_lines: 664
                dx: 12500
                dy: 12500
        """

        # Save the arguments
        self.lons = lons
        self.lats = lats
        self.grid_lons = grid_lons
        self.grid_lats = grid_lats
        self.griddef = griddef

        # Compute image coordinates
        self._set_projection()
        self._get_track_image_coordinates()

    def _set_projection(self):
        self.p = Proj(**self.griddef.projection)

    def _get_track_image_coordinates(self):
        """ Computes the image coordinates that will be used for the m"""

        # Convert track coordinates to grid projection coordinates
        tr_x, tr_y = self.p(self.lons, self.lats)

        # Convert grid coordinates to grid projection coordinates
        x, y = self.p(self.grid_lons, self.grid_lats)

        # Convert track projection coordinates to image coordinates
        # x: 0 < n_lines; y: 0 < n_cols
        dim = self.griddef.dimension
        x_min, y_min = np.nanmin(x), np.nanmin(y)
        self.ix, self.iy = (tr_x-x_min)/dim.dx, (tr_y-y_min)/dim.dy

    def get_from_grid_variable(self, gridvar, order=0, flipud=False):
        """ Returns a along-track data from a grid variable"""
        if flipud:
            gridvar = np.flipud(gridvar)
        trackvar = ndimage.map_coordinates(gridvar, [self.iy, self.ix], order=order)
        return trackvar

    def debug_map(self, *args, **kwargs):
        raise NotImplementedError()
        # track_var = self.get_from_grid_variable(*args, **kwargs)





