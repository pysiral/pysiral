# -*- coding: utf-8 -*-
"""
Created on Fri May 19 14:42:35 2017

@author: Stefan
"""

__all__ = ["mss", "icechart", "rio", "sic", "sitype", "snow", "region"]


import os
import re

import numpy as np

import scipy.ndimage as ndimage


from pyproj import Proj

from pysiral.config import options_from_dictionary
from pysiral.errorhandler import ErrorStatus


class AuxdataBaseClass(object):
    """
    Base class for all sub-type auxdata base classes (e.g. SICBaseClass).
    This class defines the mandatory set of methods and properties for all
    auxdata classes
    """

    def __init__(self, auxclass_cfg):
        """ This class should not be called directly, only its subclasses. auxclass_cfg needs to be of type
        AuxClassConfig """

        # Error handler
        self.error = ErrorStatus(self.pyclass)

        # Auxiliary class options
        if not isinstance(auxclass_cfg, AuxClassConfig):
            msg = "Invalid config object: %s (needs to be of type pysiral.auxdata.AuxClassConfig"
            msg = msg % str(auxclass_cfg)
            self.error.add_error("invalid-auxclasscfg-type", msg)
            self.error.raise_on_error()
        self._cfg = auxclass_cfg

        # General messages
        self.msgs = []

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
        self._auxvars = []

    def reset_handler_messages(self):
        """ Empties the message list. To be executed during class initialization and
        before retrieving data (e.g. since the Level-2 processor calls this instance repeatedly) """
        self.msgs = []

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
        self.reset_handler_messages()

        # Call the mandatory track extraction method. Each subclass should register its output via the
        # `register_auxvar` method of the parent class
        self.get_l2_track_vars(l2)

        # Check on errors
        if self.error.status and self.exception_on_error:
            self.error.raise_on_error()

        # Update the Level-2 object
        try:
            self.update_l2(l2)
        except KeyError:
            msg = "Invalid auxiliary parameter return from class %s" % self.pyclass
            self.error.add_error("invalid-auxvar-return", msg)
            self.error.raise_on_error()

    def register_auxvar(self, var_id, var_name, value, uncertainty=None):
        """ Register an auxiliary variable. The different parameters are necessary for the L2 data object.
        When it will be added to the l2 object in self.update_l2, the variable will be accessible from the l2 with
        the following expressions:

            value = l2.%var_id%
            uncertainty = l2.%var_id%.uncertainty

        or

            value = l2.get_parameter_by_name(%var_name%)
            uncertainty = l2.get_parameter_by_name(%var_name%_uncertainty)
        """
        auxvar_dict = dict(id=var_id, name=var_name, value=value, uncertainty=uncertainty)
        self._auxvars.append(auxvar_dict)

    def add_handler_message(self, msg):
        self.msgs.append(msg)

    def get_empty_array(self, l2, empty_val=np.nan):
        return np.full((l2.n_records), empty_val)

    def update_external_data(self):
        """ This method will check if the requested date matches current data
        and call the subclass data loader method if not """
        # Check if data for day is already loaded
        if self._requested_date != self._current_date:
            # NOTE: The implementation of this method needs to be in the subclass
            self.load_requested_auxdata()
            self._current_date = self._requested_date
            if self.has_data_loaded:
                self.add_handler_message(self.__class__.__name__ + ": Load "+self.requested_filepath)
        else:
            if self.has_data_loaded:
                self.add_handler_message(self.__class__.__name__+": Data already present")
            else:
                msg = ": No Data: Loading failed in an earlier attempt"
                self.add_handler_message(self.__class__.__name__ + msg)

    def update_l2(self, l2):
        """ Automatically add all auxiliary variables to a Level-2 data object"""
        for auxvar in self._auxvars:
            uncertainty = auxvar.get("uncertainty", None)
            l2.set_auxiliary_parameter(auxvar["id"], auxvar["name"], auxvar["value"], uncertainty)

    @property
    def pyclass(self):
        return self.__class__.__name__

    @property
    def cfg(self):
        return self._cfg

    @property
    def has_data_loaded(self):
        if not hasattr(self, "_data"):
            return False
        return self._data is not None

    @property
    def exception_on_error(self):
        if self.cfg.options.has_key("exception_on_error"):
            exception_on_error = self.cfg.options.exception_on_error
        else:
            exception_on_error = False
        return exception_on_error

    @property
    def requested_filepath(self):
        """ Returns the local file path for the requested date"""

        # Main directory
        path = self.cfg.local_repository

        # Add the subfolders
        for subfolder_tag in self.cfg.subfolders:
            subfolder = getattr(self, subfolder_tag)
            path = os.path.join(path, subfolder)

        # Get the period dict (will be constructed from filenaming)
        period_dict = {}
        attrs = re.findall("{.*?}", self.cfg.filenaming)
        for attr_def in attrs:
            attr_name = attr_def[1:-1]
            period_dict[attr_name] = getattr(self, attr_name)
        filename = self.cfg.filenaming.format(**period_dict)
        path = os.path.join(path, filename)
        return path

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


class AuxClassConfig(object):
    """ A container for configuration data for any auxilary data handler class"""

    def __init__(self):

        # General options
        self.options = None
        self.local_repository = None
        self.filename = None
        self.filenaming = None
        self.subfolders = []
        self.long_name = ""

    def set_long_name(self, docstr):
        """ Set a description of the auxdata source """
        self.long_name = docstr

    def set_options(self, **opt_dict):
        """  Pass a dictionary with options """
        if self.options is None:
            self.options = options_from_dictionary(**opt_dict)
        else:
            self.options.update(options_from_dictionary(**opt_dict))

    def set_local_repository(self, path):
        """ Set the path the local auxdata repository """
        self.local_repository = path

    def set_filename(self, filename):
        """ Set a constant filename (e.g. for mss) """
        self.filename = filename

    def set_filenaming(self, filenaming):
        """ Set the filenaming of the auxdata files """
        self.filenaming = filenaming

    def set_subfolder(self, subfolder_list):
        """ Set a list of folders (year, [month, [day]]) where auxdata files
        can be found """
        self.subfolders = subfolder_list


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



