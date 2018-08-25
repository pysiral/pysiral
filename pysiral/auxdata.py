# -*- coding: utf-8 -*-
"""
Created on Fri May 19 14:42:35 2017

@author: Stefan
"""

import os
import re

import numpy as np

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
        self.msg = ""

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
        self._auxvars = {}

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

    def update_external_data(self):
        """ This method will check if the requested date matches current data
        and call the subclass data loader method if not """
        # Check if data for day is already loaded
        if self._requested_date != self._current_date:
            # NOTE: The implementation of this method needs to be in the subclass
            self.load_requested_auxdata()
            self._current_date = self._requested_date
            self._msg = self.__class__.__name__ + ": Load "+self.requested_filename
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
        attr_defs = re.findall("{.*?}", self._filenaming)

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
