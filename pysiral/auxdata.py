# -*- coding: utf-8 -*-
"""
Created on Fri May 19 14:42:35 2017

@author: Stefan
"""

from pysiral.config import options_from_dictionary
from pysiral.errorhandler import ErrorStatus


class AuxdataBaseClass(object):
    """ Base class for all sub-type auxdata base classes (e.g. SICBaseClass).
    This class defines the mandatory set of methods and properties for
    all auxdata classes
    """

    def __init__(self):
        self._options = None
        self._local_repository = None
        self._filename = None
        self._filenaming = None
        self._subfolders = []
        self._long_name = ""
        self.error = ErrorStatus()

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

    def initialize(self, *args, **kwargs):
        """ Initialize before Level-2 processing """
        self._initialize(*args, **kwargs)

    @property
    def pyclass(self):
        return self.__class__.__name__

    @property
    def filename(self):
        return self._filename

    @property
    def longname(self):
        return self._long_name

    @property
    def options(self):
        return self._options.makeReport()
