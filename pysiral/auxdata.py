# -*- coding: utf-8 -*-
"""
Created on Fri May 19 14:42:35 2017

@author: Stefan
"""

from pysiral.config import ConfigInfo, options_from_dictionary
from pysiral.errorhandler import ErrorStatus, PYSIRAL_ERROR_CODES
from pysiral.sic import get_l2_sic_handler
from pysiral.sitype import get_l2_sitype_handler
from pysiral.snow import get_l2_snow_handler
from pysiral.mss import get_l2_ssh_class
import os


class DefaultAuxdataHandler(object):
    """ Class for retrieving handler classes for auxiliary data
    (mss, sic, sitype, snow). The classes are initialized with directory
    information from the local machine definition and the auxdata information
    from `auxdata.yaml` configuration file.
    """

    valid_auxdata_classes = ["mss", "sic", "sitype", "snow"]
    getcls = {"mss": get_l2_ssh_class, "sic": get_l2_sic_handler,
              "sitype": get_l2_sitype_handler, "snow": get_l2_snow_handler}

    def __init__(self):
        self.pysiral_config = ConfigInfo()
        self.error = ErrorStatus(caller_id=self.__class__.__name__)

    def get_pyclass(self, auxdata_class, auxdata_id):
        """
        Returns a class for handling auxiliary data files, that is initialized
        with auxdata settings in `config/auxdata_def.yaml` and with the
        directory specified in `local_machine_def.yaml`

        Args:
            auxdata_class (str): Auxdata class (e.g. mss, sic, sitype, snow)
            auxdata_id (str): Auxdata class identifier (e.g. osisaf)

        Returns:
            class: The initialized auxdata handler class
        """

        # Clear errors
        self.error.reset()

        # perform sanity check
        if auxdata_class not in self.valid_auxdata_classes:
            error_id = "auxdata_invalid_class"
            error_message = PYSIRAL_ERROR_CODES[error_id] % auxdata_class
            self.error.add_error(error_id, error_message)
            return None

        # Initialize the class with information from auxdata_def.yaml
        auxdata_def = self.get_auxdata_def(auxdata_class, auxdata_id)
        if auxdata_def is None:
            error_id = "auxdata_missing_definition"
            error_message = PYSIRAL_ERROR_CODES[error_id] % (
                    auxdata_class, auxdata_id)
            self.error.add_error(error_id, error_message)
            return None

        # retrieve the getter class
        auxdata_handler = self.getcls[auxdata_class](auxdata_def.pyclass)
        if auxdata_handler is None:
            error_id = "auxdata_invalid_class_name"
            self.error.add_error(PYSIRAL_ERROR_CODES[error_id])
            return None

        # connect to repository on local machine
        if "local_repository" in auxdata_def:
            local_repository_id = auxdata_def.local_repository
            local_repo = self.get_local_repository(
                    auxdata_class, local_repository_id)
            if local_repo is None:
                error_id = "auxdata_missing_localrepo_def"
                error_message = PYSIRAL_ERROR_CODES[error_id] % (
                    auxdata_class, auxdata_id)
                self.error.add_error(error_id, error_message)
                return None
            auxdata_handler.set_local_repository(local_repo)

        # set doc str (should be mandatory for all auxdata handlers)
        if "long_name" in auxdata_def:
            auxdata_handler.set_long_name(auxdata_def.long_name)

        # set filename (e.g. for mss)
        if "file" in auxdata_def:
            local_repository_id = auxdata_def.local_repository
            local_repo = self.get_local_repository(
                    auxdata_class, local_repository_id)
            filename = os.path.join(local_repo, auxdata_def.file)
            auxdata_handler.set_filename(filename)

        # set filenaming (e.g. for sic, sitype, snow)
        if "filenaming" in auxdata_def:
            auxdata_handler.set_filenaming(auxdata_def.filenaming)

        # set subfolders (e.g. for sic, sitype, snow)
        if "subfolders" in auxdata_def:
            auxdata_handler.set_subfolder(auxdata_def.subfolders)

        if "options" in auxdata_def:
            auxdata_handler.set_options(**auxdata_def.get("options", {}))

        return auxdata_handler

    def get_local_repository(self, auxdata_class, auxdata_id):
        """ Get the local repository for the the auxdata type and id """
        local_repo = self.pysiral_config.local_machine[auxdata_class]
        return local_repo.get(auxdata_id, None)

    def get_auxdata_def(self, auxdata_class, auxdata_id):
        """ Returns the definition in `config/auxdata_def.yaml` for
        specified auxdata class and id """
        auxdata_class_def = self.pysiral_config.auxdata[auxdata_class]
        return auxdata_class_def.get(auxdata_id, None)


def AuxdataBaseClass(object):
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
        self.long_name = ""
        self.error = ErrorStatus()

    def set_long_name(self, docstr):
        """ Set a description of the auxdata source """
        self.long_name = ""

    def set_options(self, **opt_dict):
        """  Pass a dictionary with options """
        self._options = options_from_dictionary(**opt_dict)

    def set_local_repository(self, path):
        """ Set the path the local auxdata repository """
        self._local_repository = path

    def set_filenaming(self, filenaming):
        """ Set the filenaming of the auxdata files """
        self._filenaming = filenaming

    def set_subfolders(self, subfolder_list):
        """ Set a list of folders (year, [month, [day]]) where auxdata files
        can be found """
        self._subfolders = subfolder_list

    def initialize(self, *args, **kwargs):
        """ Initialize before Level-2 processing """
        self._initialize(*args, **kwargs)

    @property
    def pyclass(self):
        return self.__class__.__name___

    @property
    def filename(self):
        return self._filename
