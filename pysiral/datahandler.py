# -*- coding: utf-8 -*-
"""
Created on Fri May 19 18:16:09 2017

@author: Stefan
"""

from pysiral.config import ConfigInfo
from pysiral.logging import DefaultLoggingClass
from pysiral.errorhandler import ErrorStatus, PYSIRAL_ERROR_CODES
from pysiral.sic import get_l2_sic_handler
from pysiral.sitype import get_l2_sitype_handler
from pysiral.snow import get_l2_snow_handler
from pysiral.mss import get_l2_ssh_class
import os


class DefaultAuxdataHandler(DefaultLoggingClass):
    """ Class for retrieving handler classes for auxiliary data
    (mss, sic, sitype, snow). The classes are initialized with directory
    information from the local machine definition and the auxdata information
    from `auxdata.yaml` configuration file.
    """

    valid_auxdata_classes = ["mss", "sic", "sitype", "snow"]
    getcls = {"mss": get_l2_ssh_class, "sic": get_l2_sic_handler,
              "sitype": get_l2_sitype_handler, "snow": get_l2_snow_handler}

    def __init__(self):
        super(DefaultAuxdataHandler, self).__init__(self.__class__.__name__)
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
            if local_repo is None and local_repository_id is not None:
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
            options = auxdata_def.get("options", None)
            if options is not None:
                auxdata_handler.set_options(**options)

        return auxdata_handler

    def get_local_repository(self, auxdata_class, auxdata_id):
        """ Get the local repository for the the auxdata type and id """
        if auxdata_id is None:
            return None
        local_repo = self.pysiral_config.local_machine.auxdata_repository
        return local_repo[auxdata_class].get(auxdata_id, None)

    def get_auxdata_def(self, auxdata_class, auxdata_id):
        """ Returns the definition in `config/auxdata_def.yaml` for
        specified auxdata class and id """
        auxdata_class_def = self.pysiral_config.auxdata[auxdata_class]
        return auxdata_class_def.get(auxdata_id, None)
