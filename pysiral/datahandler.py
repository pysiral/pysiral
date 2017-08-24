# -*- coding: utf-8 -*-
"""
Created on Fri May 19 18:16:09 2017

@author: Stefan
"""

from pysiral.config import ConfigInfo
from pysiral.logging import DefaultLoggingClass
from pysiral.errorhandler import ErrorStatus, PYSIRAL_ERROR_CODES
from pysiral.iotools import get_local_l1bdata_files
from pysiral.sic import get_l2_sic_handler
from pysiral.sitype import get_l2_sitype_handler
from pysiral.snow import get_l2_snow_handler
from pysiral.mss import get_l2_ssh_class

from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import glob
import os
import re


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

        # Allow l2 settings to not specify a data handler
        # In this case a dummy data handler is returned
        if auxdata_id is None:
            return self.getcls[auxdata_class]("NoneHandler")

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


class DefaultL1bDataHandler(DefaultLoggingClass):
    """ Class for retrieving default l1b directories and filenames """

    def __init__(self, mission_id, hemisphere, version="default"):
        super(DefaultL1bDataHandler, self).__init__(self.__class__.__name__)
        self._mission_id = mission_id
        self._hemisphere = hemisphere
        self._version = version

    def get_files_from_time_range(self, time_range):
        return get_local_l1bdata_files(
                self._mission_id, time_range, self._hemisphere,
                version=self._version)


class L2iDataHandler(DefaultLoggingClass):
    """ Class for retrieving default l1b directories and filenames """

    l2i_file_pattern = "l2i*.nc"

    def __init__(self, base_directory, force_l2i_subfolder=True):
        super(L2iDataHandler, self).__init__(self.__class__.__name__)
        self.error = ErrorStatus(caller_id=self.__class__.__name__)
        self._base_directory = base_directory
        self._force_l2i_subfolder = force_l2i_subfolder
        self._subdirectory_list = self.get_subdirectory_list()
        self._validate_base_directory()

    def get_files_from_time_range(self, time_range):
        """ XXX: This currently hard-coded to monthly periods """
        lookup_directory = self.get_lookup_directory(time_range.start)
        if not os.path.isdir(lookup_directory):
            return []
        l2i_files = glob.glob(os.path.join(lookup_directory,
                                           self.l2i_file_pattern))
        return sorted(l2i_files)

    def get_files_for_day(self, day_dt):
        """ Retrieve a list of l2i files with data points for a given day.
        Also specifically looks for files with had a start time on the
        previous day """

        # Get the lookup directory
        lookup_directory = self.get_lookup_directory(day_dt)

        # XXX: We are not evaluating the netCDF attributes at this point
        #      but assuming that the filename contains start and stop
        #      time. This is a pretty safe assumption, but this approach
        #      should be replaced as soon as a proper inspection tool is
        #      available
        day_search = "*%s*" % day_dt.strftime("%Y%m%d")
        day_search = self.l2i_file_pattern.replace("*", day_search)
        search_str = os.path.join(lookup_directory, day_search)
        l2i_files = glob.glob(search_str)

        # Check if day is the first day of the month
        # yes -> check last file of previous month which might have data
        #        for the target day
        if day_dt.day == 1:
            previous_day = day_dt - timedelta(days=1)
            lookup_directory = self.get_lookup_directory(previous_day)
            search_str = os.path.join(lookup_directory, day_search)
            additional_l2i_files = glob.glob(search_str)
            l2i_files.extend(additional_l2i_files)

        # All done, return sorted output
        return sorted(l2i_files)

    def _validate_base_directory(self):
        """ Performs sanity checks and enforces the l2i subfolder """
        # 1. Path must exist
        if not os.path.isdir(self._base_directory):
            msg = "Invalid l2i product directory: %s"
            msg = msg % str(self._base_directory)
            self.error.add_error("invalid-l2i-productdir", msg)
            self.error.raise_on_error()

    def get_lookup_directory(self, dt):
        """ Return the sub folders for a given time (datetime object) """
        subfolders = ["%4g" % dt.year, "%02g" % dt.month]
        lookup_directory = os.path.join(self.product_basedir, *subfolders)
        return lookup_directory

    def get_subdirectory_list(self):
        """ Returns a list of all subdirectories of type yyyy/mm """
        subdirectory_list = list()
        years = sorted(next(os.walk(self.product_basedir))[1])
        # filter any invalid directories
        years = [y for y in years if re.match(r'[1-3][0-9]{3}', y)]
        for year in years:
            subdir_year = os.path.join(self.product_basedir, year)
            months = sorted(next(os.walk(subdir_year))[1])
            # filter any invalid directories
            months = [m for m in months if re.match(r'[0-1][0-9]', m)]
            subdirectory_list.extend([[year, m] for m in months])
        return subdirectory_list

    @property
    def product_basedir(self):
        return self._base_directory

    @property
    def subdirectory_list(self):
        return self._subdirectory_list

    @property
    def start_month(self):
        """ Returns a date time object for the first month of the l2i
        product repository """
        first_month = self.subdirectory_list[0]
        return datetime(int(first_month[0]), int(first_month[1]), 1)

    @property
    def stop_month(self):
        """ Returns a date time object for the last month of the l2i
        product repository """
        last_month = self.subdirectory_list[-1]
        return datetime(int(last_month[0]), int(last_month[1]), 1) + \
            relativedelta(months=1, microseconds=-1)
