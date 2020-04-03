# -*- coding: utf-8 -*-
"""
Created on Fri May 19 18:16:09 2017

@author: Stefan
"""

from pysiral import get_cls
from pysiral.auxdata import AuxClassConfig
from pysiral.config import ConfigInfo
from pysiral.logging import DefaultLoggingClass
from pysiral.errorhandler import ErrorStatus, PYSIRAL_ERROR_CODES
from pysiral.iotools import get_local_l1bdata_files

from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta

import glob
from pathlib import Path
import re


class DefaultAuxdataClassHandler(DefaultLoggingClass):
    """ Class for retrieving handler classes for auxiliary data
    (mss, sic, sitype, snow). The classes are initialized with directory
    information from the local machine definition and the auxdata information
    from `auxdata.yaml` configuration file.
    """

    def __init__(self):
        super(DefaultAuxdataClassHandler, self).__init__(self.__class__.__name__)
        self.pysiral_config = ConfigInfo()
        self.error = ErrorStatus(caller_id=self.__class__.__name__)

    def get_pyclass(self, auxdata_class, auxdata_id, l2_procdef_opt):
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

        # Initialize the class with information from auxdata_def.yaml
        auxdata_def = self.get_auxdata_def(auxdata_class, auxdata_id)
        if auxdata_def is None:
            error_id = "auxdata_missing_definition"
            error_message = PYSIRAL_ERROR_CODES[error_id] % (auxdata_class, auxdata_id)
            self.error.add_error(error_id, error_message)
            self.error.raise_on_error()

        # Set the auxdata config
        cfg = AuxClassConfig()

        # connect to repository on local machine
        if "local_repository" in auxdata_def:
            local_repository_id = auxdata_def.local_repository
            local_repo = self.get_local_repository(auxdata_class, local_repository_id)
            if local_repo is None and local_repository_id is not None:
                error_id = "auxdata_missing_localrepo_def"
                error_message = PYSIRAL_ERROR_CODES[error_id] % (auxdata_class, auxdata_id)
                self.error.add_error(error_id, error_message)
                self.error.raise_on_error()
            cfg.set_local_repository(local_repo)

        # set doc str (should be mandatory for all auxdata handlers)
        if "long_name" in auxdata_def:
            cfg.set_long_name(auxdata_def.long_name)

        # set filename (e.g. for mss)
        if "filename" in auxdata_def:
            local_repository_id = auxdata_def.local_repository
            local_repo = self.get_local_repository(auxdata_class, local_repository_id)
            filename = Path(local_repo) / auxdata_def.filename
            cfg.set_filename(filename)

        # set filenaming (e.g. for sic, sitype, snow)
        if "filenaming" in auxdata_def:
            cfg.set_filenaming(auxdata_def.filenaming)

        # set subfolders (e.g. for sic, sitype, snow)
        if "subfolders" in auxdata_def:
            cfg.set_subfolder(auxdata_def.subfolders)

        # Set the default options from the auxiliary definition file
        if "options" in auxdata_def:
            options = auxdata_def.get("options", None)
            if options is not None:
                cfg.set_options(**options)

        # Override option with definition from the l2 processor settings
        if l2_procdef_opt is not None:
            cfg.set_options(**l2_procdef_opt)

        # Get the auxiliary data class
        module_name, class_name = "pysiral.auxdata.%s" % (auxdata_class), auxdata_def["pyclass"]
        auxclass = get_cls(module_name, class_name)
        if auxclass is None:
            error_id = "auxdata_invalid_class_name"
            msg = "Invalid Auxdata class: %s.%s" % (module_name, class_name)
            self.error.add_error(PYSIRAL_ERROR_CODES[error_id], msg)
            self.error.raise_on_error()

        # Init the auxiliary class
        # Note: This will trigger any action defined in the subclasses, such as reading static background files
        auxdata_handler = auxclass(cfg)

        # All done, return
        return auxdata_handler

    def get_local_repository(self, auxdata_class, auxdata_id):
        """ Get the local repository for the the auxdata type and id """
        if auxdata_id is None:
            return None
        aux_repo_defs = self.pysiral_config.local_machine.auxdata_repository
        try:
            local_repo_auxclass = aux_repo_defs[auxdata_class]
        except KeyError:
            msg = "Missing auxdata definition in local_machine_def.yaml: auxdata_repository.%s" % auxdata_class
            self.error.add_error("missing-localmachinedef-tag", msg)
            self.error.raise_on_error()
        return local_repo_auxclass.get(auxdata_id, None)

    def get_auxdata_def(self, auxdata_class, auxdata_id):
        """ Returns the definition in `config/auxdata_def.yaml` for
        specified auxdata class and id """
        try:
            auxdata_class_def = self.pysiral_config.auxdata[auxdata_class]
        except KeyError:
            msg = "Invalid auxdata class [%s] in auxdata_def.yaml" % auxdata_class
            self.error.add_error("invalid-auxdata-class", msg)
            self.error.raise_on_error()
        return auxdata_class_def.get(auxdata_id, None)


class DefaultL1bDataHandler(DefaultLoggingClass):
    """ Class for retrieving default l1b directories and filenames """

    def __init__(self, mission_id, hemisphere, version="default"):
        super(DefaultL1bDataHandler, self).__init__(self.__class__.__name__)
        self._mission_id = mission_id
        self._hemisphere = hemisphere
        self._version = version
        self._last_directory = None

    def get_files_from_time_range(self, time_range):
        files, search_directory = get_local_l1bdata_files(
                self._mission_id, time_range, self._hemisphere,
                version=self._version)
        self._last_directory = search_directory
        return files

    @property
    def last_directory(self):
        return str(self._last_directory)


class L2iDataHandler(DefaultLoggingClass):
    """ Class for retrieving default l1b directories and filenames """

    def __init__(self, base_directory, force_l2i_subfolder=True):
        super(L2iDataHandler, self).__init__(self.__class__.__name__)
        self.error = ErrorStatus(caller_id=self.__class__.__name__)
        self._base_directory = base_directory
        self._force_l2i_subfolder = force_l2i_subfolder
        self._subdirectory_list = self.get_subdirectory_list()
        self._validate_base_directory()

    def get_files_from_time_range(self, time_range):
        """ Get all files that fall into time range (May be spread over
        the different year/ month subfolders """
        l2i_files = []
        for year, month, day in time_range.days_list:
            lookup_directory = self.get_lookup_directory(year, month)
            if not Path(lookup_directory).is_dir():
                continue
            l2i_pattern = self.get_l2i_search_str(year=year, month=month, day=day)
            result = Path(lookup_directory).glob(l2i_pattern)
            l2i_files.extend(sorted(result))
        return l2i_files

    def get_files_for_day(self, day_dt):
        """ Retrieve a list of l2i files with data points for a given day.
        Also specifically looks for files with had a start time on the
        previous day """

        # Get the lookup directory
        lookup_directory = self.get_lookup_directory(day_dt.year, day_dt.month)

        # XXX: We are not evaluating the netCDF attributes at this point
        #      but assuming that the filename contains start and stop
        #      time. This is a pretty safe assumption, but this approach
        #      should be replaced as soon as a proper inspection tool is
        #      available
        day_search = self.get_l2i_search_str(year=day_dt.year, month=day_dt.month, day=day_dt.day)
        l2i_files = Path(lookup_directory).glob(day_search)

        # Check if day is the first day of the month
        # yes -> check last file of previous month which might have data
        #        for the target day
        if day_dt.day == 1:
            previous_day = day_dt - timedelta(days=1)
            lookup_directory = self.get_lookup_directory(previous_day.year, previous_day.month)
            additional_l2i_files = Path(lookup_directory).glob(day_search)
            l2i_files.extend(additional_l2i_files)

        # All done, return sorted output
        return sorted(l2i_files)

    def _validate_base_directory(self):
        """ Performs sanity checks and enforces the l2i subfolder """
        # 1. Path must exist
        if not Path(self._base_directory).is_dir():
            msg = "Invalid l2i product directory: %s"
            msg = msg % str(self._base_directory)
            self.error.add_error("invalid-l2i-productdir", msg)
            self.error.raise_on_error()

    def get_lookup_directory(self, year, month):
        """ Return the sub folders for a given time (datetime object) """
        subfolders = ["%4g" % year, "%02g" % month]
        lookup_directory = Path(self.product_basedir) / "/".join(subfolders)
        return lookup_directory

    def get_subdirectory_list(self):
        """ Returns a list of all subdirectories of type yyyy/mm """
        subdirectory_list = list()
        try:
            years = sorted([f for f in Path(self.product_basedir.iterdir()) if f.is_dir()])
        except StopIteration:
            self.log.warning("No subdirectories in %s" % self.product_basedir)
            return []
        # filter any invalid directories
        years = [y for y in years if re.match(r'[1-3][0-9]{3}', y)]
        for year in years:
            subdir_year = Path(self.product_basedir) / year
            months = [f for f in subdir_year.iterdir() if f.is_dir()]
            # filter any invalid directories
            months = [m for m in months if re.match(r'[0-1][0-9]', m)]
            subdirectory_list.extend([[year, m] for m in months])
        return subdirectory_list

    def get_l2i_search_str(self, year=None, month=None, day=None):
        """ Returns a search pattern for l2i files with optional refined
        search for year, month, day. Note: month & day can only be set,
        if the year & year + month respectively is set
        Examples:
            *l2i*.nc
            *l2i*2017*.nc
            *l2i*201704*.nc
            *l2i*20170401*.nc
        """
        date_str = "*"
        if year is not None:
            date_str += "%04g" % year
        if month is not None and year is not None:
            date_str += "%02g" % month
        else:
            raise ValueError("year must be set if month is set")
        if day is not None and month is not None:
            date_str += "%02g" % day
        else:
            raise ValueError("year & month must be set if day is set")
        if len(date_str) > 1:
            date_str += "*"
        l2i_file_pattern = "*l2i%s.nc" % date_str
        return l2i_file_pattern

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
        return datetime(int(last_month[0]), int(last_month[1]), 1) + relativedelta(months=1, microseconds=-1)
