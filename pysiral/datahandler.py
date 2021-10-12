# -*- coding: utf-8 -*-
"""
Created on Fri May 19 18:16:09 2017

@author: Stefan
"""

import re
from typing import List, Union
from pathlib import Path
from loguru import logger
from attrdict import AttrDict

from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta

from dateperiods import DatePeriod

from pysiral import get_cls, psrlcfg
from pysiral.auxdata import AuxClassConfig
from pysiral._class_template import DefaultLoggingClass
from pysiral.errorhandler import ErrorStatus, PYSIRAL_ERROR_CODES
from pysiral.output import PysiralOutputFilenaming


class DefaultAuxdataClassHandler(DefaultLoggingClass):
    """ Class for retrieving handler classes for auxiliary data
    (mss, sic, sitype, snow). The classes are initialized with directory
    information from the local machine definition and the auxdata information
    from `auxdata.yaml` configuration file.
    """

    def __init__(self):
        super(DefaultAuxdataClassHandler, self).__init__(self.__class__.__name__)
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
                error_message = f"Missing entry `auxdata_repository.{auxdata_class}.{auxdata_id}` in " + \
                    f"local_machine_def ({psrlcfg.local_machine_def_filepath})"
                self.error.add_error(error_id, error_message)
                self.error.raise_on_error()
            empty_str = len(local_repo) == 0 if local_repo is not None else False
            if empty_str:
                msg = "Path definition for {}.{} exists in local_machine_def.yaml, but is empty string"
                msg = msg.format(auxdata_class, auxdata_id)
                logger.warning(msg)
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
        module_name, class_name = f"pysiral.auxdata.{auxdata_class}", auxdata_def["pyclass"]
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
        aux_repo_defs = psrlcfg.local_machine.auxdata_repository
        try:
            local_repo_auxclass = aux_repo_defs[auxdata_class]
        except KeyError:
            local_repo_auxclass = {}
            msg = "Missing auxdata definition in local_machine_def.yaml: auxdata_repository.%s" % auxdata_class
            self.error.add_error("missing-localmachinedef-tag", msg)
            self.error.raise_on_error()
        return local_repo_auxclass.get(auxdata_id, None)

    def get_auxdata_def(self, auxdata_class: str, auxdata_id: str) -> "AttrDict":
        """
        Returns the definition in `config/auxdata_def.yaml` for specified auxdata class and id.
        Raises an error if the entry is not found.
        :param auxdata_class: The code for auxiliary data type (sic, mss, sitype, snow, ...)
        :param auxdata_id: The id of a specific data set for the auxiliary data class (e.g. sic:osisaf-operational)
        :return: The configuration dictionary
        """

        auxdata_def = psrlcfg.auxdef.get_definition(auxdata_class, auxdata_id)
        if auxdata_def is None:
            msg = f"Cannot find entry for auxiliary data set {auxdata_class}:{auxdata_id} in auxdata_def.yaml"
            self.error.add_error("invalid-auxdata-class", msg)
            self.error.raise_on_error()
        return auxdata_def.attrdict


class L1PDataHandler(DefaultLoggingClass):
    """ Class for querying L1P data files """

    def __init__(self,
                 platform: str,
                 hemisphere: str,
                 source_version: str = None,
                 file_version: str = None,
                 ):
        """
        Init the class
        :param platform: pysiral compliant platform id (cryosat2, envisat, ...)
        :param hemisphere: hemishphere code (north, nh, ...)
        :param source_version: Input version, e.g. CryoSat-2 baseline-d, ...
        :param: file_version: The file version of the l1p data (e.g. v1p1).
            NOTE: Only newer l1p files use the version subfolder. Thus, this is optional.
        """
        super(L1PDataHandler, self).__init__(self.__class__.__name__)

        # Save args
        self._platform = platform
        self._hemisphere = hemisphere
        self._source_version = source_version
        self._file_version = file_version if file_version is not None else self._autodetect_file_version()
        self._last_directory = None

    def get_files_from_time_range(self, time_range: DatePeriod) -> List[str]:
        """
        Query l1p files for a a given time range.
        :param time_range: a dateperiods.DatePeriod instance
        :return:
        """

        # Validate time_range (needs to be of type DatePeriod)
        if not isinstance(time_range, DatePeriod):
            error = ErrorStatus()
            msg = "Invalid type of time_range, required: dateperiods.DatePeriod, was %s" % (type(time_range))
            error.add_error("invalid-timerange-type", msg)
            error.raise_on_error()

        # 1) get list of all files for monthly folders
        yyyy, mm = "%04g" % time_range.tcs.year, "%02g" % time_range.tcs.month
        directory = Path(self.l1p_base_dir)
        if self._file_version is not None:
            directory = directory / self._file_version
        directory = directory / self._hemisphere / yyyy / mm
        all_l1p_files = sorted(list(directory.rglob("*.nc")))

        # 3) Check if files are in requested time range
        # This serves two purposes: a) filter out files with timestamps that do
        # not belong in the directory. b) get a subset if required
        l1p_filepaths = [l1p_file for l1p_file in all_l1p_files if self.l1p_in_trange(l1p_file, time_range)]

        # Save last search directory
        self._last_directory = directory

        # Done
        return l1p_filepaths

    def _autodetect_file_version(self) -> Union[str, None]:
        """
        Autodetect l1p file version from the directory structure
        :return:
        """

        # Get all sub-folders of the l1p directory
        all_subfolders = [d.name for d in Path(self.l1p_base_dir).iterdir() if d.is_dir()]

        # Filter hemisphere codes from the sub-folder names. These will be present as sub-folders
        # if the file version is not part of the directory structure
        invalid_versions = ["north", "south", "global", "nh", "sh"]
        file_versions = [d for d in all_subfolders if d not in invalid_versions]

        if len(file_versions) == 0:
            logger.info("l1p file version autodetect: No file version subfolder")
            return None
        elif len(file_versions) == 1:
            logger.info(f"l1p file version autodetect: {file_versions[0]}")
            return file_versions[0]
        else:
            msg = f"l1p file version autodetect failure: Multiple versions found [{file_versions}]"
            msg = msg + " -> specify file_version"
            self.error.add_error("l1p-input-error", msg)
            self.error.raise_on_error()
        return None

    @staticmethod
    def l1p_in_trange(fn: str, tr: DatePeriod) -> bool:
        """ Returns flag if filename is within time range """
        # Parse infos from l1bdata filename
        fnattr = PysiralOutputFilenaming()
        fnattr.parse_filename(fn)
        # Compute overlap between two start/stop pairs
        is_overlap = fnattr.start <= tr.tce.dt and fnattr.stop >= tr.tcs.dt
        return is_overlap

    @property
    def l1p_base_dir(self) -> str:
        """
        Returns the the l1p base base
        # TODO: Evaluate use of function (properties shouldn't raise errors?)
        :return:
        """
        l1p_base_dir = None
        try:
            l1p_base_dir = psrlcfg.local_machine.l1b_repository[self._platform][self._source_version]["l1p"]
        except (AttributeError, KeyError):
            msg = f"missing entry `l1bdata.{self._platform}.{self._source_version}.l1p` in local_machine_def.yaml"
            msg = msg + f" ({psrlcfg.local_machine_def_filepath})"
            error_id = "l1bdata_missing_localrepo_def"
            self.error.add_error(error_id, msg)
            self.error.raise_on_error()
        return l1p_base_dir

    @property
    def last_directory(self) -> str:
        return str(self._last_directory)


class L2iDataHandler(DefaultLoggingClass):
    """ Class for retrieving default l1b directories and filenames """

    def __init__(self, base_directory, force_l2i_subfolder=True, search_str="l2i"):
        super(L2iDataHandler, self).__init__(self.__class__.__name__)
        self.error = ErrorStatus(caller_id=self.__class__.__name__)
        self._base_directory = base_directory
        self._force_l2i_subfolder = force_l2i_subfolder
        self._subdirectory_list = self.get_subdirectory_list()
        self._search_str = search_str
        self._validate_base_directory()

    def get_files_from_time_range(self, time_range):
        """ Get all files that fall into time range (May be spread over
        the different year/ month subfolders """
        l2i_files = []

        for daily_periods in time_range.get_segments("day"):
            year, month, day = daily_periods.tcs.year, daily_periods.tcs.month, daily_periods.tcs.day
            lookup_directory = self.get_lookup_directory(year, month)
            if not Path(lookup_directory).is_dir():
                continue
            l2i_pattern = self.get_l2i_search_str(year=year, month=month, day=day)
            result = list(Path(lookup_directory).glob(l2i_pattern))
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
        l2i_files = list(l2i_files)

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
            msg = "Invalid %s product directory: %s"
            msg = msg % (self._search_str, str(self._base_directory))
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

        # 1. Check if directory exists
        if not self.product_basedir.is_dir():
            msg = "Directory {} does not exist".format(str(self.product_basedir))
            self.error.add_error("invalid-l2i-product-dir", msg)
            self.error.raise_on_error()

        try:
            years = sorted([f.parts[-1] for f in Path(self.product_basedir).iterdir() if f.is_dir()])
        except StopIteration:
            logger.warning("No subdirectories in %s" % self.product_basedir)
            return []
        # filter any invalid directories
        years = [y for y in years if re.match(r'[1-3][0-9]{3}', y)]
        for year in years:
            subdir_year = Path(self.product_basedir) / year
            months = [f.parts[-1] for f in subdir_year.iterdir() if f.is_dir()]
            # filter any invalid directories
            months = [m for m in months if re.match(r'[0-1][0-9]', m)]
            subdirectory_list.extend([[year, m] for m in months])
        return subdirectory_list

    @staticmethod
    def get_l2i_search_str(year: int = None,
                           month: int = None,
                           day: int = None) -> str:
        """
        Returns a search pattern for l2i files with optional refined search for year, month, day.
        NOTE: month & day can only be set, if the year & year + month respectively is set

        Examples:
            *l2i*.nc
            *l2i*2017*.nc
            *l2i*201704*.nc
            *l2i*20170401*.nc
        :param year:
        :param month:
        :param day:
        :return:
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
        l2i_file_pattern = "%s.nc" % date_str
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
