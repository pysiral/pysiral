# -*- coding: utf-8 -*-

import os
import re
import numpy as np
from pathlib import Path
from loguru import logger
from collections import deque

from pysiral.errorhandler import ErrorStatus
from pysiral.logging import DefaultLoggingClass


class CryoSat2MonthlyFileListAllModes(DefaultLoggingClass):
    """
    Class for the construction of a list of CryoSat-2 SAR/SIN files
    sorted by acquisition time
    """

    def __init__(self):

        name = self.__class__.__name__
        super(CryoSat2MonthlyFileListAllModes, self).__init__(name)

        self.folder_sar = None
        self.folder_sin = None
        self.year = None
        self.month = None
        self.day_list = None
        self.time_range = None
        self.pattern = ".DBL"
        self._list = deque([])
        self._sorted_list = None

    def search(self, time_range):

        self.year = time_range.start.year
        self.month = time_range.start.month

        # Create a list of day if not full month is required
        if not time_range.is_full_month:
            self.day_list = np.arange(time_range.start.day, time_range.stop.day+1)

        # Search all sar/sin product files per month
        for radar_mode in ["sar", "sin"]:
            self._search_specific_mode_files(radar_mode)

        # Sort sar/sin files by acquisition date
        self._sort_mixed_mode_file_list()

        # Limit the date range (if necessary)
        self._limit_to_time_range()

    @property
    def sorted_list(self):
        return [item[0] for item in self._sorted_list]

    def _search_specific_mode_files(self, mode):

        search_toplevel_folder = self._get_toplevel_search_folder(mode)

        # walk through files
        for dirpath, dirnames, filenames in os.walk(search_toplevel_folder):

            logger.info("Searching folder: %s" % dirpath)

            # Get the list of all dbl files
            cs2files = [fn for fn in filenames if self.pattern in fn]
            logger.info("Found %g %s level-1b files" % (len(cs2files), mode))

            # reform the list that each list entry is of type
            # [full_path, identifier (start_date)] for later sorting
            # of SAR and SIN files
            sublist = [self._get_list_item(fn, dirpath) for fn in cs2files]
            self._list.extend(sublist)

    def _limit_to_time_range(self):

        # self.day_list is only set if time_range is not a full month
        if self.day_list is None:
            return

        # Cross-check the data label and day list
        self._sorted_list = [fn for fn in self._sorted_list if int(fn[1][6:8]) in self.day_list]

        logger.info("%g files match time range of this month" % (len(self._sorted_list)))

    def _get_toplevel_search_folder(self, mode):
        folder = Path(getattr(self, "folder_"+mode))
        if self.year is not None:
            folder = folder / "{:04g}".format(self.year)
        if self.month is not None:
            folder = folder / "{:02g}".format(self.month)
        return folder

    @staticmethod
    def _get_list_item(filename, dirpath):
        return Path(dirpath) / filename, filename.split("_")[6]

    def _sort_mixed_mode_file_list(self):
        dtypes = [('path', object), ('start_time', object)]
        self._sorted_list = np.array(self._list, dtype=dtypes)
        self._sorted_list.sort(order='start_time')


class BaselineDFileDiscovery(DefaultLoggingClass):

    def __init__(self, cfg):
        cls_name = self.__class__.__name__
        super(BaselineDFileDiscovery, self).__init__(cls_name)
        self.error = ErrorStatus(caller_id=cls_name)

        # Save config
        self.cfg = cfg

        # Properties
        self._sorted_list = []

        # Init empty file lists
        self._reset_file_list()

    def get_file_for_period(self, period):
        """ Return a list of sorted files """
        # Make sure file list are empty
        self._reset_file_list()
        for mode in self.cfg.lookup_modes:
            self._append_files(mode, period)
        return self.sorted_list

    def _reset_file_list(self):
        self._list = deque([])
        self._sorted_list = []

    def _append_files(self, mode, period):
        lookup_year, lookup_month = period.tcs.year, period.tcs.month
        lookup_dir = self._get_lookup_dir(lookup_year, lookup_month, mode)
        logger.info("Search directory: %s" % lookup_dir)
        n_files = 0
        for daily_period in period.get_segments("day"):
            # Search for specific day
            year, month, day = daily_period.tcs.year, daily_period.tcs.month, daily_period.tcs.day
            file_list = self._get_files_per_day(lookup_dir, year, month, day)
            tcs_list = self._get_tcs_from_filenames(file_list)
            n_files += len(file_list)
            for file, tcs in zip(file_list, tcs_list):
                self._list.append((file, tcs))
        logger.info(" Found %g %s files" % (n_files, mode))

    def _get_files_per_day(self, lookup_dir, year, month, day):
        """ Return a list of files for a given lookup directory """
        # Search for specific day
        filename_search = self.cfg.filename_search.format(year=year, month=month, day=day)
        return sorted(Path(lookup_dir).glob(filename_search))

    def _get_lookup_dir(self, year, month, mode):
        yyyy, mm = "%04g" % year, "%02g" % month
        return Path(self.cfg.lookup_dir[mode]) / yyyy / mm

    def _get_tcs_from_filenames(self, files):
        """
        Extract the part of the filename that indicates the time coverage start (tcs)
        :param files: a list of files
        :return: tcs: a list with time coverage start strings of same length as files
        """
        tcs = []
        for filename in files:
            filename_segments = re.split(r"_+|\.", str(Path(filename).name))
            tcs.append(filename_segments[self.cfg.tcs_str_index])
        return tcs

    @property
    def sorted_list(self):
        dtypes = [('path', object), ('start_time', object)]
        self._sorted_list = np.array(self._list, dtype=dtypes)
        self._sorted_list.sort(order='start_time')
        return [item[0] for item in self._sorted_list]
