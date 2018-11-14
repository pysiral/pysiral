# -*- coding: utf-8 -*-

import os

import numpy as np
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
            self.day_list = np.arange(
                time_range.start.day, time_range.stop.day+1)

        # Search all sar/sin product files per month
        self._search_specific_mode_files("sar")
        self._search_specific_mode_files("sin")

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

            self.log.info("Searching folder: %s" % dirpath)

            # Get the list of all dbl files
            cs2files = [fn for fn in filenames if self.pattern in fn]
            self.log.info("Found %g %s level-1b files" % (
                len(cs2files), mode))

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
        self._sorted_list = [fn for fn in self._sorted_list if
                             int(fn[1][6:8]) in self.day_list]

        self.log.info("%g files match time range of this month" % (
            len(self._sorted_list)))

    def _get_toplevel_search_folder(self, mode):
        folder = getattr(self, "folder_"+mode)
        if self.year is not None:
            folder = os.path.join(folder, "%4g" % self.year)
        if self.month is not None:
            folder = os.path.join(folder, "%02g" % self.month)
        return folder

    def _get_list_item(self, filename, dirpath):
        return (os.path.join(dirpath, filename), filename.split("_")[6])

    def _sort_mixed_mode_file_list(self):
        dtypes = [('path', object), ('start_time', object)]
        self._sorted_list = np.array(self._list, dtype=dtypes)
        self._sorted_list.sort(order='start_time')


class BaselineDFileDiscovery(DefaultLoggingClass):

    def __init__(self, cfg):
        cls_name = self.__class__.__name__
        super(BaselineDFileDiscovery, self).__init__(cls_name)
        self.error = ErrorStatus(caller_id=cls_name)

        # Init empty file lists
        self._reset_file_list()

    def get_file_for_period(self, period):
        """ Return a list of sorted files """

        # Make sure file list are empty
        self._reset_file_list()

        for mode in self.cfg.lookup_modes:
            self._add_input_files(period, mode)

    def _reset_file_list(self):
        self._list = deque([])
        self._sorted_list = None

    def _append_files(self, mode, period):

        for year, month, day in period.dayslist:
            lookup_dir = self._get_lookup_dir(year, month, mode)

    def _get_lookup_dir(self, year, month, mode):
        return None


    # @property