# -*- coding: utf-8 -*-

import os
import numpy as np

from collections import deque
from pysiral.logging import DefaultLoggingClass


class EnvisatFileList(DefaultLoggingClass):
    """
    Class for the construction of a list of Envisat N1 files
    sorted by acquisition time
    XXX: Currently only support order by month/date and not cycle
    """

    def __init__(self):

        super(EnvisatFileList, self).__init__(self.__class__.__name__)

        self.folder = None
        self.year = None
        self.month = None
        self.pattern = ".N1"
        self.day_list = None
        self.time_range = None
        self._list = deque([])
        self._sorted_list = []

    def search(self, time_range):

        # Only per month search possible at this moment
        self.year = time_range.start.year
        self.month = time_range.start.month

        # Create a list of day if not full month is required
        if not time_range.is_full_month:
            self.day_list = np.arange(
                time_range.start.day, time_range.stop.day+1)

        # List all files for the specific month
        self._get_file_listing()

        # Limit the date range (if necessary)
        self._limit_to_time_range()

    @property
    def sorted_list(self):
        return [item[0] for item in self._sorted_list]

    def _get_file_listing(self):
        search_toplevel_folder = self._get_toplevel_search_folder()
        # walk through files
        for dirpath, dirnames, filenames in os.walk(search_toplevel_folder):

            self.log.info("Searching folder: %s" % dirpath)

            # Get the list of all .N1 files
            sgdr_files = [fn for fn in filenames if self.pattern in fn]
            sgdr_files = sorted(sgdr_files)

            # Report total number of files
            self.log.info("Found %g level-1b SGDR files" % len(sgdr_files))

            sgdr_list = [self._get_list_item(fn, dirpath) for fn in sgdr_files]
            self._sorted_list.extend(sgdr_list)

    def _get_toplevel_search_folder(self):
        folder = self.folder
        if self.year is not None:
            folder = os.path.join(folder, "%4g" % self.year)
        if self.month is not None and self.year is not None:
            folder = os.path.join(folder, "%02g" % self.month)
        return folder

    def _get_list_item(self, filename, dirpath):
        """ Return full path and date str """
        return (os.path.join(dirpath, filename),
                filename.split("_")[2].split("-")[1][1:])

    def _limit_to_time_range(self):

        # self.day_list is only set if time_range is not a full month
        if self.day_list is None:
            return

        # Cross-check the data label and day list
        self._sorted_list = [fn for fn in self._sorted_list if
                             int(fn[1][-2:]) in self.day_list]

        self.log.info("%g files match time range of this month" % (
            len(self._sorted_list)))
