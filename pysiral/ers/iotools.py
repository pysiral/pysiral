# -*- coding: utf-8 -*-

import os
import numpy as np

from collections import deque
from pysiral.logging import DefaultLoggingClass


class ERSFileList(DefaultLoggingClass):
    """
    Class for the construction of a list of Envisat N1 files
    sorted by acquisition time
    XXX: Currently only support order by month/date and not cycle
    """

    def __init__(self):
        super(ERSFileList, self).__init__(self.__class__.__name__)
        self.folder = None
        self.year = None
        self.month = None
        self.time_range = None
        self.day_list = []
        self.pattern = ".NC"
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

        self._get_file_listing()

    @property
    def sorted_list(self):
        return self._sorted_list

    def _get_file_listing(self):
        search_toplevel_folder = self._get_toplevel_search_folder()
        # walk through files
        for dirpath, dirnames, filenames in os.walk(search_toplevel_folder):

            # Envisat file structure specific:
            # additional folders for days, therefore ignore top level folder
            if dirpath == search_toplevel_folder:
                continue

            # Get the day from the directory name
            current_day = int(os.path.split(dirpath)[-1])

            # Check if in day list
            if current_day not in self.day_list:
                continue

            self.log.info("Searching folder: %s" % dirpath)
            # Get the list of all dbl files
            files = [os.path.join(dirpath, fn) for fn in filenames
                     if self.pattern in fn]
            self.log.info("Found %g level-1b SGDR files" % len(files))
            # reform the list that each list entry is of type
            # [full_path, identifier (start_date)] for later sorting
            # of SAR and SIN files
            sublist = sorted(files)
            self._sorted_list.extend(sublist)

    def _get_toplevel_search_folder(self):
        folder = self.folder
        if self.year is not None:
            folder = os.path.join(folder, "%4g" % self.year)
        if self.month is not None and self.year is not None:
            folder = os.path.join(folder, "%02g" % self.month)
        return folder
