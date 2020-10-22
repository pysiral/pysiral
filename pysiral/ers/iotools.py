# -*- coding: utf-8 -*-


import os
import re
import numpy as np

import dateutil
from pathlib import Path
from loguru import logger
from datetime import timedelta
from parse import compile

from collections import deque
from pysiral.errorhandler import ErrorStatus
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
                time_range.start.day, time_range.stop.day + 1)

        self._get_file_listing()

    @property
    def sorted_list(self):
        return self._sorted_list

    def _get_file_listing(self):
        """
        Look for files
        :return:
        """
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

            logger.info("Searching folder: %s" % dirpath)
            # Get the list of all dbl files
            files = [os.path.join(dirpath, fn) for fn in filenames
                     if self.pattern in fn]
            logger.info("Found %g level-1b SGDR files" % len(files))
            # reform the list that each list entry is of type
            # [full_path, identifier (start_date)] for later sorting
            # of SAR and SIN files
            sublist = sorted(files)
            self._sorted_list.extend(sublist)

    def _get_toplevel_search_folder(self):
        folder = Path(self.folder)
        if self.year is not None:
            folder = folder / "%4g" % self.year
        if self.month is not None and self.year is not None:
            folder = folder / "%02g" % self.month
        return folder


class ERSCycleBasedSGDR(DefaultLoggingClass):

    def __init__(self, cfg):
        """
        File discovery for a cycle based order
        :param cfg:
        """

        cls_name = self.__class__.__name__
        super(ERSCycleBasedSGDR, self).__init__(cls_name)
        self.error = ErrorStatus(caller_id=cls_name)

        # Save config
        self.cfg = cfg

        # Establish a lookup table that maps cycles to days from the cycle folder names
        self._create_date_lookup_table()

        # Init empty file lists
        self._reset_file_list()

    def get_file_for_period(self, period):
        """
        Query for Sentinel Level-2 files for a specific period.
        :param period: dateperiods.DatePeriod
        :return: sorted list of filenames
        """
        # Make sure file list are empty
        self._reset_file_list()
        self._query(period)
        return self.sorted_list

    def _create_date_lookup_table(self):
        """
        Get a look up table based on the folder names for each cycle. The content will be stored in
        self._date_lookup_table
        :return: None
        """

        # --- Parameters ---
        lookup_dir = self.cfg.lookup_dir  # The main repository directory
        regex = self.cfg.cycle_folder_regex  # The regex expression used to identify cycle folders
        folder_parser = compile(self.cfg.folder_parser)  # a parser string indicating parameters in folder name

        # --- Find all cycle folders ---
        # Note: the regex only applies to the subfolder name, not the full path
        cycle_folders = [x[0] for x in os.walk(lookup_dir) if re.match(regex, os.path.split(x[0])[-1])]

        # --- Construct lookup table from each subfolder name ---
        self._lookup_table = []
        for cycle_folder in cycle_folders:

            # Parse the folder name
            result = folder_parser.parse(os.path.split(cycle_folder)[-1])

            # Get start and end coverage as datetimes
            tcs, tce = dateutil.parser.parse(result["tcs"]), dateutil.parser.parse(result["tce"])

            # Compute list of dates between two datetimes
            delta = tce - tcs
            dates = []
            for i in range(delta.days + 1):
                date = tcs + timedelta(days=i)
                dates.append("%04g-%02g-%02g" % (date.year, date.month, date.day))

            # Add entry to lookup table
            result_dict = {}
            for item in result.named.items():
                result_dict[item[0]] = item[1]
            self._lookup_table.append(dict(dir=cycle_folder, dates=dates, **result_dict))

    def _reset_file_list(self):
        """ Resets the result of previous file searches """
        self._list = deque([])
        self._sorted_list = []

    def _query(self, period):
        """
        Searches for files in the given period and stores result in property _sorted_list
        :param period: dateperiods.DatePeriod
        :return: None
        """

        # Loop over all months in the period
        daily_periods = period.get_segments("day")
        for daily_period in daily_periods:

            year, month, day = daily_period.tcs.year, daily_period.tcs.month, daily_period.tcs.day

            # The date format in the lookup table is yyyy-mm-dd
            datestr = "%04g-%02g-%02g" % (year, month, day)

            # Get a list of cycle folders
            cycle_folders = self._get_cycle_folders_from_lookup_table(datestr)
            if len(cycle_folders) > 2:
                raise IOError("Date %s in more than 2 cycle folders (this should not happen)" % datestr)

            # Query each cycle folder
            filename_search = self.cfg.filename_search.format(year=year, month=month, day=day)
            for cycle_folder in cycle_folders:
                sgdr_files = Path(cycle_folder).glob(filename_search)
                self._sorted_list.extend(sorted(sgdr_files))

    def _get_cycle_folders_from_lookup_table(self, datestr):
        """
        Return a list of cycle folders that contain a date. Should be between 0 and 2
        :param datestr:
        :return:
        """

        cycle_folders = []
        for entry in self._lookup_table:
            if datestr in entry["dates"]:
                cycle_folders.append(entry["dir"])
        return cycle_folders

    @property
    def sorted_list(self):
        return list(self._sorted_list)
