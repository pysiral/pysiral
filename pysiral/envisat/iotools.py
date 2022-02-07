# -*- coding: utf-8 -*-

import os
import numpy as np
from loguru import logger
from pathlib import Path
from typing import Tuple, List
from collections import deque

from dateperiods import DatePeriod
from pysiral.errorhandler import ErrorStatus
from pysiral.logging import DefaultLoggingClass


class EnvisatSGDRNC(DefaultLoggingClass):

    def __init__(self, cfg) -> None:
        """
        File discovery for Envisat SGDG netcdf files with yyyy/mm/dd subfolder structure
        :param cfg:
        """

        cls_name = self.__class__.__name__
        super(EnvisatSGDRNC, self).__init__(cls_name)
        self.error = ErrorStatus(caller_id=cls_name)

        # Save config
        self.cfg = cfg

        # Init empty file lists
        self._reset_file_list()

    def get_file_for_period(self, period: "DatePeriod") -> List[str]:
        """
        Query for Sentinel Level-2 files for a specific period.
        :param period: dateperiods.DatePeriod
        :return: sorted list of filenames
        """
        # Make sure file list are empty
        self._reset_file_list()
        self._query(period)
        return self.sorted_list

    def _reset_file_list(self) -> None:
        """ Resets the result of previous file searches """
        self._list = deque([])
        self._sorted_list = []

    def _query(self, period: "DatePeriod") -> None:
        """
        Searches for files in the given period and stores result in property _sorted_list
        :param period: dateperiods.DatePeriod
        :return: None
        """

        # Loop over all months in the period
        daily_periods = period.get_segments("day")
        for daily_period in daily_periods:

            year, month, day = daily_period.tcs.year, daily_period.tcs.month, daily_period.tcs.day

            # Search folder
            lookup_folder = self._get_lookup_folder(year, month, day)
            if not Path(lookup_folder).is_dir():
                continue

            # Query the daily folder
            filename_search = self.cfg.filename_search.format(year=year, month=month, day=day)
            sgdr_files = list(Path(lookup_folder).glob(filename_search))

            # Add files to result list
            if not sgdr_files:
                continue
            self._sorted_list.extend(sorted(sgdr_files))

    def _get_lookup_folder(self, year, month, day) -> Path:
        yyyy, mm, dd = "%04g" % year, "%02g" % month, "%02g" % day
        return Path(self.cfg.lookup_dir) / yyyy / mm / dd

    @property
    def sorted_list(self) -> List[str]:
        return list(self._sorted_list)
