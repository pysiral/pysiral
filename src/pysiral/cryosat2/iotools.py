# -*- coding: utf-8 -*-


import re
from collections import deque
from pathlib import Path

import numpy as np
from loguru import logger
from typing import List, Dict, Union
from parse import parse

from pydantic import BaseModel


class FileDiscoveryConfig(BaseModel):

    local_machine_def_tag: str
    lookup_modes: List[str]
    # Prefilled by pysiral.l1preproc.Level1PreProcJobDef._get_local_input_directory
    # TODO: Move here during Level-1 Preprocessor refactoring
    lookup_dir: Dict[str, Union[str, Path]]
    baseline: str
    platform: str = "cryosat2"

    def filename_search(self, radar_mode) -> str:
        return f"CS_*_SIR_{radar_mode.upper()}_*1B_{{year:04d}}{{month:02d}}{{day:02d}}*_{self.baseline}*.nc"

    @property
    def filename_parser(self) -> str:
        return r"CS_{data_record_type}__SIR_{radar_mode}_{processing_level}_{time_coverage_start}_{time_coverage_end}_{baseline}{file_version}.nc"


class ESACryoSat2ICEL1bProductsFileDiscovery(object):

    def __init__(self, cfg_dict: dict) -> None:
        """ Initialize the file discovery class with a configuration dictionary """

        # Save config
        self.cfg = FileDiscoveryConfig(**cfg_dict)

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

    def _append_files(self, mode, period) -> None:
        lookup_year, lookup_month = period.tcs.year, period.tcs.month
        lookup_dir = self._get_lookup_dir(lookup_year, lookup_month, mode)
        logger.info("Search directory: %s" % lookup_dir)
        n_files = 0
        for daily_period in period.get_segments("day"):
            # Search for specific day
            year, month, day = daily_period.tcs.year, daily_period.tcs.month, daily_period.tcs.day
            file_list = self._get_files_per_day(lookup_dir, year, month, day, mode)
            tcs_list = self._get_tcs_from_filenames(file_list)
            n_files += len(file_list)
            for file, tcs in zip(file_list, tcs_list):
                self._list.append((file, tcs))
        logger.info(" Found %g %s files" % (n_files, mode))

    def _get_files_per_day(self, lookup_dir, year, month, day, mode):
        """ Return a list of files for a given lookup directory """
        # Search for specific day

        file_name_template = self.cfg.filename_search(mode)
        filename_search = file_name_template.format(year=year, month=month, day=day)
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
            # filename_segments = re.split(r"_+|\.", str(Path(filename).name))
            result = parse(self.cfg.filename_parser, str(Path(filename).name))
            if result is not None:
                tcs.append(result['time_coverage_start'])
            else:
                logger.warning(f"Could not parse filename: {filename}")
        return tcs

    @property
    def sorted_list(self):
        dtypes = [('path', object), ('start_time', object)]
        self._sorted_list = np.array(self._list, dtype=dtypes)
        self._sorted_list.sort(order='start_time')
        return [item[0] for item in self._sorted_list]
