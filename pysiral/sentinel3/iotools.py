# -*- coding: utf-8 -*-

import os
import re
import pandas as pd
from dataclasses import dataclass, field
from collections import deque
from datetime import datetime
from pathlib import Path

from dateperiods import DatePeriod
from loguru import logger
from parse import parse

from pysiral.core import DefaultLoggingClass
from pysiral.core.clocks import StopWatch
from pysiral.core.errorhandler import ErrorStatus


class Sentinel3FileList(DefaultLoggingClass):
    """
    Class for the construction of a list of Sentinel-3 SRAL L2 files
    sorted by acquisition time
    XXX: Based on the file and directory structure of the early access data
    """

    def __init__(self):

        super(Sentinel3FileList, self).__init__(self.__class__.__name__)
        self.folder = None
        self.time_range = None
        self.target = "enhanced_measurement.nc"
        self._sorted_list = []

    def search(self, time_range):
        """ Find all files falling in a defined time range """
        # Reset search result and save time range
        self._sorted_list = []
        self.time_range = time_range
        self._get_file_listing()

    @property
    def sorted_list(self):
        return list(self._sorted_list)

    def _get_file_listing(self):
        """ List all files in time range """

        monthly_periods = self.time_range.get_segments("month")
        for monthly_period in monthly_periods:

            year, month = monthly_period.tcs.year, monthly_period.tcs.month

            # Get the file list for each month
            toplevel_folder = self._get_toplevel_search_folder(year, month)
            l2_file_list = get_sentinel3_l1b_filelist(
                    toplevel_folder, self.target)

            # Get list of days for particular year/month
            days = self.time_range.get_days_for_month(year, month)
            for day in days:
                daystr = "%04g%02g%02g" % (year, month, day)
                match = [f[0] for f in l2_file_list if re.search(daystr, f[1])]
                self._sorted_list.extend(sorted(match))

    def _get_toplevel_search_folder(self, year, month):
        return Path(self.folder) / "%4g" % year / "%02g" % month


class CodaL2SralFileDiscovery(DefaultLoggingClass):
    """ Class to retrieve Sentinel-3 SRAL files from the Copernicus Online Data Archive """

    def __init__(self, cfg):
        """

        :param cfg: dict/treedict configuration options (see l1proc config file)
        """
        cls_name = self.__class__.__name__
        super(CodaL2SralFileDiscovery, self).__init__(cls_name)
        self.error = ErrorStatus(caller_id=cls_name)

        # Save config
        self.cfg = cfg

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

    def _query(self, period):
        """
        Searches for files in the given period and stores result in property _sorted_list
        :param period: dateperiods.DatePeriod
        :return: None
        """

        # Loop over all months in the period
        monthly_periods = period.get_segments("month")
        for monthly_period in monthly_periods:

            year, month = monthly_period.tcs.year, monthly_period.tcs.month

            # Get the file list for each month
            toplevel_folder = self._get_toplevel_search_folder(year, month)
            l2_file_list = get_sentinel3_l1b_filelist(toplevel_folder, self.cfg.filename_search)

            # Get list of days for particular year/month
            days = period.get_days_for_month(year, month)
            for day in days:
                daystr = "%04g%02g%02g" % (year, month, day)
                match = [f[0] for f in l2_file_list if re.search(daystr, f[1])]
                self._sorted_list.extend(sorted(match))

    def _get_toplevel_search_folder(self, year, month):
        """ Get the folder for the file search """
        return Path(self.cfg.lookup_dir) / "%4g" % year / "%02g" % month

    def _reset_file_list(self):
        """ Resets the result of previous file searches """
        self._list = deque([])
        self._sorted_list = []

    @property
    def sorted_list(self):
        """ Return the search result """
        return self._sorted_list


class L2SeaIceFileDiscovery(object):
    """ Class to retrieve Sentinel-3 SRAL files from the Copernicus Online Data Archive """

    def __init__(self, cfg):
        """

        :param cfg: dict/treedict configuration options (see l1proc config file)
        """

        # Save config
        self.cfg = cfg

        # Create inventory
        logger.info(f"Sentinel-3 source directory: {cfg.lookup_dir}")
        timer = StopWatch().start()
        self.catalogue = self._get_dataset_catalogue()
        logger.info(f"Found {self.n_catalogue_files} files ({self.cfg.filename_search})")
        timer.stop()
        logger.debug(f"Created Sentinel-3 file catalogue in {timer.get_seconds():.04f} seconds")

        # Init empty file lists
        self._reset_file_list()

    def _get_dataset_catalogue(self) -> pd.DataFrame:
        """
        Create a catalogues of file properties for all Sentinel-3 files in the lookup directory

        :return: catalogue as pandas DataFrame
        """

        # Simple catalogue format
        # [(datetime, filepath), (datetime, filepath), ...]
        nc_filepaths = sorted(Path(self.cfg.lookup_dir).glob(f"**/{self.cfg.filename_search}"))
        return pd.DataFrame([S3FileNaming(nc_filepath) for nc_filepath in nc_filepaths])

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

    def _query(self, period: DatePeriod) -> None:
        """
        Searches for files in the given period and stores result in property _sorted_list
        :param period: dateperiods.DatePeriod
        :return: None
        """

        # Find all files within the target period
        subset_df = self.catalogue[
            (self.catalogue.time_coverage_start >= period.tcs.dt) &
            (self.catalogue.time_coverage_end <= period.tce.dt)
        ]

        # Find possible time coverage duplicates and only keep the one with the latest creation time
        # NOTE: This is a workaround for some issues found with data downloaded from
        #       the Copernicus Data Space Ecosystem. Any more eloquent solution should be
        #       implemented in the download process.
        subset_cleaned = subset_df.groupby("time_coverage_start").apply(lambda x: x.loc[x.creation_time.idxmax()])
        if subset_df.shape[0] != subset_cleaned.shape[0]:
            logger.warning(f"Removed {subset_df.shape[0] - subset_cleaned.shape[0]} duplicate(s) in time coverage")

        self._sorted_list = subset_cleaned["filepath"].tolist()

    def _reset_file_list(self):
        """ Resets the result of previous file searches """
        self._list = deque([])
        self._sorted_list = []

    @property
    def sorted_list(self):
        """ Return the search result """
        return self._sorted_list

    @property
    def n_catalogue_files(self) -> int:
        """
        Return the number of files in the file catalogue
        :return:
        """
        return len(self.catalogue)


def get_sentinel3_l1b_filelist(folder, target_nc_filename):
    """ Returns a list with measurement.nc files for given month """
    s3_file_list = []
    for root, dirs, files in os.walk(folder):
        for name in files:
            if name == target_nc_filename:
                # Get the start datetime from the folder name
                datestr = os.path.split(root)[-1].split("_")[7]
                s3_file_list.append((os.path.join(root, name), datestr))
    return s3_file_list


def get_sentinel3_sral_l1_from_l2(l2_filename, target="enhanced_measurement.nc"):
    """ Returns the corresponding sral l1 file to a given l2 filename """

    # XXX: This is based on the data structure of the early access
    #      for expert users

    # Step 1: Replace product tag in folder structure
    l1nc_filename = l2_filename.replace("SR_2_LAN", "SR_1_SRA")

    # folder name for L1 and L2 data files are different, need to replace
    # one date tag with asterisk and search for match

    # split the directories
    directories = Path(l1nc_filename).parent.parts

    # split the dates with asterisk, orbit number should provide unique
    # match for the test data set
    s3_orbit_dir = directories[-1]
    s3_orbit_dir_components = s3_orbit_dir.split("_")
    s3_orbit_dir_components[7] = "*"     # Start time
    s3_orbit_dir_components[8] = "*"     # Stop time
    s3_orbit_dir_components[9] = "*"     # product creation time
    s3_orbit_dir_components[10] = "*"    # some version code

    # Compile the folder again with search pattern
    search_s3_l1b_folder = "_".join(s3_orbit_dir_components)
    sral_l1_folder = Path(*directories[:-1]).glob(search_s3_l1b_folder)

    if len(sral_l1_folder) > 0:
        l1nc_filename = sral_l1_folder[0] / target
    else:
        return None
    return l1nc_filename


@dataclass
class S3FileNaming:
    """
    Deciphering the Sentinel-3 filenaming convention
    (source: Sentinel 3 PDGS File Naming Convention (EUM/LEO-SEN3/SPE/10/0070, v1D, 24 June 2016)
    """
    filepath: Path
    is_valid_file: bool = field(init=False)
    file_id: str = field(init=False)
    mission_id: str = field(init=False)
    data_source: str = field(init=False)
    processing_level: str = field(init=False)
    data_type_id: str = field(init=False)
    time_coverage_start: datetime = field(init=False)
    time_coverage_end: datetime = field(init=False)
    creation_time: datetime = field(init=False)
    instance_id: str = field(init=False)
    product_generation_center: str = field(init=False)
    product_platform: str = field(init=False)
    timeliness: str = field(init=False)
    baseline: str = field(init=False)
    extension: str = field(init=False)

    def __post_init__(self):

        # Define the filenaming convention
        filenaming_convention = "{mission_id:3}_{data_source:2}_{processing_level:1}_" \
                                "{data_type_id:6}_{time_coverage_start:15}_{time_coverage_end:15}_" \
                                "{creation_time:15}_{instance_id:17}_{product_generation_center:3}_" \
                                "{product_platform:1}_{timeliness:2}_{baseline:3}.{extension}"

        # Properties to be stored as strings
        str_keys = ["mission_id", "data_source", "processing_level", "data_type_id",
                    "instance_id", "product_generation_center", "product_platform",
                    "timeliness", "baseline", "extension"]

        # Properties to be stored as datetime objects
        dt_keys = ["time_coverage_start", "time_coverage_end", "creation_time"]

        # The file is the name of the last part of the path
        self.file_id = self.filepath.parent.parts[-1]

        # Parse the file id
        elements = parse(filenaming_convention, self.file_id)
        if elements is None:
            logger.error(f"{self.file_id} is not a valid sentinel3 filename [{filenaming_convention}")
            self.is_valid_file = False
            return
        self.is_valid_file = True

        # Assign the parsed string properties
        for key in str_keys:
            setattr(self, key, elements.named[key])

        # Assign the parsed datetime properties
        for key in dt_keys:
            setattr(self, key, datetime.strptime(elements.named[key], "%Y%m%dT%H%M%S"))
