# -*- coding: utf-8 -*-

import os
import re
from pathlib import Path
from collections import deque

from pysiral.errorhandler import ErrorStatus
from pysiral.core import DefaultLoggingClass


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
