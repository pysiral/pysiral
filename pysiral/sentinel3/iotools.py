# -*- coding: utf-8 -*-

from pysiral.logging import DefaultLoggingClass

import os
import glob
import numpy as np


class Sentinel3FileList(DefaultLoggingClass):
    """
    Class for the construction of a list of Sentinel-3 SRAL L2 files
    sorted by acquisition time
    XXX: Based on the file and directory structure of the early access data
    """

    def __init__(self):

        super(Sentinel3FileList, self).__init__(self.__class__.__name__)
        self.folder = None
        self.year = None
        self.month = None
        self.target = "measurement.nc"
        self._sorted_list = []

    def search(self, time_range):

        self.year = time_range.start.year
        self.month = time_range.start.month

        # Create a list of day if not full month is required
        if not time_range.is_full_month:
            self.day_list = np.arange(
                time_range.start.day, time_range.stop.day+1)

        # Create a list of all files in folder
        self._get_file_listing()

        # Limit the date range (if necessary)
        self._limit_to_time_range()

    @property
    def sorted_list(self):
        return [item[0] for item in self._sorted_list]

    def _get_file_listing(self):
        search_toplevel_folder = self._get_toplevel_search_folder()
        s3_l2_file_list = get_sentinel3_l1b_filelist(
            search_toplevel_folder, target=self.target)
        self._sorted_list = sorted(s3_l2_file_list)

    def _get_toplevel_search_folder(self):
        folder = self.folder
        if self.year is not None:
            folder = os.path.join(folder, "%4g" % self.year)
        if self.month is not None and self.year is not None:
            folder = os.path.join(folder, "%02g" % self.month)
        return folder

    def _limit_to_time_range(self):

        # self.day_list is only set if time_range is not a full month
        if self.day_list is None:
            return

        # Cross-check the data label and day list
        self._sorted_list = [fn for fn in self._sorted_list if
                             int(fn[1][6:8]) in self.day_list]

        self.log.info("%g files match time range of this month" % (
            len(self._sorted_list)))


def get_sentinel3_l1b_filelist(folder, target="measurement.nc"):
    """ Returns a list with measurement.nc files for given month """
    s3_l1b_file_list = []
    for root, dirs, files in os.walk(folder):
        for name in files:
            if name == target:
                # Get the start datetime from the folder name
                datestr = os.path.split(root)[-1].split("_")[7]
                s3_l1b_file_list.append((os.path.join(root, name), datestr))
    return s3_l1b_file_list


def get_sentinel3_sral_l1_from_l2(l2_filename, target="measurement.nc"):
    """ Returns the corresponding sral l1 file to a given l2 filename """

    # XXX: This is based on the data structure of the early access
    #      for expert users

    # Step 1: Replace product tag in folder structure
    l1nc_filename = l2_filename.replace("SR_2_LAN", "SR_1_SRA")

    # folder name for L1 and L2 data files are different, need to replace
    # one date tag with asterisk and search for match

    # split the directories
    folder, filename = os.path.split(l1nc_filename)
    directories = folder.split(os.sep)

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
    main_dir = os.sep.join(directories[:-1])
    sral_l1_search = os.path.join(main_dir, search_s3_l1b_folder)

    # Search and return first match. If no match is found return none
    sral_l1_folder = glob.glob(sral_l1_search)

    if len(sral_l1_folder) > 0:
        l1nc_filename = os.path.join(sral_l1_folder[0], target)
    else:
        return None
    return l1nc_filename
