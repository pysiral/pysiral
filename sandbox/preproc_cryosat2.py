# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 18:54:15 2016

@author: Stefan

Testbed for CryoSat-2 specific pre-processor that
- creates l1bncfiles from the l1b source files
- concatenates adjacent files (e.g. at borders of SAR and SIN modes)
- trims to ocean data in Arctic and Antarctic

Filenaming:
    l1bdat_v$VERS_$REG_$MISSION_$ORBIT_$YYYYMMDDHHMISS_$YYYYMMDDHHMISS.nc

$VERS            l1bdata version       00 (beta)
$REG             region                [north | south]
$MISSION         mission short name    => cryosat2
$ORBIT           orbit/cylce number    e.g. 00026000
$YYYYMMDDHHMISS  start and end time

Definition of Arctic:
    latitude >= 50.0

Definition of Antarctic:
    latitude <= -50.0


Concatinating orbit files (CryoSat-2 specific):

- Loop over L1b files (month wise, sar and sin)
- parse l1b source files
- check if ocean content
- start a stack of l1bdata objectes
- check of adjacent to previous file
    yes: add to l1bdata stack
    no: - concatenate files in stack
        - write to l1bdata nc file
        - clear stack
        - reopen stack with current file

"""

from pysiral.config import ConfigInfo
from pysiral.l1bdata import L1bConstructor

from collections import deque
from operator import itemgetter
from logbook import Logger, StreamHandler

import numpy as np
import os
import sys


StreamHandler(sys.stdout).push_application()
log = Logger('Logbook')


""" Definitions """
POLAR_THRESHOLD = 50.


def preproc_cryosat2():

    """ Get the configuration data for handling CryoSat-2 data """
    config = ConfigInfo()
    cryosat_l1b_repository = config.local_machine.l1b_repository.cryosat2

    """ Get the list of files (SAR and SIN in chronological order) """
    # for the case of this test only for one month of data
    cryosat2_files = CryoSat2FileListAllModes()
    cryosat2_files.folder_sar = cryosat_l1b_repository.sar
    cryosat2_files.folder_sin = cryosat_l1b_repository.sin
    cryosat2_files.year = 2015
    cryosat2_files.month = 3
    cryosat2_files.search()

    """ Start the pre-processing """
    cryosat2_file_list = cryosat2_files.sorted_list
    # Initialize the stack with the first file
    l1bdata_stack = []
    # now loop over all files
    for cryosat2_l1b_file in cryosat2_file_list:
        # get the current file
        log.info("Parsing file: %s" % cryosat2_l1b_file)
        l1b = get_cryosat2_l1bdata(cryosat2_l1b_file, config)
        log.info("... done")
        # test if data is in polar regions at all
        if not region_is_arctic_or_antarctic_ocean(l1b):
            continue
        # test if current l1b data objects needs to be connected to
        # previous one


def get_cryosat2_l1bdata(filename, config):
    l1b = L1bConstructor(config)
    l1b.mission = "cryosat2"
    l1b.filename = filename
    l1b.construct()
    return l1b


def region_is_arctic_or_antarctic_ocean(l1b):
    """
    Test if a CryoSat-2 l1b file has data over ocean in either Arctic
    of Antartic
    """
    # 1) test if either minimum or maximum latitude is in polar regions
    lat_range = np.abs([l1b.info.lat_min, l1b.info.lat_max])
    is_polar = np.amin(lat_range) >= POLAR_THRESHOLD
    log.info("... is_polar: %s" % str(is_polar))
    # 2) test if there is any ocean data at all
    ocean_flag = l1b.surface_type.get_by_name("ocean")
    has_ocean = ocean_flag.num > 0
    log.info("... has_ocean: %s" % str(is_polar))
    return is_polar and has_ocean


class CryoSat2FileListAllModes(object):

    def __init__(self):
        self.folder_sar = None
        self.folder_sin = None
        self.year = None
        self.month = None
        self.pattern = ".DBL"
        self._list = deque([])

    def search(self):
        self._search_specific_mode_files("sar")
        self._search_specific_mode_files("sin")
        self._sort_mixed_mode_file_list()

    @property
    def sorted_list(self):
        return [item[0] for item in self._list]

    def _search_specific_mode_files(self, mode):
        search_toplevel_folder = self._get_toplevel_search_folder(mode)
        # walk through files
        for dirpath, dirnames, filenames in os.walk(search_toplevel_folder):
            log.info("Search Folder: %s" % dirpath)
            # Get the list of all dbl files
            cs2files = [fn for fn in filenames if self.pattern in fn]
            log.info("Found %g CryoSat-2 l1b files" % len(cs2files))
            # reform the list that each list entry is of type
            # [full_path, identifier (start_date)] for later sorting
            # of SAR and SIN files
            sublist = [self._get_list_item(fn, dirpath) for fn in cs2files]
            self._list.extend(sublist)

    def _get_toplevel_search_folder(self, mode):
        folder = getattr(self, "folder_"+mode)
        if self.year is not None:
            folder = os.path.join(folder, "%4g" % self.year)
        if self.month is not None:
            folder = os.path.join(folder, "%02g" % self.month)
        return folder

    def _get_list_item(self, filename, dirpath):
        return [os.path.join(dirpath, filename), filename.split("_")[6]]

    def _sort_mixed_mode_file_list(self):
        log.info("Sorting files")
        self._list = sorted(self._list, key=itemgetter(1))
        log.info("... done")


if __name__ == "__main__":
    preproc_cryosat2()
