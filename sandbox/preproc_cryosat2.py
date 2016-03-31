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

file where strange things are happening:

over southern ocean but no has_ocean => false
D:\awi\altim\data\altimetry\cryosat2\baseline-c\SIR_SAR_L1\2015\03\...
 ...CS_LTA__SIR_SAR_1B_20150301T025222_20150301T025526_C001.DBL

potential wrong lat/lon:
D:\awi\altim\data\altimetry\cryosat2\baseline-c\SIR_SAR_L1\2015\03\...
 ...CS_LTA__SIR_SAR_1B_20150301T032811_20150301T033251_C001.DBL

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
log = Logger('cryosat2-preproc')


""" Definitions """
POLAR_THRESHOLD = 50.
SEARCH_LENGTH_INITIAL_STACK_FILE = 10
MAX_CONNECTED_FILES_TIMEDELTA_SECONDS = 10
MAX_INNER_NONOCEAN_SEGMENT_NRECORDS = 1000

# jump to specific file indices for debug purposes
test_continous_landmass_index = 20


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

    """ loop over remaining files in file list """
    n = len(cryosat2_files.sorted_list)
    l1bdata_stack = []
    for i, cryosat2_l1b_file in enumerate(cryosat2_files.sorted_list):
        if i < 15:
            continue
        # Parse the current file and split into polar ocean segments
        log.info("Parsing file %g of %g: %s" % (i+1, n, cryosat2_l1b_file))
        l1b_segments = get_cryosat2_l1bdata_ocean_segments(
            cryosat2_l1b_file, config)
        # Skip if no relevant data was found
        if l1b_segments is None:
            log.info("no polar ocean data, skipping ...")
            continue
        log.info("%g polar ocean data segments" % len(l1b_segments))
        # XXX: Debug
        debug_stack_orbit_plot(l1bdata_stack, l1b_segments)
        # Loop over segments and check connectivity
        for l1b_segment in l1b_segments:
            if not l1b_is_connected_to_stack(l1b_segment, l1bdata_stack):
                # => break criterium, save existing stack and start over
                log.info("segement unconnected, exporting current stack")
                concatenate_and_export_stack_files(l1bdata_stack)
                # Reset the l1bdata stack
                l1bdata_stack = [l1b_segment]
                continue
            # polar ocean data and connected => add to stack
            l1bdata_stack.append(l1b_segment)
    """ Export the final stack => done """
    concatenate_and_export_stack_files(l1bdata_stack)


def get_cryosat2_l1bdata_ocean_segments(filename, config):
    """
    Returns the source CryoSat-2 l1b data as a list of
    pysiral.L1bdata objects.

    """
    l1b = L1bConstructor(config)
    l1b.mission = "cryosat2"
    l1b.filename = filename
    l1b.construct()
    # Extract relevant segments over ocean
    l1b_list = extract_polar_ocean_segments(l1b)
    return l1b_list


def extract_polar_ocean_segments(l1b):
    """
    Extract segments of continous ocean data not separated by larger
    landmasses

    """
    # 1) Check if file has ocean data in the polar regions at all
    has_polar_ocean = region_is_arctic_or_antarctic_ocean(l1b)
    if not has_polar_ocean:
        return None

    # 2) Trim l1b data to polar region (known that list won't be emtpy)
    #    CryoSat-2 specific: l1b data does not cover both polar regions
    log.info("... trim to polar ocean subset")
    is_polar = np.abs(l1b.time_orbit.latitude) >= POLAR_THRESHOLD
    polar_subset = np.where(is_polar)[0]
    l1b.trim_to_subset(polar_subset)

    # 3) Trim l1b data from the first to the last ocean data record
    log.info("... trim outer non-ocean regions")
    ocean = l1b.surface_type.get_by_name("ocean")
    first_ocean_index = get_first_array_index(ocean.flag, True)
    last_ocean_index = get_last_array_index(ocean.flag, True)
    n = l1b.info.n_records-1
    is_full_ocean = first_ocean_index == 0 and last_ocean_index == n
    if not is_full_ocean:
        ocean_subset = np.arange(first_ocean_index, last_ocean_index+1)
        log.info("... ocean subset [%g:%g]" % (
            first_ocean_index, last_ocean_index))
        l1b.trim_to_subset(ocean_subset)

    # 4) Identify larger landmasses and split orbit into segments if necessary
    log.info("... test for inner larger landmasses")
    ocean = l1b.surface_type.get_by_name("ocean")
    not_ocean_flag = np.logical_not(ocean.flag)
    segments_len, segments_start, not_ocean = rle(not_ocean_flag)
    landseg_index = np.where(not_ocean)[0]
    log.info("... number of landmasses: %g" % len(landseg_index))
    if len(landseg_index) == 0:
        # no land segements, return single segment
        return [l1b]
    large_landsegs_index = np.where(
        segments_len[landseg_index] > MAX_INNER_NONOCEAN_SEGMENT_NRECORDS)[0]
    large_landsegs_index = landseg_index[large_landsegs_index]
    log.info("... number of large landmasses: %g" % len(large_landsegs_index))
    if len(large_landsegs_index) == 0:
        # no large land segments, return single segment
        return [l1b]
    # Large land segments exist, split the l1b data object
    # first and last is always ocean, since already trimmed before
    # (see step 3)
    # start
    l1b_segments = []
    start_index = 0
    for index in large_landsegs_index:
        stop_index = segments_start[index]
        subset_list = np.arange(start_index, stop_index)
        log.debug("... ocean segment: [%g:%g]" % (start_index, stop_index))
        l1b_segments.append(l1b.extract_subset(subset_list))
        start_index = segments_start[index+1]
    # extract the last subset
    log.debug("... ocean segment: [%g:%g]" % (start_index, len(ocean.flag)-1))
    last_subset_list = np.arange(start_index, len(ocean.flag))
    l1b_segments.append(l1b.extract_subset(last_subset_list))
    return l1b_segments


def region_is_arctic_or_antarctic_ocean(l1b):
    """
    Test if a CryoSat-2 l1b file has data over ocean in either Arctic
    of Antartic
    """
    # 1) test if either minimum or maximum latitude is in polar regions
    lat_range = np.abs([l1b.info.lat_min, l1b.info.lat_max])
    log.debug("... lat_min: %14.7f" % l1b.info.lat_min)
    log.debug("... lat_max: %14.7f" % l1b.info.lat_max)
    is_polar = np.amax(lat_range) >= POLAR_THRESHOLD
    log.info("... is_in_polar_region: %s" % str(is_polar))
    # 2) test if there is any ocean data at all
    ocean = l1b.surface_type.get_by_name("ocean")
    has_ocean = ocean.num > 0
    log.info("... has_ocean: %s" % str(has_ocean))
    # Return a flag and the ocean flag list for later use
    return is_polar and has_ocean


def l1b_is_connected_to_stack(l1b, l1b_stack):
    """
    Check if the start time of file i and the stop time if file i-1
    indicate neighbouring orbit segments (e.g. due to radar mode change)

    """
    if len(l1b_stack) == 0:  # Stack is empty
        return True
    last_l1b_in_stack = l1b_stack[-1]
    timedelta = l1b.info.start_time - last_l1b_in_stack.info.stop_time
    is_connected = timedelta.seconds <= MAX_CONNECTED_FILES_TIMEDELTA_SECONDS
    log.debug("... last_l1b_in_stack.info.stop_time: %s" %
              str(last_l1b_in_stack.info.stop_time))
    log.debug("... l1b.info.start_time: %s" % str(l1b.info.start_time))
    log.info("... is connected: %s (timedelta = %g sec)" % (
        str(is_connected), timedelta.seconds))
    return is_connected


#def has_large_landmass_in_between(l1b):
#    """
#    Test if content of file is of type ocean - landmass - ocean
#    and the landmass is of a certain size. In this case, the file need to
#    splitted in a separate process.
#
#    """
#    # last index must be ocean
#    ocean = l1b.surface_type.get_by_name("ocean")
#    if not ocean.flag[-1]:
#        return False
#    # get land flag (for simplicity everything which is not ocean)
#    not_ocean = np.logical_not(ocean.flag)
#
#
#    first_land_occurance = list(is_land).index(True)
#    last_land_
#    stop


def get_first_array_index(array, value):
    """ Get the index in array of the first occurance of ``value`` """
    try:
        index = list(array).index(value)
    except:
        index = None
    return index


def get_last_array_index(array, value):
    """ Get the index in array of the last occurance of ``value`` """
    listarray = list(array)
    try:
        index = (len(listarray) - 1) - listarray[::-1].index(value)
    except:
        index = None
    return index


def rle(inarray):
    """
    run length encoding. Partial credit to R rle function.
    Multi datatype arrays catered for including non Numpy
    returns: tuple (runlengths, startpositions, values)

    from: http://stackoverflow.com/questions/1066758/find-length-of-sequences-
                 of-identical-values-in-a-numpy-array
    """
    ia = np.array(inarray)                   # force numpy
    n = len(ia)
    if n == 0:
        return (None, None, None)
    else:
        y = np.array(ia[1:] != ia[:-1])      # pairwise unequal (string safe)
        i = np.append(np.where(y), n - 1)    # must include last element posi
        z = np.diff(np.append(-1, i))        # run lengths
        p = np.cumsum(np.append(0, z))[:-1]  # positions
        return(z, p, ia[i])


def concatenate_and_export_stack_files(l1bdata_stack):
    log.debug("Length of l1bdata_stack: %g" % len(l1bdata_stack))
    debug_stack_export_orbit_plot(l1bdata_stack)


def debug_stack_orbit_plot(l1b_stack, l1b_segments):

    from mpl_toolkits.basemap import Basemap
    import matplotlib.pyplot as plt

    grid_keyw = {"dashes": (None, None), "color": "#bcbdbf",
                 "linewidth": 0.5, "latmax": 88}

    gridb_keyw = {"dashes": (None, None), "color": "#003e6e",
                  "linewidth": 2, "latmax": 88}

    lon_0 = l1b_segments[-1].info.lon_max
    lat_0 = l1b_segments[-1].info.lat_max

    plt.figure("Stack Debug Map", figsize=(12, 12), facecolor="#ffffff")
    m = Basemap(projection='ortho', lon_0=lon_0, lat_0=lat_0, resolution='l')
    m.fillcontinents(color='#00ace5', lake_color='#00ace5')

    # draw parallels and meridians.
    m.drawparallels(np.arange(-90., 120., 10.), **grid_keyw)
    m.drawparallels(np.arange(-50., 51., 100.), **gridb_keyw)
    m.drawmeridians(np.arange(0., 420., 30.), **grid_keyw)

    # Draw the contents of the l1b stack
    for l1b in l1b_stack:
        x, y = m(l1b.time_orbit.longitude, l1b.time_orbit.latitude)
        m.plot(x, y, color="#003e6e", linewidth=2.0)
        m.scatter(x[0], y[0], color="#003e6e")

    # Draw the segments from the current l1b file
    for l1b in l1b_segments:
        ocean = l1b.surface_type.get_by_name("ocean")
        ocean_list = np.where(ocean.flag)[0]
        x, y = m(l1b.time_orbit.longitude, l1b.time_orbit.latitude)
        m.plot(x, y, color="#aa0000", linewidth=2.0, zorder=100)
        m.scatter(x[ocean_list], y[ocean_list], color="#76FF7A", s=15,
                  zorder=99)
        m.scatter(x[0], y[0], color="#aa0000", zorder=101, s=20)
    plt.show(block=True)


def debug_stack_export_orbit_plot(l1b_stack):

    from mpl_toolkits.basemap import Basemap
    import matplotlib.pyplot as plt

    grid_keyw = {"dashes": (None, None), "color": "#bcbdbf",
                 "linewidth": 0.5, "latmax": 88}

    gridb_keyw = {"dashes": (None, None), "color": "#003e6e",
                  "linewidth": 2, "latmax": 88}

    plt.figure("Stack Export Debug Map", figsize=(18, 9), facecolor="#ffffff")
    m = Basemap(projection='cyl', resolution='l')
    m.fillcontinents(color='#00ace5', lake_color='#00ace5')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-90., 120., 10.), **grid_keyw)
    m.drawparallels(np.arange(-50., 51., 100.), **gridb_keyw)
    m.drawmeridians(np.arange(0., 420., 30.), **grid_keyw)
    for l1b in l1b_stack:
        x, y = m(l1b.time_orbit.longitude, l1b.time_orbit.latitude)
        m.scatter(x, y, color="#ff0000", s=2, zorder=100)
        m.scatter(x[0], y[0], color="#aa0000", s=20, zorder=101)

    plt.show(block=True)


class CryoSat2FileListAllModes(object):
    """
    Class for the construction of a list of CryoSat-2 SAR/SIN files
    sorted by acquisition time
    """

    def __init__(self):
        self.folder_sar = None
        self.folder_sin = None
        self.year = None
        self.month = None
        self.pattern = ".DBL"
        self._list = deque([])
        self._sorted_list = None

    def search(self):
        self._search_specific_mode_files("sar")
        self._search_specific_mode_files("sin")
        self._sort_mixed_mode_file_list()

    @property
    def sorted_list(self):
        return [item[0] for item in self._sorted_list]

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
        return (os.path.join(dirpath, filename), filename.split("_")[6])

    def _sort_mixed_mode_file_list(self):
        log.info("Sorting files")
        dtypes = [('path', object), ('start_time', object)]
        self._sorted_list = np.array(self._list, dtype=dtypes)
        self._sorted_list.sort(order='start_time')
        log.info("... done")


if __name__ == "__main__":
    preproc_cryosat2()
