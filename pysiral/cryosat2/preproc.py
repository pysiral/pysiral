# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 20:52:32 2016

@author: Stefan
"""

from pysiral.cryosat2.iotools import CryoSat2FileListAllModes
from pysiral.l1bdata import L1bConstructor
from pysiral.config import ConfigInfo
from pysiral.helper import (get_first_array_index, get_last_array_index, rle)
from pysiral.path import validate_directory
from pysiral.output import (l1bnc_filenaming, L1bDataNC)
from pysiral.logging import stdout_logger

import numpy as np


class CryoSat2PreProc(object):

    def __init__(self):
        self.month = None
        self.year = None

    def execute(self):

        """ Get the logging instance """
        log = stdout_logger("cryosat2-preproc")

        """ Get the configuration data for handling CryoSat-2 data """
        config = ConfigInfo()
        cryosat_l1b_repository = config.local_machine.l1b_repository.cryosat2

        """ Get the list of files (SAR and SIN in chronological order) """
        # for the case of this test only for one month of data
        cryosat2_files = CryoSat2FileListAllModes()
        cryosat2_files.log = log
        cryosat2_files.folder_sar = cryosat_l1b_repository.sar
        cryosat2_files.folder_sin = cryosat_l1b_repository.sin
        cryosat2_files.year = self.year
        cryosat2_files.month = self.month
        cryosat2_files.search()

        """ Start the CryoSat-2 pre-processor """
        job = CryoSat2PreProcJob()
        job.config = config
        job.log = log
        job.files = cryosat2_files.sorted_list
        job.merge_and_export_polar_ocean_subsets()


class CryoSat2PreProcJob(object):

    def __init__(self):
        self.log = None
        self.config = ConfigInfo()
        self.options = None
        self.files = []
        self._get_default_options()

    def merge_and_export_polar_ocean_subsets(self):
        """ loop over remaining files in file list """
        log = self.log
        n = len(self.files)
        if n == 0:
            return
        l1bdata_stack = []
        for i, cryosat2_l1b_file in enumerate(self.files):
#            if i < 20:
#                continue
            # Parse the current file and split into polar ocean segments
            log.info("Parsing file %g of %g: %s" % (
                i+1, n, cryosat2_l1b_file))
            l1b_segments = self.get_cryosat2_l1bdata_ocean_segments(
                cryosat2_l1b_file, self.config)
            # Skip if no relevant data was found
            if l1b_segments is None:
                log.info("no polar ocean data, skipping ...")
                continue
            log.info("%g polar ocean data segments" % len(l1b_segments))
            # XXX: Debug
            # debug_stack_orbit_plot(l1bdata_stack, l1b_segments)
            # Loop over segments and check connectivity
            for l1b_segment in l1b_segments:
                if not self.l1b_is_connected_to_stack(
                        l1b_segment, l1bdata_stack):
                    # => break criterium, save existing stack and start over
                    log.info("segement unconnected, exporting current stack")
                    self.concatenate_and_export_stack_files(l1bdata_stack)
                    # Reset the l1bdata stack
                    l1bdata_stack = [l1b_segment]
                    continue
                # polar ocean data and connected => add to stack
                l1bdata_stack.append(l1b_segment)

        """ Export the final stack => done """
        self.concatenate_and_export_stack_files(l1bdata_stack)

    def get_cryosat2_l1bdata_ocean_segments(self, filename, config):
        """
        Returns the source CryoSat-2 l1b data as a list of
        pysiral.L1bdata objects.

        """
        l1b = L1bConstructor(config)
        l1b.mission = "cryosat2"
        l1b.filename = filename
        l1b.construct()
        # Extract relevant segments over ocean
        l1b_list = self.extract_polar_ocean_segments(l1b)
        return l1b_list

    def extract_polar_ocean_segments(self, l1b):
        """
        Extract segments of continous ocean data not separated by larger
        landmasses

        """
        log = self.log
        # 1) Check if file has ocean data in the polar regions at all
        has_polar_ocean = self.region_is_arctic_or_antarctic_ocean(l1b)
        if not has_polar_ocean:
            return None

        # 2) Trim l1b data to polar region (known that list won't be emtpy)
        #    CryoSat-2 specific: l1b data does not cover both polar regions
        log.info("... trim to polar ocean subset")
        polar_threshold = self.options.polar_threshold
        is_polar = np.abs(l1b.time_orbit.latitude) >= polar_threshold
        polar_subset = np.where(is_polar)[0]
        l1b.trim_to_subset(polar_subset)

        # 3) Trim l1b data from the first to the last ocean data record
        log.info("... trim outer non-ocean regions")
        ocean = l1b.surface_type.get_by_name("ocean")
        first_ocean_index = get_first_array_index(ocean.flag, True)
        last_ocean_index = get_last_array_index(ocean.flag, True)
        if first_ocean_index is None or last_ocean_index is None:
            return None
        n = l1b.info.n_records-1
        is_full_ocean = first_ocean_index == 0 and last_ocean_index == n
        if not is_full_ocean:
            ocean_subset = np.arange(first_ocean_index, last_ocean_index+1)
            log.info("... ocean subset [%g:%g]" % (
                 first_ocean_index, last_ocean_index))
            l1b.trim_to_subset(ocean_subset)

        # 4) Identify larger landmasses and split orbit into segments
        #    if necessary
        log.info("... test for inner larger landmasses")
        ocean = l1b.surface_type.get_by_name("ocean")
        not_ocean_flag = np.logical_not(ocean.flag)
        segments_len, segments_start, not_ocean = rle(not_ocean_flag)
        landseg_index = np.where(not_ocean)[0]
        log.info("... number of landmasses: %g" % len(landseg_index))
        if len(landseg_index) == 0:
            # no land segements, return single segment
            return [l1b]
        treshold = self.options.max_inner_nonocean_segment_nrecords
        large_landsegs_index = np.where(
            segments_len[landseg_index] > treshold)[0]
        large_landsegs_index = landseg_index[large_landsegs_index]
        log.info("... number of large landmasses: %g" % len(
            large_landsegs_index))
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
        log.debug("... ocean segment: [%g:%g]" % (
            start_index, len(ocean.flag)-1))
        last_subset_list = np.arange(start_index, len(ocean.flag))
        l1b_segments.append(l1b.extract_subset(last_subset_list))
        return l1b_segments

    def region_is_arctic_or_antarctic_ocean(self, l1b):
        """
        Test if a CryoSat-2 l1b file has data over ocean in either Arctic
        of Antartic
        """
        log = self.log
        # 1) test if either minimum or maximum latitude is in polar regions
        lat_range = np.abs([l1b.info.lat_min, l1b.info.lat_max])
        polar_threshold = self.options.polar_threshold
        is_polar = np.amax(lat_range) >= polar_threshold
        log.info("... is_in_polar_region: %s" % str(is_polar))
        # 2) test if there is any ocean data at all
        ocean = l1b.surface_type.get_by_name("ocean")
        has_ocean = ocean.num > 0
        log.info("... has_ocean: %s" % str(has_ocean))
        # Return a flag and the ocean flag list for later use
        return is_polar and has_ocean

    def l1b_is_connected_to_stack(self, l1b, l1b_stack):
        """
        Check if the start time of file i and the stop time if file i-1
        indicate neighbouring orbit segments (e.g. due to radar mode change)

        """
        if len(l1b_stack) == 0:  # Stack is empty
            return True
        last_l1b_in_stack = l1b_stack[-1]
        timedelta = l1b.info.start_time - last_l1b_in_stack.info.stop_time
        threshold = self.options.max_connected_files_timedelta_seconds
        is_connected = timedelta.seconds <= threshold
    #    log.debug("... last_l1b_in_stack.info.stop_time: %s" %
    #              str(last_l1b_in_stack.info.stop_time))
    #    log.debug("... l1b.info.start_time: %s" % str(l1b.info.start_time))
    #    log.info("... is connected: %s (timedelta = %g sec)" % (
    #        str(is_connected), timedelta.seconds))
        return is_connected

    def concatenate_and_export_stack_files(self, l1bdata_stack):
        log = self.log
        # log.debug("Length of l1bdata_stack: %g" % len(l1bdata_stack))
        # debug_stack_export_orbit_plot(l1bdata_stack)
        # Concatenate the files
        l1b_merged = l1bdata_stack[0]
        if l1b_merged.info.radar_mode == "sin":
            l1b_merged.reduce_waveform_bin_count(256)
        l1bdata_stack.pop(0)
        for orbit_segment in l1bdata_stack:
            # Rebin
            if orbit_segment.info.radar_mode == "sin":
                orbit_segment.reduce_waveform_bin_count(256)
            l1b_merged.append(orbit_segment)
        # stop
        # Prepare data export
        config = self.config
        export_folder, export_filename = l1bnc_filenaming(l1b_merged, config)
        log.info("Create file: %s" % export_filename)
        log.info("Folder: %s" % export_folder)
        validate_directory(export_folder)
        # Export the data object
        ncfile = L1bDataNC()
        ncfile.l1b = l1b_merged
        ncfile.config = config
        ncfile.output_folder = export_folder
        ncfile.filename = export_filename
        ncfile.export()

    def _get_default_options(self):
        self.options = self.config.mission.cryosat2.preproc.options


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
