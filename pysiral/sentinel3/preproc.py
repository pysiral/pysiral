# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 20:52:32 2016

@author: Stefan
"""

from pysiral.sentinel3.iotools import Sentinel3FileList
from pysiral.l1bdata import L1bConstructor
from pysiral.config import ConfigInfo
from pysiral.helper import (get_first_array_index, get_last_array_index, rle)
from pysiral.path import validate_directory
from pysiral.output import (l1bnc_filenaming, L1bDataNC)
from pysiral.logging import stdout_logger
from pysiral.flag import FlagContainer

import numpy as np


class Sentinel3PreProc(object):

    def __init__(self):
        self.month = None
        self.year = None
        self.mission = None
        self.skip = 0

    def execute(self):

        """ Get the logging instance """
        log = stdout_logger("sentinel3-preproc")

        """ Get the configuration data for handling Sentinel3X data """
        config = ConfigInfo()
        l1b_repository = config.local_machine.l1b_repository
        s3_l1b_repository = l1b_repository[self.mission].sral

        """ Get the list of files """
        # for the case of this test only for one month of data
        s3_files = Sentinel3FileList()
        s3_files.log = log
        s3_files.folder = s3_l1b_repository
        s3_files.year = self.year
        s3_files.month = self.month
        s3_files.search()

        """ Start the CryoSat-2 pre-processor """
        job = S3PreProcJob(self.mission)
        job.config = config
        job.log = log
        job.files = s3_files.sorted_list[self.skip:]
        job.merge_and_export_polar_ocean_subsets()


class S3PreProcJob(object):

    def __init__(self, mission):
        self.log = None
        self.config = ConfigInfo()
        self.options = None
        self.files = []
        self.mission = mission
        self.new_segment_delta_seconds = 10
        self._get_default_options()

    def merge_and_export_polar_ocean_subsets(self):
        """ loop over remaining files in file list """
        log = self.log
        n = len(self.files)
        if n == 0:
            return
        l1bdata_stack = []
        for i, sentinel3_l2_file in enumerate(self.files):

            # Parse the current file and split into polar ocean segments
            log.info("+ Parsing file %g of %g: %s" % (
                i+1, n, sentinel3_l2_file))
            l1b_segments = self.get_sentinel3_l1bdata_ocean_segments(
                sentinel3_l2_file, self.config)

            # Skip if no relevant data was found
            if l1b_segments is None:
                log.info("- no polar ocean data, skipping")
                continue
            log.info("- %g polar ocean data segments" % len(l1b_segments))

            # XXX: Debug
            # debug_stack_orbit_plot(l1bdata_stack, l1b_segments)

            # Loop over segments and check connectivity
            for l1b_segment in l1b_segments:
                if not self.l1b_is_connected_to_stack(
                        l1b_segment, l1bdata_stack):
                    # => break criterium, save existing stack and start over
                    log.info("- segment unconnected, exporting current stack")
                    self.concatenate_and_export_stack_files(l1bdata_stack)
                    # Reset the l1bdata stack
                    l1bdata_stack = [l1b_segment]
                    continue
                # polar ocean data and connected => add to stack
                l1bdata_stack.append(l1b_segment)

        """ Export the final stack => done """
        self.concatenate_and_export_stack_files(l1bdata_stack)

    def get_sentinel3_l1bdata_ocean_segments(self, filename, config):
        """
        Returns the source ERS l1b data as a list of
        pysiral.L1bdata objects.
        """

        # Read the sentinel3 file
        l1b = L1bConstructor(config)
        l1b.mission = self.mission
        l1b.filename = filename
        l1b.get_header_info()

        self.log.info(
            "- sar_mode_percent: %.1f" % l1b.info.sar_mode_percent)
        self.log.info(
            "- open_ocean_percent: %.1f" % l1b.info.open_ocean_percent)

        # check sar mode percentage and open ocean data
        if l1b.info.sar_mode_percent <= 0.:
            self.log.info("- no sar mode data, skipping file")
            return None

        # check sar mode percentage and open ocean data
        if l1b.info.open_ocean_percent <= 0.:
            self.log.info("- no open ocean data, skipping file")
            return None

        # check sar mode percentage and open ocean data
        if l1b.error_status:
            self.log.info("- l1b error status raised, skipping file")
            return None

        # Read the full data set
        l1b.construct()

        # Extract relevant segments over ocean
        l1b_single_file_subsets = self.extract_polar_segments(l1b)

        ocean_segments = []
        for l1b_single_file_subset in l1b_single_file_subsets:
            l1b_ocean_sections = self.extract_polar_ocean_segments(
                l1b_single_file_subset)
            ocean_segments.extend(l1b_ocean_sections)

        return ocean_segments

    def extract_polar_segments(self, l1b):
        """
        Sentinel3 early data specific: There might be gaps in the data files
        with no SAR coverage
        """

        # The threshold is the same for arctic and antarctic (50°)
        polar_threshold = self.options.polar_threshold

        is_polar = FlagContainer(
            np.abs(l1b.time_orbit.latitude > polar_threshold))
        polar_subset = l1b.extract_subset(is_polar.indices)

        from matplotlib.dates import date2num

        timestamp = polar_subset.time_orbit.timestamp
        delta_seconds = np.diff(86400.*date2num(timestamp))
        threshold = self.options.max_connected_files_timedelta_seconds
        new_segments_start_index = np.where(delta_seconds > threshold)[0]+1
        n_segments = len(new_segments_start_index)+1

        # 1 segment -> return full orbit
        if n_segments == 1:
            return [polar_subset]

        # Create two lists describing the start and end indeces of
        # orbit segments
        segments_start_index = [0]
        segments_start_index.extend(new_segments_start_index)
        segments_stop_index = [index-1 for index in new_segments_start_index]
        segments_stop_index.append(polar_subset.n_records-1)

        # Extract the individual data segments
        l1bdata_segments = []
        for i in np.arange(n_segments):
            subset_list = np.arange(segments_start_index[i],
                                    segments_stop_index[i]+1)
            l1bdata_segment = polar_subset.extract_subset(subset_list)
            l1bdata_segments.append(l1bdata_segment)

        return l1bdata_segments

#        import matplotlib.pyplot as plt
#
#        plt.figure()
#        xref = np.arange(l1b.n_records)
#        plt.plot(xref, timestamp)
#        plt.scatter(xref[new_segments_index], timestamp[new_segments_index])
#
#        plt.figure()
#        plt.plot(np.diff(86400.*date2num(l1b.time_orbit.timestamp)))
#        plt.show()
#
#        stop

    def extract_polar_ocean_segments(self, l1b):
        """
        Extract segments of continous ocean data not separated by larger
        landmasses.
        """

        # 1) Trim l1b data from the first to the last ocean data record
        self.log.info("- trim outer non-ocean regions")
        ocean = l1b.surface_type.get_by_name("ocean")
        first_ocean_index = get_first_array_index(ocean.flag, True)
        last_ocean_index = get_last_array_index(ocean.flag, True)
        if first_ocean_index is None or last_ocean_index is None:
            return None
        n = l1b.info.n_records-1
        is_full_ocean = first_ocean_index == 0 and last_ocean_index == n
        if not is_full_ocean:
            ocean_subset = np.arange(first_ocean_index, last_ocean_index+1)
            self.log.info("- ocean subset [%g:%g]" % (
                 first_ocean_index, last_ocean_index))
            l1b.trim_to_subset(ocean_subset)

        # 2) Identify larger landmasses and split orbit into segments
        #    if necessary
        self.log.info("- test for inner larger landmasses")
        ocean = l1b.surface_type.get_by_name("ocean")
        not_ocean_flag = np.logical_not(ocean.flag)
        segments_len, segments_start, not_ocean = rle(not_ocean_flag)
        landseg_index = np.where(not_ocean)[0]
        self.log.info("- number of landmasses: %g" % len(landseg_index))
        if len(landseg_index) == 0:
            # no land segements, return single segment
            return [l1b]
        treshold = self.options.max_inner_nonocean_segment_nrecords
        large_landsegs_index = np.where(
            segments_len[landseg_index] > treshold)[0]
        large_landsegs_index = landseg_index[large_landsegs_index]
        self.log.info("- number of large landmasses: %g" % len(
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
            self.log.debug("- ocean segment: [%g:%g]" % (
                start_index, stop_index))
            l1b_segments.append(l1b.extract_subset(subset_list))
            start_index = segments_start[index+1]
        # extract the last subset
        self.log.debug("- ocean segment: [%g:%g]" % (
            start_index, len(ocean.flag)-1))
        last_subset_list = np.arange(start_index, len(ocean.flag))
        l1b_segments.append(l1b.extract_subset(last_subset_list))
        return l1b_segments

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
        l1bdata_stack.pop(0)
        for orbit_segment in l1bdata_stack:
            l1b_merged.append(orbit_segment)

        # Prepare data export
        config = self.config
        export_folder, export_filename = l1bnc_filenaming(l1b_merged, config)
        log.info("+ Creating l1bdata netCDF: %s" % export_filename)
        log.info("- in folder: %s" % export_folder)
        validate_directory(export_folder)

        # Export the data object
        ncfile = L1bDataNC()
        ncfile.l1b = l1b_merged
        ncfile.config = config
        ncfile.output_folder = export_folder
        ncfile.filename = export_filename
        ncfile.export()

    def _get_default_options(self):
        self.options = self.config.mission[self.mission].preproc.options


def debug_stack_orbit_plot(l1b_stack, l1b_segments):

    from mpl_toolkits.basemap import Basemap
    import matplotlib.pyplot as plt

    grid_keyw = {"dashes": (None, None), "color": "#bcbdbf",
                 "linewidth": 0.5, "latmax": 88}

    gridb_keyw = {"dashes": (None, None), "color": "#003e6e",
                  "linewidth": 2, "latmax": 88}

#    lon_0 = l1b_segments[-1].info.lon_max
#    lat_0 = l1b_segments[-1].info.lat_max

    plt.figure("Stack Debug Map", figsize=(12, 6), facecolor="#ffffff")
    m = Basemap(projection='cyl', resolution='l')
    # m = Basemap(projection='ortho', lon_0=lon_0, lat_0=lat_0, resolution='l')
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
