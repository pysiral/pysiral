# -*- coding: utf-8 -*-
"""
Created on Wed Aug 03 17:16:56 2016

@author: shendric
"""

from pysiral.config import (ConfigInfo, DefaultCommandLineArguments,
                            TimeRangeRequest)
from pysiral.errorhandler import ErrorStatus
from pysiral.flag import FlagContainer
from pysiral.logging import DefaultLoggingClass
from pysiral.output import (l1bnc_filenaming, get_l1bdata_export_folder,
                            L1bDataNC)
from pysiral.path import validate_directory

import numpy as np
import argparse
import os
import glob


class L1bPreProc(DefaultLoggingClass):

    def __init__(self, name):

        super(DefaultLoggingClass, self).__init__(self.__class__.__name__)
        self.error = ErrorStatus()
        self._jobdef = None
        self._l1b_file_list = []
        self._pysiral_config = ConfigInfo()

    def settings(self, jobdef):
        """ jobdef needs to be of type L1bPreProcJobSettings """
        self._jobdef = jobdef

    def execute(self):
        self.merge_and_export_polar_ocean_subsets()

    def remove_old_l1bdata(self):
        """
        Remove all l1b data files in current l1bdata export directory
        defined by the job definition
        """

        # Test if job definition is present
        if self._jobdef is None:
            return

        # Get the l1bdata diretory
        export_folder = get_l1bdata_export_folder(
            self._pysiral_config,
            self._jobdef.mission_id,
            self._jobdef.input_version,
            self._jobdef.hemisphere,
            self._jobdef.time_range.start.year,
            self._jobdef.time_range.start.month)

        # Get list of netcdf files
        search_pattern = os.path.joint(export_folder, "*.nc")
        l1bdata_files = glob.glob(search_pattern)

        # Delete files
        self.log.info("Removing %g l1bdata files in %s" % (
            len(l1bdata_files), export_folder))
        for l1bdata_filename in l1bdata_files:
            os.remove(l1bdata_filename)

    def merge_and_export_polar_ocean_subsets(self):
        """ loop over remaining files in file list """
        log = self.log
        n = len(self._l1b_file_list)
        if n == 0:
            return
        l1bdata_stack = []
        for i, l1b_file in enumerate(self._l1b_file_list):

            # Parse the current file and split into polar ocean segments
            log.info("+ Parsing file %g of %g: %s" % (i+1, n, l1b_file))
            l1b_segments = self.get_l1bdata_ocean_segments(l1b_file)

            # Skip if no relevant data was found
            if l1b_segments is None:
                log.info("- no %s polar ocean data, skipping" % (
                     self._jobdef.hemisphere))
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

    def extract_polar_segments(self, l1b):
        """
        Envisat specific: Extract Arctic and Antarctic subsets and
        return both in the correct sequence
        """

        # The threshold is the same for arctic and antarctic (50°)
        polar_threshold = self.options.polar_threshold

        # Target hemisphere: north, south or global
        target_hemisphere = self._jobdef.hemisphere

        # Extract Arctic (if target hemisphere, else empty list)
        if target_hemisphere == "north" or target_hemisphere == "global":
            is_arctic = FlagContainer(
                l1b.time_orbit.latitude > polar_threshold)
            arctic_subset = l1b.extract_subset(is_arctic.indices)
            self.log.info("- Extracted Arctic subset (%g records)" %
                          is_arctic.num)
        else:
            arctic_subset = None

        # Extract Antarctic (if target hemisphere, else empty list)
        if target_hemisphere == "south" or target_hemisphere == "global":
            is_antarctic = FlagContainer(
                l1b.time_orbit.latitude < -1.*polar_threshold)
            antarctic_subset = l1b.extract_subset(is_antarctic.indices)
            self.log.info("- Extracted Antarctic subset (%g records)" % (
                is_antarctic.num))
        else:
            antarctic_subset = None

        # Segments may be empty, test all cases
        arctic_is_valid = True
        antarctic_is_valid = True
        try:
            arctic_subset.info.start_time
        except AttributeError:
            arctic_is_valid = False

        try:
            antarctic_subset.info.start_time
        except AttributeError:
            antarctic_is_valid = False

        # Return both unsorted, Empty segment will be discarded later
        if not arctic_is_valid or not antarctic_is_valid:
            return [arctic_subset, antarctic_subset]

        # Order in sequence depeding on start time
        if arctic_subset.info.start_time < antarctic_subset.info.start_time:
            return [arctic_subset, antarctic_subset]
        else:
            return [antarctic_subset, arctic_subset]

    def concatenate_and_export_stack_files(self, l1bdata_stack):

        # XXX: This is hard coded since CryoSat-2 is the only
        #      mission that requires this
        sar_bin_count = {"baseline-b": 128, "baseline-c": 256}

        log = self.log

        # log.debug("Length of l1bdata_stack: %g" % len(l1bdata_stack))
        # debug_stack_export_orbit_plot(l1bdata_stack)

        # Concatenate the files
        l1b_merged = l1bdata_stack[0]
        bin_count = sar_bin_count[l1b_merged.info.mission_data_version.lower()]

        if l1b_merged.radar_modes == "sin":
            l1b_merged.reduce_waveform_bin_count(bin_count)

        l1bdata_stack.pop(0)
        for orbit_segment in l1bdata_stack:
            # Rebin
            if orbit_segment.radar_modes == "sin":
                orbit_segment.reduce_waveform_bin_count(bin_count)
            l1b_merged.append(orbit_segment)
        # stop
        # Prepare data export
        config = self._pysiral_config
        export_folder, export_filename = l1bnc_filenaming(
            l1b_merged, config, self._jobdef.input_version)
        log.info("Creating l1bdata netCDF: %s" % export_filename)
        log.info(". in folder: %s" % export_folder)
        validate_directory(export_folder)

        # Export the data object
        ncfile = L1bDataNC()
        ncfile.l1b = l1b_merged
        ncfile.config = config
        ncfile.output_folder = export_folder
        ncfile.filename = export_filename
        ncfile.export()

    def region_is_arctic_or_antarctic_ocean(self, l1b):
        """
        Test if a l1b file has data over ocean in either Arctic or Antartic
        """
        # 1) test if either minimum or maximum latitude is in polar regions
        lat_range = np.abs([l1b.info.lat_min, l1b.info.lat_max])
        polar_threshold = self.options.polar_threshold
        is_polar = np.amax(lat_range) >= polar_threshold
        self.log.info("- hemisphere: %s" % l1b.info.region_name)
        self.log.info("- above polar threshold: %s" % str(is_polar))
#        self.log.debug("- l1b.info.lat_min: %.7g" % l1b.info.lat_min)
#        self.log.debug("- l1b.info.lat_max: %.7g" % l1b.info.lat_max)
        # 2) test if there is any ocean data at all
        has_ocean = l1b.info.open_ocean_percent > 0
        self.log.info("- has ocean data: %s" % str(has_ocean))
#        self.log.debug("- l1b.info.open_ocean_percent: %.7g" % (
#            l1b.info.open_ocean_percent))
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

    @property
    def has_empty_file_list(self):
        return len(self._l1b_file_list) == 0


class L1bPreProcJobSettings(DefaultLoggingClass):

    def __init__(self, config):

        super(L1bPreProcJobSettings, self).__init__(self.__class__.__name__)

        # Save pointer to pysiral configuration
        self.pysiral_config = config

        # Initialize the time range and set to monthly per default
        self.time_range = TimeRangeRequest()

        # Error Status
        self.error = ErrorStatus()

        # Initialize job parameter
        self.args = None
        self.iterations = []

    def parse_command_line_arguments(self):

        # use python module argparse to parse the command line arguments
        # (first validation of required options and data types)
        self.args = self.parser.parse_args()

    def generate_preprocessor_iterations(self):

        # The input data is organized in folder per month, therefore
        # the processing period is set accordingly
        self.time_range.set_period("monthly")
        self.log.info("Pre-processor Period is monthly")
        self.time_range.set_exclude_month(self.args.exclude_month)
        self.log.info("Excluding month: %s" % str(self.args.exclude_month))
        self.iterations = self.time_range.get_iterations()
        self.log.info("Number of iterations: %g" % len(self.iterations))

    def process_requested_time_range(self):

        config = self.pysiral_config
        mission_info = config.get_mission_info(self.args.mission_id)

        # Set and validate the time range
        start_date, stop_date = self.args.start_date, self.args.stop_date

        self.time_range.set_range(start_date, stop_date)

        # Check if any errors in definitions
        self.time_range.error.raise_on_error()

        self.log.info("Requested time range is: %s" % self.time_range.label)

        # Clip time range to mission time range
        is_clipped = self.time_range.clip_to_range(
            mission_info.data_period.start, mission_info.data_period.stop)

        if is_clipped:
            self.log.info("Clipped to mission time range: %s till %s" % (
                mission_info.data_period.start, mission_info.data_period.stop))
            # Check if range is still valid
            self.time_range.raise_if_empty()

    def get_mission_preprocessor(self):

        from pysiral.cryosat2.preproc import CryoSat2PreProc
        from pysiral.envisat.preproc import EnvisatPreProc
        from pysiral.ers.preproc import ERSPreProc
        from pysiral.sentinel3.preproc import Sentinel3PreProc

        if self.mission_id == "cryosat2":
            return CryoSat2PreProc
        elif self.mission_id == "envisat":
            return EnvisatPreProc
        elif self.mission_id == "ers2":
            return ERSPreProc
        elif self.mission_id == "sentinel3a":
            return Sentinel3PreProc
        else:
            error_code = self.__class__.__name__+" (1)"
            error_message = "Invalid mission_id: %s" % self.mission_id
            self.error.add_error(error_code, error_message)

    @property
    def parser(self):
        # XXX: Move back to caller

        # Take the command line options from default settings
        # -> see config module for data types, destination variables, etc.
        clargs = DefaultCommandLineArguments()

        # List of command line option required for pre-processor
        # (argname, argtype (see config module), destination, required flag)
        options = [
            ("-mission", "mission", "mission_id", True),
            ("-start", "date", "start_date", True),
            ("-stop", "date", "stop_date", True),
            ("-hemisphere", "hemisphere", "hemisphere", False),
            ("-exclude-month", "exclude-month", "exclude_month", False),
            ("-input-version", "input-version", "input_version", False),
            ("--remove-old", "remove-old", "remove_old_l1bdata", False)]

        # create the parser
        parser = argparse.ArgumentParser()
        for option in options:
            argname, argtype, destination, required = option
            argparse_dict = clargs.get_argparse_dict(
                argtype, destination, required)
            parser.add_argument(argname, **argparse_dict)

        return parser

    @property
    def mission_id(self):
        return self.args.mission_id

    @property
    def hemisphere(self):
        return self.args.hemisphere

    @property
    def input_version(self):
        return self.args.input_version

    @property
    def remove_old(self):
        return self.args.remove_old_l1bdata


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
