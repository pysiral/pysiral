# -*- coding: utf-8 -*-
"""
Created on Wed Aug 03 17:16:56 2016

@author: shendric

"""

from pysiral.config import (ConfigInfo, TimeRangeRequest)
from pysiral.errorhandler import ErrorStatus
from pysiral.flag import FlagContainer
from pysiral.helper import (get_first_array_index, get_last_array_index, rle)
from pysiral.logging import DefaultLoggingClass
from pysiral.output import (PysiralOutputFilenaming, PysiralOutputFolder,
                            L1bDataNC)
from pysiral.path import filename_from_path

import numpy as np
import os
import glob


class L1bPreProc(DefaultLoggingClass):
    """
    This is a parent class that bundles functionality for
    mission-specific pre-processors that are in the modules
    pysiral.$mission$.preproc

    This class alone is not functional
    """

    def __init__(self, name):

        # Enable logging capability (self.log)
        super(L1bPreProc, self).__init__(name)

        # Error handler
        self.error = ErrorStatus()

        # Job definition ( class L1bPreProcJob)
        self._jobdef = None

        # Mission Options
        self._mdef = None

        # List of l1b input files
        # Needs to be filled by the mission specific classes
        self._l1b_file_list = []

        # pysiral configuration
        self._pysiral_config = ConfigInfo()

    def set_job_definition(self, jobdef):
        """ jobdef needs to be of type L1bPreProcJobSettings """
        self._jobdef = jobdef
        self._get_mission_options(jobdef.mission_id)

    def get_input_files_local_machine_def(self, time_range, version="default"):
        self.log.info("Getting input file list from local_machine_def.yaml")
        self._get_input_files_local_machine_def(time_range, version)
        self.log.info("Total number of files: %g" % len(self._l1b_file_list))

    def get_ocean_segments_from_input_file(self, filename):
        """
        Parse mission specific L1b data and return polar ocean
        segments
        """
        # Redirect to mission-specific pre-processor
        return self._get_l1bdata_ocean_segments(filename)

    def execute(self):
        """ Runs the l1b pre-processing """

        # Test if job definition is present
        self._validate_jobdef()
        if self.error.status:
            return

        # Loops over file in self._l1b_file_list
        # -> have to be set by the mission specific processor
        self.merge_and_export_polar_ocean_subsets()

    def remove_old_l1bdata(self):
        """
        Remove all l1b data files in current l1bdata export directory
        defined by the job definition
        """

        # Test if job definition is present
        self._validate_jobdef()
        if self.error.status:
            return

        if self._jobdef.hemisphere == "global":
            hemisphere_targets = ["north", "south"]
        else:
            hemisphere_targets = [self._jobdef.hemisphere]

        for hemisphere in hemisphere_targets:

            # Get the l1bdata diretory
            export_folder = PysiralOutputFolder(config=self._pysiral_config)
            export_folder.l1bdata_from_list(
                self._jobdef.mission_id,
                self._jobdef.input_version,
                hemisphere,
                self._jobdef.time_range.start_dt.year,
                self._jobdef.time_range.start_dt.month)

            # Get list of netcdf files
            search_pattern = os.path.join(export_folder.path, "*.nc")
            l1bdata_files = glob.glob(search_pattern)

            # Delete files
            self.log.info("Removing %g l1bdata files in %s" % (
                len(l1bdata_files), export_folder.path))
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
            log.info("+ [ Processing %g of %g ]" % (i+1, n))
            log.info("- Input file: %s" % filename_from_path(l1b_file))
            l1b_segments = self.get_ocean_segments_from_input_file(l1b_file)

            # Skip if no relevant data was found
            if l1b_segments is None:
                log.info("- %s polar ocean data: False -> skip file" % (
                     self._jobdef.hemisphere))
                continue
            else:
                log.info("- %s polar ocean data: True" % (
                     self._jobdef.hemisphere))
            log.info("- %g polar ocean data segments" % len(l1b_segments))

            # XXX: Debug
            # debug_stack_orbit_plot(l1bdata_stack, l1b_segments)
            # Loop over segments and check connectivity
            for l1b_segment in l1b_segments:

                if not self.l1b_is_connected_to_stack(
                        l1b_segment, l1bdata_stack):

                    # => break criterium, save existing stack and start over
                    log.info("- segment unconnected, exporting current stack")
                    l1b_merged = self.merge_l1b_stack(l1bdata_stack)
                    self.export_l1b_to_netcdf(l1b_merged)

                    # Reset the l1bdata stack
                    l1bdata_stack = [l1b_segment]
                    continue

                # polar ocean data and connected => add to stack
                l1bdata_stack.append(l1b_segment)

    def trim_single_hemisphere_segment_to_polar_region(self, l1b):
        """
        Extract polar region of interest from a segment that is either
        north or south (not global)
        """
        self.log.info("- trimming to polar ocean subset")
        polar_threshold = self._mdef.polar_threshold
        is_polar = np.abs(l1b.time_orbit.latitude) >= polar_threshold
        polar_subset = np.where(is_polar)[0]
        if len(polar_subset) != l1b.n_records:
            l1b.trim_to_subset(polar_subset)
        return l1b

    def trim_non_ocean_data(self, l1b):
        """ Remove leading and trailing data that is not if type ocean """

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

        return l1b

    def split_at_large_non_ocean_segments(self, l1b):
        """
        Identify larger segments that are not ocean (land, land ice)
        and split the segments if necessary
        """

        # Find connected nonocean segments
        self.log.info("- test for inner larger landmasses")
        ocean = l1b.surface_type.get_by_name("ocean")
        not_ocean_flag = np.logical_not(ocean.flag)
        segments_len, segments_start, not_ocean = rle(not_ocean_flag)
        landseg_index = np.where(not_ocean)[0]
        self.log.info("- number of landmasses: %g" % len(landseg_index))

        # no land segements, return segment
        if len(landseg_index) == 0:
            return [l1b]

        # Test if nonocean segments above threshold
        treshold = self._mdef.max_inner_nonocean_segment_nrecords
        large_landsegs_index = np.where(
            segments_len[landseg_index] > treshold)[0]
        large_landsegs_index = landseg_index[large_landsegs_index]
        self.log.info("- number of large landmasses: %g" % len(
            large_landsegs_index))

        # no large land segments, return segment
        if len(large_landsegs_index) == 0:
            return [l1b]

        # Large land segments exist, split the l1b data object
        l1b_segments = []
        start_index = 0
        for index in large_landsegs_index:
            stop_index = segments_start[index]
            subset_list = np.arange(start_index, stop_index)
            l1b_segments.append(l1b.extract_subset(subset_list))
            start_index = segments_start[index+1]
        # extract the last subset
        last_subset_list = np.arange(start_index, len(ocean.flag))
        l1b_segments.append(l1b.extract_subset(last_subset_list))

        # Return a list of segments
        return l1b_segments

    def extract_polar_segments_from_halforbit(self, l1b):
        """
        specific method for mission which data sets are organized in
        half orbits (e.g. ERS-1/2, Envisat)

        Extract Arctic and Antarctic subsets and return both in
        the correct sequence. Returns always two segments, but entries
        might be None if no data for specific hemisphere is found.

        XXX: Do not use for full orbits
        """

        # The threshold is the same for arctic and antarctic (50°)
        lat_threshold = self._mdef.polar_threshold

        # Target hemisphere: north, south or global
        target_hemisphere = self._jobdef.hemisphere

        # Measurement position
        position = l1b.time_orbit

        log_entry = "- Extracted %s subset (%g records)"

        # Extract Arctic (if target hemisphere, else empty list)
        if target_hemisphere == "north" or target_hemisphere == "global":
            is_arctic = FlagContainer(position.latitude > lat_threshold)
            arctic_subset = l1b.extract_subset(is_arctic.indices)
            self.log.info(log_entry % ("Arctic", is_arctic.num))
        else:
            arctic_subset = None

        # Extract Antarctic (if target hemisphere, else empty list)
        if target_hemisphere == "south" or target_hemisphere == "global":
            is_antactic = FlagContainer(position.latitude < -1.*lat_threshold)
            antarctic_subset = l1b.extract_subset(is_antactic.indices)
            self.log.info(log_entry % ("Antarctic", is_antactic.num))
        else:
            antarctic_subset = None

        # Segments may be empty, test all cases
        arctic_is_valid = arctic_subset is not None
        antarctic_is_valid = antarctic_subset is not None

        # Return both unsorted, Empty segment will be discarded later
        if not arctic_is_valid or not antarctic_is_valid:
            return [arctic_subset, antarctic_subset]

        # Order in sequence depeding on start time
        if arctic_subset.info.start_time < antarctic_subset.info.start_time:
            return [arctic_subset, antarctic_subset]
        else:
            return [antarctic_subset, arctic_subset]

    def merge_l1b_stack(self, l1bdata_stack):
        """ Merge a stack of connected l1b objects """

        # Merge the stack to a single l1b instance
        l1b_merged = l1bdata_stack[0]
        l1bdata_stack.pop(0)
        for orbit_segment in l1bdata_stack:
            l1b_merged.append(orbit_segment)

        return l1b_merged

    def export_l1b_to_netcdf(self, l1b):
        """ Write l1b object to netcdf file """

        # Get and create export folder
        export_folder = PysiralOutputFolder(config=self._pysiral_config)
        export_folder.l1bdata_from_l1b(l1b, version=self._jobdef.input_version)
        export_folder.create()

        # Get export filename
        filenaming = PysiralOutputFilenaming()
        export_filename = filenaming.from_l1b(l1b)

        # Export the data object
        ncfile = L1bDataNC()
        ncfile.l1b = l1b
        ncfile.config = self._pysiral_config
        ncfile.output_folder = export_folder.path
        ncfile.filename = export_filename
        ncfile.export()

        log_entry = "- Exported l1bdata: %s"
        self.log.info(log_entry % (export_filename))

    def region_is_arctic_or_antarctic_ocean(self, l1b):
        """
        Test if a l1b file has data over ocean in either Arctic or Antartic
        """

        if l1b is None:
            return False

        # 1) test if either minimum or maximum latitude is in polar regions
        try:
            lat_range = np.abs([l1b.info.lat_min, l1b.info.lat_max])
        except:
            return False

        polar_threshold = self._mdef.polar_threshold
        is_polar = np.amax(lat_range) >= polar_threshold
        self.log.info("- hemisphere: %s" % l1b.info.region_name)
        self.log.info("- above polar threshold: %s" % str(is_polar))

        # 2) test if there is any ocean data at all
        has_ocean = l1b.info.open_ocean_percent > 0
        # self.log.info("- has ocean data: %s" % str(has_ocean))

        # Return a flag and the ocean flag list for later use
        return is_polar and has_ocean

    def l1b_is_connected_to_stack(self, l1b, l1b_stack):
        """
        Check if the start time of file i and the stop time if file i-1
        indicate neighbouring orbit segments (e.g. due to radar mode change)

        """

        # Stack is empty (return True -> create a new stack)
        if len(l1b_stack) == 0:
            return True

        # Test if segments are adjacent based on time gap between them
        last_l1b_in_stack = l1b_stack[-1]
        timedelta = l1b.info.start_time - last_l1b_in_stack.info.stop_time
        threshold = self._mdef.max_connected_files_timedelta_seconds
        is_connected = timedelta.seconds <= threshold

        return is_connected

    def _get_mission_options(self, mission_id):
        """
        Get mission specific preproc options
        (see config/mission.def.yaml
        """
        self._mdef = self._pysiral_config.mission[mission_id].settings

    def _validate_jobdef(self):
        """ Test if job definition is present """
        if self._jobdef is None:
            error_code = self.__class__.__name__+" (01)"
            error_message = "Job definition not set"
            self.error.add_error(error_code, error_message)

    @property
    def has_empty_file_list(self):
        return len(self._l1b_file_list) == 0


class L1bPreProcJob(DefaultLoggingClass):
    """ Container for the definition and handling of a pre-processor job """

    def __init__(self):

        super(L1bPreProcJob, self).__init__(self.__class__.__name__)

        # Save pointer to pysiral configuration
        self.pysiral_config = ConfigInfo()

        # Initialize the time range and set to monthly per default
        self.time_range = TimeRangeRequest()

        # Error Status
        self.error = ErrorStatus()

        # Initialize job parameter
        self.options = L1bPreProcJobOptions()

        # List for iterations (currently only month-wise)
        self.iterations = []

    def generate_preprocessor_iterations(self):
        """ Break the requested time range into monthly iterations """

        # The input data is organized in folder per month, therefore
        # the processing period is set accordingly
        self.time_range.set_period("monthly")
        self.log.info("Pre-processor Period is monthly")
        self.time_range.set_exclude_month(self.options.exclude_month)
        self.log.info("Excluding month: %s" % str(self.options.exclude_month))
        self.iterations = self.time_range.get_iterations()
        self.log.info("Number of iterations: %g" % len(self.iterations))

    def process_requested_time_range(self):
        """
        Verify time range with mission data availability and create
        datetime objects
        """

        config = self.pysiral_config
        mission_info = config.get_mission_info(self.options.mission_id)

        # Set and validate the time range
        start_date, stop_date = self.options.start_date, self.options.stop_date
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
        """ Return the mission specific pre-processor class """

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
            error_code = self.__class__.__name__+" (01)"
            error_message = "Invalid mission_id: %s" % self.mission_id
            self.error.add_error(error_code, error_message)

    @property
    def mission_id(self):
        return self.options.mission_id

    @property
    def hemisphere(self):
        return self.options.hemisphere

    @property
    def input_version(self):
        return self.options.input_version

    @property
    def remove_old(self):
        return self.options.remove_old


class L1bPreProcJobOptions(object):
    """ Simple container for preprocessor options """

    def __init__(self):
        self.mission_id = None
        self.hemisphere = "global"
        self.input_version = "default"
        self.remove_old = False
        self.start_date = None
        self.stop_date = None
        self.exclude_month = []

    def from_dict(self, options_dict):
        for parameter in options_dict.keys():
            if hasattr(self, parameter):
                setattr(self, parameter, options_dict[parameter])


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
