
import os
import numpy as np
from operator import attrgetter
from datetime import timedelta

from pysiral import get_cls
from pysiral.clocks import StopWatch
from pysiral.config import ConfigInfo, TimeRangeRequest, get_yaml_config
from pysiral.helper import (ProgressIndicator, get_first_array_index, get_last_array_index, rle)
from pysiral.errorhandler import ErrorStatus
from pysiral.path import validate_directory
from pysiral.logging import DefaultLoggingClass
from pysiral.output import L1bDataNC


def get_preproc(type, input_adapter, output_handler, cfg):
    """
    A function returning the pre-processor class corresponding the type definition
    :param type: type of the pre-processor (`orbit_segment`)
    :param input_adapter: A class that return a L1bData object for a given input product file
    :param output_handler: A class that creates a pysiral l1p product from the merged L1bData object
    :param cfg: a treedict of options for the pre-processor
    :return: Initialized pre-processor class
    """

    # A lookup dictionary for the appropriate class
    preproc_class_lookup_dict = {"custom_orbit_segment": L1PreProcCustomOrbitSegment}

    # Try the get the class
    cls = preproc_class_lookup_dict.get(type, None)

    # Error handling
    if cls is None:
        msg = "Unrecognized Level-1 Pre-Processor class type: %s" % (str(type))
        msg += "\nKnown types:"
        for key in preproc_class_lookup_dict.keys():
            msg += "\n - %s" % key
        error = ErrorStatus(caller_id="Level1PreProcessor")
        error.add_error("invalid-l1preproc-class", msg)
        error.raise_on_error()

    # Return the initialized class
    return cls(input_adapter, output_handler, cfg)


class L1PreProcBase(DefaultLoggingClass):

    def __init__(self, cls_name, input_adapter, output_handler, cfg):

        # Make sure the logger/error handler has the name of the parent class
        super(L1PreProcBase, self).__init__(cls_name)
        self.error = ErrorStatus(caller_id=cls_name)

        # The class that translates a given input file into an L1BData object
        self.input_adapter = input_adapter

        # Output data handler that creates l1p netCDF files from l1 data objects
        self.output_handler = output_handler

        # The configuration for the pre-processor
        self.cfg = cfg

        # The stack of Level-1 objects is a simple list
        self.l1_stack = []

    def process_input_files(self, input_file_list):
        """
        Main entry point for the Level-Preprocessor.
        :param input_file_list: A list full filepath for the pre-processor
        :return: None
        """

        # Validity Check
        n_input_files = len(input_file_list)
        if n_input_files == 0:
            self.log.warning("Passed empty input file list to process_input_files()")
            return

        # Init helpers
        prgs = ProgressIndicator(n_input_files)

        # A class that is passed to the input adapter to check if the pre-processsor wants the
        # content of the current file
        polar_ocean_check = L1PreProcPolarOceanCheck(self.__class__.__name__, self.polar_ocean_props)

        # orbit segments may or may not be connected, therefore the list of input file
        # needs to be processed sequentially.
        for i, input_file in enumerate(input_file_list):

            # Step 1: Read Input
            # Map the entire orbit segment into on Level-1 data object. This is the task
            # of the input adaptor. The input handler gets only the filename and the target
            # region to assess whether it is necessary to parse and transform the file content
            # for the sake of computational efficiency.
            self.log.info("+ Process input file %s" % prgs.get_status_report(i))
            l1 = self.input_adapter.get_l1(input_file, polar_ocean_check)
            if l1 is None:
                self.log.info("- No polar ocean data for curent job -> skip file")
                continue

            # Step 2: Extract and subset
            # The input files may contain unwanted data (low latitude/land segments). It is the job of the
            # L1PReProc children class to return only the relevant segments over polar ocean as a list of l1 objects.
            l1_segments = self.extract_polar_ocean_segments(l1)

            # Step 3: Post-processing
            # Computational expensive post-processing (e.g. computation of waveform shape parameters) can now be
            # executed as the the Level-1 segments are cropped to the minimal length.
            self.l1_post_processing(l1_segments)

            # Step 4: Merge orbit segments
            # Add the list of orbit segments to the l1 data stack and merge those that are connected
            # (e.g. two half orbits connected at the pole) into a single l1 object. Orbit segments that
            # are unconnected from other segments in the stack will be exported to netCDF files.
            self.l1_stack_merge_and_export(l1_segments)

        # Step : Export the last item in the stack
        l1_merged = self.l1_get_merged_stack()
        self.l1_export_to_netcdf(l1_merged)

    def l1_post_processing(self, l1_segments):
        """
        Apply the post-processing procedures defined in the l1p processor definition file.

        :param l1_segments: A list of Level-1 data objects
        :return: None, the l1_segments are changed in place
        """

        # Get the post processing options
        pre_processing_items = self.cfg.get("pre_processing_items", None)
        if pre_processing_items is None:
            self.log.info("No pre processing items defined")
            return

        # Measure time for the different post processors
        timer = StopWatch()

        # Get the list of post-processing items
        for pp_item in pre_processing_items:
            timer.start()
            pp_class = get_cls(pp_item["module_name"], pp_item["class_name"], relaxed=False)
            post_processor = pp_class(**pp_item["options"])
            for l1 in l1_segments:
                post_processor.apply(l1)
            timer.stop()
            msg = "- L1 pre-processing item `%s` applied in %.3f seconds" % (pp_item["label"], timer.get_seconds())
            self.log.info(msg)

    def l1_stack_merge_and_export(self, l1_segments):
        """
        Add the input Level-1 segments to the l1 stack and export the unconnected ones as l1p netCDF products
        :param l1_segments:
        :return: None
        """

        # Loop over all input segments
        for l1 in l1_segments:

            # Test if l1 segment is connected to stack
            is_connected = self.l1_is_connected_to_stack(l1)

            # Case 1: Segment is connected
            # -> Add the l1 segment to the stack and check the next segment.
            if is_connected:
                self.log.info("- L1 segment connected -> add to stack")
                self.l1_stack.append(l1)

            # Case 2: Segment is not connected
            # -> In this case all items in the l1 stack will be merged and the merged l1 object will be
            #    exported to a l1p netCDF product. The current l1 segment that was unconnected to the stack
            #    will become the next stack
            else:
                self.log.info("- L1 segment unconnected -> exporting current stack")
                l1_merged = self.l1_get_merged_stack()
                self.l1_export_to_netcdf(l1_merged)
                self.l1_stack = [l1]

    def l1_is_connected_to_stack(self, l1):
        """
        Check if the start time of file i and the stop time if file i-1 indicate neighbouring orbit segments
        (e.g. due to radar mode change, or two half-orbits
        :param l1:
        :return: Flag if l1 is connected (True of False)
        """

        # Stack is empty (return True -> create a new stack)
        if self.stack_len == 0:
            return True

        # Test if segments are adjacent based on time gap between them
        timedelta = l1.info.start_time - self.last_stack_item.info.stop_time
        threshold = self.cfg.orbit_segment_connectivity.max_connected_segment_timedelta_seconds
        is_connected = timedelta.seconds <= threshold

        return is_connected

    def l1_get_merged_stack(self):
        """
        Concatenates all items in the l1 stack and returns the merged Level-1 data object.
        Note: This operation leaves the state of the Level-1 stack untouched
        :return: Level-1 data object
        """
        l1_merged = self.l1_stack[0]
        for l1 in self.l1_stack[1:]:
            l1_merged.append(l1)
        return l1_merged

    def l1_export_to_netcdf(self, l1):
        """
        Exports the Level-1 object as as l1p netCDF
        :param l1_merged: The Level-1 object to exported
        :return:
        """
        self.output_handler.export_to_netcdf(l1)
        self.log.info("- Written l1p product: %s" % self.output_handler.last_written_file)


    def trim_single_hemisphere_segment_to_polar_region(self, l1):
        """
        Extract polar region of interest from a segment that is either north or south (not global)

        :param l1: Input Level-1 object
        :return: Trimmed Input Level-1 object
        """
        polar_threshold = self.cfg.polar_ocean.polar_latitude_threshold
        is_polar = np.abs(l1.time_orbit.latitude) >= polar_threshold
        polar_subset = np.where(is_polar)[0]
        if len(polar_subset) != l1.n_records:
            l1.trim_to_subset(polar_subset)
        return l1

    def trim_two_hemisphere_segment_to_polar_regions(self, l1):
        """
        Extract polar regions of interest from a segment that is either north, south or both. The method will
        preserve the order of the hemispheres

        :param l1: Input Level-1 object
        :return: List of Trimmed Input Level-1 objects
        """

        polar_threshold = self.cfg.polar_ocean.polar_latitude_threshold
        l1_list = []

        # Loop over the two hemispheres
        for hemisphere in self.cfg.polar_ocean.target_hemisphere:

            if hemisphere == "north":
                is_polar = l1.time_orbit.latitude >= polar_threshold

            elif hemisphere == "south":
                is_polar = l1.time_orbit.latitude <= (-1.0*polar_threshold)

            else:
                msg = "Unknown hemisphere: %s [north|south]" % hemisphere
                self.error.add_error("invalid-hemisphere", msg)
                self.error.raise_on_error()

            # Extract the subset (if applicable)
            polar_subset = np.where(is_polar)[0]
            n_records_subset = len(polar_subset)

            # is true subset -> add subset to output list
            if n_records_subset != l1.n_records and n_records_subset > 0:
                l1_segment = l1.extract_subset(polar_subset)
                l1_list.append(l1_segment)

            # entire segment in polar region -> add full segment to output list
            elif n_records_subset == l1.n_records:
                l1_list.append(l1)

            # no coverage in target hemisphere -> remove segment from list
            else:
                pass

        # Last step: Sort the list to maintain temporal order
        # (only if more than 1 segment)
        if len(l1_list) > 1:
            l1_list = sorted(l1_list, key=attrgetter("tcs"))

        return l1_list

    def filter_small_ocean_segments(self, l1):
        """
        This method sets the surface type flag of very small ocean segments to land. This action should prevent
        large portions of land staying in the l1 segment is a small fjord etc is crossed. It should also filter
        out smaller ocean segments that do not have a realistic chance of freeboard retrieval.

        :param l1: A pysiral.l1bdata.Level1bData instance
        :return: filtered l1 object
        """

        # Minimum size for valid ocean segments
        ocean_mininum_size_nrecords = self.cfg.polar_ocean.ocean_mininum_size_nrecords

        # Get the clusters of ocean parts in the l1 object
        ocean_flag = l1.surface_type.get_by_name("ocean").flag
        land_flag = l1.surface_type.get_by_name("land").flag
        segments_len, segments_start, not_ocean = rle(ocean_flag)

        # Find smaller than threshold ocean segments
        small_cluster_indices = np.where(segments_len < ocean_mininum_size_nrecords)[0]

        # Do not mess with the l1 object if not necessary
        if len(small_cluster_indices == 0):
            return l1

        # Set land flag -> True for small ocean segments
        for small_cluster_index in small_cluster_indices:
            i0 = segments_start[small_cluster_index]
            i1 = i0 + segments_len[small_cluster_index]
            land_flag[i0:i1] = True

        # Update the l1 surface type flag by re-setting the land flag
        l1.surface_type.add_flag(land_flag, "land")

        # All done
        return l1

        # import matplotlib.pyplot as plt
        # import sys
        #
        # print segments_len
        #
        # plt.figure()
        # plt.plot(ocean_flag, alpha=0.5)
        # plt.plot(land_flag, alpha=0.5)
        # plt.show()
        # sys.exit()

    def trim_non_ocean_data(self, l1):
        """
        Remove leading and trailing data that is not if type ocean. 
        :param l1: The input Level-1 objects
        :return: The subsetted Level-1 objects. (Segments with no ocean data are removed from the list)
        """""" """

        ocean = l1.surface_type.get_by_name("ocean")
        first_ocean_index = get_first_array_index(ocean.flag, True)
        last_ocean_index = get_last_array_index(ocean.flag, True)
        if first_ocean_index is None or last_ocean_index is None:
            return None
        n = l1.info.n_records-1
        is_full_ocean = first_ocean_index == 0 and last_ocean_index == n
        if not is_full_ocean:
            ocean_subset = np.arange(first_ocean_index, last_ocean_index+1)
            l1.trim_to_subset(ocean_subset)
        return l1

    def split_at_large_non_ocean_segments(self, l1):
        """
        Identify larger segments that are not ocean (land, land ice) and split the segments if necessary.
        The return value will always be a list of Level-1 object instances, even if no non-ocean data
        segment is present in the input data file
        :param l1: Input Level-1 object
        :return: a list of Level-1 objects.
        """

        # Identify connected non-ocean segments within the orbit
        ocean = l1.surface_type.get_by_name("ocean")
        not_ocean_flag = np.logical_not(ocean.flag)
        segments_len, segments_start, not_ocean = rle(not_ocean_flag)
        landseg_index = np.where(not_ocean)[0]

        # no non-ocean segments, return full segment
        if len(landseg_index) == 0:
            return [l1]

        # Test if non-ocean segments above the size threshold that will require a split of the segment.
        # The motivation behind this step to keep l1p data files as small as possible, while tolerating
        # smaller non-ocean sections
        treshold = self.cfg.polar_ocean.allow_nonocean_segment_nrecords
        large_landsegs_index = np.where(segments_len[landseg_index] > treshold)[0]
        large_landsegs_index = landseg_index[large_landsegs_index]

        # no segment split necessary, return full segment
        if len(large_landsegs_index) == 0:
            return [l1]

        # Split of orbit segment required, generate individual Level-1 segments from the ocean segments
        l1_segments = []
        start_index = 0
        for index in large_landsegs_index:
            stop_index = segments_start[index]
            subset_list = np.arange(start_index, stop_index)
            l1_segments.append(l1.extract_subset(subset_list))
            start_index = segments_start[index+1]

        # Extract the last subset
        last_subset_list = np.arange(start_index, len(ocean.flag))
        l1_segments.append(l1.extract_subset(last_subset_list))

        # Return a list of segments
        return l1_segments

    def split_at_time_discontinuities(self, l1_list):
        """
        Split l1 object(s) at discontinuities of the timestamp value and return the expanded list with l1 segments.

        :param l1_list: [list] a list of l1b_files
        :return: expanded list
        """

        # Prepare input (should always be list)
        seconds_threshold = self.cfg.timestamp_discontinuities.split_at_time_gap_seconds
        dt_threshold = timedelta(seconds=seconds_threshold)

        # Output (list with l1b segments)
        l1_segments = []

        for l1 in l1_list:

            # Get timestamp discontinuities (if any)
            time = l1.time_orbit.timestamp

            # Get start start/stop indices pairs
            segments_start = np.array([0])
            segments_start_indices = np.where(np.ediff1d(time) > dt_threshold)[0]+1
            segments_start = np.append(segments_start, segments_start_indices)

            segments_stop = segments_start[1:]-1
            segments_stop = np.append(segments_stop, len(time)-1)

            # Check if only one segment found
            if len(segments_start) == 1:
                l1_segments.append(l1)
                continue

            # Extract subsets
            segment_indices = zip(segments_start, segments_stop)
            for start_index, stop_index in segment_indices:
                subset_indices = np.arange(start_index, stop_index+1)
                l1_segment = l1.extract_subset(subset_indices)
                l1_segments.append(l1_segment)

        return l1_segments

    @property
    def target_region_def(self):
        if not self.cfg.has_key("polar_ocean"):
            msg = "Missing configuration key `polar_ocean` in Level-1 Pre-Processor Options"
            self.error.add_error("l1preproc-missing-option", msg)
            self.error.raise_on_error()
        return self.cfg.polar_ocean

    @property
    def polar_ocean_props(self):
        if not self.cfg.has_key("polar_ocean"):
            msg = "Missing configuration key `polar_ocean` in Level-1 Pre-Processor Options"
            self.error.add_error("l1preproc-missing-option", msg)
            self.error.raise_on_error()
        return self.cfg.polar_ocean

    @property
    def orbit_segment_connectivity_props(self):
        if not self.cfg.has_key("orbit_segment_connectivity"):
            msg = "Missing configuration key `orbit_segment_connectivity` in Level-1 Pre-Processor Options"
            self.error.add_error("l1preproc-missing-option", msg)
            self.error.raise_on_error()
        return self.cfg.orbit_segment_connectivity

    @property
    def stack_len(self):
        return len(self.l1_stack)

    @property
    def last_stack_item(self):
        return self.l1_stack[-1]


class L1PreProcCustomOrbitSegment(L1PreProcBase):
    """ A Pre-Processor for input files with arbitrary segment lenght (e.g. CryoSat-2) """

    def __init__(self, *args):
        super(L1PreProcCustomOrbitSegment, self).__init__(self.__class__.__name__, *args)
        # Override the logger name of the input adapter for better logging experience
        self.input_adapter.log.name = self.__class__.__name__

    def extract_polar_ocean_segments(self, l1):
        """
        Splits the input Level-1 object into the polar ocean segments (e.g. by trimming land at the edges
        or by splitting into several parts if there are land masses with the orbit segment). The returned
        polar ocean segments should be generally free of data over non-ocean parts of the orbit, except
        for smaller parts within the orbit.

        NOTE: This subclass of the Level-1 Pre-Processor is designed for input data type with arbitrary
              orbit segment length (e.g. data of CryoSat-2 where the orbit segments of the input data
              is controlled by the mode mask changes).

        :param l1: A Level-1 data object
        :return: A list of Level-1 data objects (subsets of polar ocean segments from input l1)
        """

        # Step: Filter small ocean segments
        # NOTE: The objective is to remove any small marine regions (e.g. in fjords) that do not have any
        #       reasonable chance of freeboard/ssh retrieval early on in the pre-processing.
        if self.cfg.polar_ocean.has_key("ocean_mininum_size_nrecords"):
            self.log.info("- filter ocean segments")
            l1 = self.filter_small_ocean_segments(l1)

        # Step: Trim the orbit segment to latitude range for the specific hemisphere
        # NOTE: There is currently no case of an input data set that is not of type half-orbit and that
        #       would have coverage in polar regions of both hemisphere. Therefore `l1_subset` is assumed to
        #       be a single Level-1 object instance and not a list of instances.  This needs to be changed if
        #      `input_file_is_single_hemisphere=False`
        self.log.info("- extracting polar region subset(s)")
        if self.cfg.polar_ocean.input_file_is_single_hemisphere:
            l1_list = [self.trim_single_hemisphere_segment_to_polar_region(l1)]
        else:
            l1_list = self.trim_two_hemisphere_segment_to_polar_regions(l1)

        # Step: Split the l1 segments at time discontinuities.
        # NOTE: This step is optional. It requires the presence of the options branch `timestamp_discontinuities`
        #       in the l1proc config file
        if self.cfg.has_key("timestamp_discontinuities"):
            self.log.info("- split at time discontinuities")
            l1_list = self.split_at_time_discontinuities(l1_list)

        # Step: Trim the non-ocean parts of the subset (e.g. land, land-ice, ...)
        # NOTE: Generally it can be assumed that the l1 object passed to this method contains polar ocean data.
        #       But there tests before only include if there is ocean data and data above the polar latitude
        #       threshold. It can therefore happen that trimming the non-ocean data leaves an empty Level-1 object.
        #       In this case an empty list is returned.
        self.log.info("- trim outer non-ocean regions")
        l1_trimmed_list = []
        for l1 in l1_list:
            l1_trimmed = self.trim_non_ocean_data(l1)
            if l1_trimmed is not None:
                l1_trimmed_list.append(l1_trimmed)

        # Step: Split the remaining subset at non-ocean parts.
        # NOTE: There is no need to split the orbit at small features. See option `allow_nonocean_segment_nrecords`
        #       in the l1p processor definition. But even if there are no segments to split, the output will always
        #       be a list per requirements of the Level-1 pre-processor workflow.
        l1_list = []
        for l1 in l1_trimmed_list:
            l1_splitted_list = self.split_at_large_non_ocean_segments(l1)
            l1_list.extend(l1_splitted_list)

        # All done, return the list of polar ocean segments
        return l1_list


class L1PreProcPolarOceanCheck(DefaultLoggingClass):
    """
    A small helper class that can be passed to input adapter to check whether the l1 segment is
    wanted or not
    """

    def __init__(self, log_name, cfg):
        cls_name = self.__class__.__name__
        super(L1PreProcPolarOceanCheck, self).__init__(log_name)
        self.error = ErrorStatus(caller_id=cls_name)

        # Save Parameter
        self.cfg = cfg

    def has_polar_ocean_segments(self, product_metadata):
        """
        Checks if there are polar oceans segments based on the metadata of a L1 data object
        :param product_metadata: A l1bdata.L1BMetaData object
        :return: Boolean Flag (true: in region of interest, false: not in region of interest)
        """

        # 1 Check: Needs ocean data
        if product_metadata.open_ocean_percent <= 1e-6:
            self.log.info("- No ocean data")
            return False

        # 2. Must be in target hemisphere
        # NOTE: the definition of hemisphere in l1 data is above or below the equator
        hemisphere = product_metadata.hemisphere
        target_hemisphere = self.cfg.get("target_hemisphere", None)
        if not hemisphere == "global" and not hemisphere in target_hemisphere:
            self.log.info("- No data in target hemishere: %s" % "".join( self.cfg.target_hemispheres))
            return False

        # 3. Must be at higher latitude than the polar latitude threshold
        lat_range = np.abs([product_metadata.lat_min, product_metadata.lat_max])
        polar_latitude_threshold = self.cfg.get("polar_latitude_threshold", None)
        if np.amax(lat_range) < polar_latitude_threshold:
            msg = "- No data above polar latitude threshold (min:%.1f, max:%.1f) [req:+/-%.1f]"
            msg = msg % (product_metadata.lat_min, product_metadata.lat_max, polar_latitude_threshold)
            self.log.info(msg)
            return False

        # 4. All tests passed
        return True









class Level1PreProcJobDef(DefaultLoggingClass):
    """ A class that contains the information for the Level-1 pre-processor JOB (not the pre-processor class!) """

    def __init__(self, l1p_settings_id_or_file, tcs, tce, exclude_month=[], hemisphere="global",
                 output_handler_cfg={}):
        """
        The settings for the Level-1 pre-processor job
        :param l1p_settings_id_or_file: An id of an proc/l1 processor config file (filename excluding the .yaml
                                        extension) or an full filepath to a yaml config file
        :param tcs: [int list] Time coverage start (YYYY MM [DD])
        :param tce: [int list] Time coverage end (YYYY MM [DD]) [int list]
        :param exclude_month: [int list] A list of month that will be ignored
        :param hemisphere: [str] The target hemisphere (`north`, `south`, `global`:default).
        :param output_handler_cfg: [dict] An optional dictionary with options of the output handler
                                   (`overwrite_protection`: [True, False], `remove_old`: [True, False])
        """

        super(Level1PreProcJobDef, self).__init__(self.__class__.__name__)
        self.error = ErrorStatus()

        # Get pysiral configuration
        # TODO: Move to global
        self._cfg = ConfigInfo()

        # Store command line options
        self._hemisphere = hemisphere

        # Parse the l1p settings file
        self.set_l1p_processor_def(l1p_settings_id_or_file)

        # Get full requested time range
        self._time_range = TimeRangeRequest(tcs, tce, exclude_month=exclude_month)
        self.log.info("Requested time range is %s" % self.time_range.label)

        # Store the data handler options
        self._output_handler_cfg = output_handler_cfg

        # Measure execution time
        self.stopwatch = StopWatch()

    @classmethod
    def from_args(cls, args):
        """ Init the Processor Definition from the pysiral-l1preproc command line argument object """

        # Optional Keywords
        kwargs = {}
        if args.exclude_month is not None:
            kwargs["exclude_month"] = args.exclude_month
        data_handler_cfg = dict()
        data_handler_cfg["overwrite_protection"] = args.overwrite_protection
        data_handler_cfg["remove_old"] = args.remove_old
        kwargs["output_handler_cfg"] = data_handler_cfg
        kwargs["hemisphere"] = args.hemisphere

        # Return the initialized class
        return cls(args.l1p_settings, args.start_date, args.stop_date, **kwargs)

    def set_l1p_processor_def(self, l1p_settings_id_or_file):
        """ Parse the content of the processor definition file """

        # 1. Resolve the absolute file path
        procdef_file_path = self.get_l1p_proc_def_filename(l1p_settings_id_or_file)

        # 2. Read the content
        self.log.info("Parsing L1P processor definition file: %s" % procdef_file_path)
        self._l1pprocdef = get_yaml_config(procdef_file_path)

        # 3. Expand info (input data lookup directories)
        self._get_local_input_directory()

        # 4. update hemisphere for input adapter
        self._l1pprocdef.level1_preprocessor.options.polar_ocean.target_hemisphere = self.target_hemisphere

    def get_l1p_proc_def_filename(self, l1p_settings_id_or_file):
        """ Query pysiral config to obtain filename for processor definition file """

        # A. Check if already filename
        if os.path.isfile(l1p_settings_id_or_file):
            return l1p_settings_id_or_file

        # B. Not a file, try to resolve filename via pysiral config
        filename = self.pysiral_cfg.get_settings_file("proc", "l1", l1p_settings_id_or_file)
        if filename is None:
            msg = "Invalid Level-1 pre-processor definition filename or id: %s\n" % l1p_settings_id_or_file
            msg = msg + " \nRecognized Level-1 pre-processor definitions ids:\n"
            ids = self.pysiral_cfg.get_setting_ids("proc", "l1")
            for id in ids:
                msg = msg + "    - " + id + "\n"
            self.error.add_error("invalid-l1p-outputdef", msg)
            self.error.raise_on_error()
        return filename

    def _get_local_input_directory(self):
        """ Replace the tag for local machine def with the actual path info """
        input_handler_cfg = self.l1pprocdef.input_handler.options
        local_machine_def_tag = input_handler_cfg.local_machine_def_tag
        primary_input_def = self.pysiral_cfg.local_machine.l1b_repository
        platform, tag = self.l1pprocdef.platform, local_machine_def_tag

        # Get the value
        branch = primary_input_def[platform][tag]

        # Sanity Checks
        expected_branch_name = "root.l1b_repository.%s.%s" % (platform, tag)
        if branch.isDangling():
            msg = "Missing definition in `local_machine_def.yaml`. Expected branch: %s"
            msg = msg % expected_branch_name
            self.error.add_error("local-machine-def-missing-tag", msg)
            self.error.raise_on_error()

        # Validity checks
        # TODO: These checks are probably better located in a separate method?
        for key in ["source", "l1p"]:

            # 1. Branch must have specific keys for input and output
            if not branch.has_key(key):
                msg = "Missing definition in `local_machine_def.yaml`. Expected value: %s.%s"
                msg = msg % (expected_branch_name, key)
                self.error.add_error("local-machine-def-missing-tag", msg)
                self.error.raise_on_error()

            # 2. The value of each branch must be a valid directory
            directory = str(branch[key])
            if not os.path.isdir(directory):
                msg = "Invalid directory in `local_machine_def.yaml`: %s is not a valid directory"
                msg = msg % directory
                self.error.add_error("local-machine-def-invalid-dir", msg)
                self.error.raise_on_error()

        # Update the lookup dir parameter
        self.l1pprocdef.input_handler.options.lookup_dir = branch.source

    @property
    def hemisphere(self):
        return self._hemisphere

    @property
    def target_hemisphere(self):
        values = {"north": ["north"], "south": ["south"], "global": ["north", "south"]}
        return values[self.hemisphere]

    @property
    def pysiral_cfg(self):
        return self._cfg

    @property
    def l1pprocdef(self):
        return self._l1pprocdef

    @property
    def time_range(self):
        return self._time_range

    @property
    def output_handler_cfg(self):
        return self._output_handler_cfg


class Level1POutputHandler(DefaultLoggingClass):
    """
    The output handler for l1p product files
    NOTE: This is not a subclass of OutputHandlerbase due to the special nature of pysiral l1p products
    """

    def __init__(self, cfg):

        cls_name = self.__class__.__name__
        super(Level1POutputHandler, self).__init__(cls_name)
        self.error = ErrorStatus(caller_id=cls_name)
        self.cfg = cfg

        self.pysiral_cfg = ConfigInfo()

        # Init class properties
        self._path = None
        self._filename = None

    def remove_old_if_applicable(self, period):
        self.log.warning("Not implemented: self.remove_old_if_applicable")
        return

    def export_to_netcdf(self, l1):
        """
        Workflow to export a Level-1 object to l1p netCDF product. The workflow includes the generation of the
        output path (if applicable).
        :param l1: The Level-1 object to be exported
        :return: None
        """

        # Get filename and path
        self.set_output_filepath(l1)

        # Check if path exists
        validate_directory(self.path)

        # Export the data object
        ncfile = L1bDataNC()
        ncfile.l1b = l1
        ncfile.output_folder = self.path
        ncfile.filename = self.filename
        ncfile.export()

    def set_output_filepath(self, l1):
        """
        Sets the class properties required for the file export
        :param l1: The Level-1 object
        :return: None
        """

        local_machine_def_tag = self.cfg.get("local_machine_def_tag", None)
        if local_machine_def_tag is None:
            msg = "Missing mandatory option %s in l1p processor definition file -> aborting"
            msg = msg % "root.output_handler.options.local_machine_def_tag"
            msg = msg + "\nOptions: \n"+self.cfg.makeReport()
            self.error.add_error("missing-option", msg)
            self.error.raise_on_error()

        # TODO: This is work in progress
        filename_template = "pysiral-l1p-{platform}-{source}-{hemisphere}-{tcs}-{tce}-{file_version}.nc"
        time_fmt = "%Y%m%dT%H%M%S"
        values = {"platform": l1.info.mission,
                  "source": self.cfg.version.source_file_tag,
                  "hemisphere": l1.info.hemisphere,
                  "tcs": l1.time_orbit.timestamp[0].strftime(time_fmt),
                  "tce": l1.time_orbit.timestamp[-1].strftime(time_fmt),
                  "file_version": self.cfg.version.version_file_tag}
        self._filename = filename_template.format(**values)

        local_repository = self.pysiral_cfg.local_machine.l1b_repository
        export_folder = local_repository[l1.info.mission][local_machine_def_tag].l1p
        yyyy = "%04g" % l1.time_orbit.timestamp[0].year
        mm = "%02g" % l1.time_orbit.timestamp[0].month
        self._path = os.path.join(export_folder, l1.info.hemisphere, yyyy, mm)

    @property
    def path(self):
        return self._path

    @property
    def filename(self):
        return self._filename

    @property
    def last_written_file(self):
        return os.path.join(self.path, self.filename)