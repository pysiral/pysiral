# -*- coding: utf-8 -*-

"""

"""

import sys
import copy
import numpy as np
from loguru import logger
from attrdict import AttrDict
from pathlib import Path
from operator import attrgetter
from datetime import timedelta
from typing import Union, List, TypeVar, Dict, overload

from dateperiods import DatePeriod, PeriodIterator

from pysiral import psrlcfg
from pysiral.l1bdata import Level1bData, L1bMetaData
from pysiral.l1preproc.procitems import L1PProcItemDef
from pysiral.clocks import StopWatch
from pysiral.config import get_yaml_config
from pysiral.helper import (ProgressIndicator, get_first_array_index, get_last_array_index, rle)
from pysiral.errorhandler import ErrorStatus
from pysiral.core import DefaultLoggingClass
from pysiral.output import L1bDataNC


class Level1PInputHandlerBase(DefaultLoggingClass):
    """
    Base class (mostly for type checking).
    """

    def __init__(self,
                 cfg: AttrDict,
                 raise_on_error: bool = False,
                 cls_name: str = None
                 ) -> None:
        """
        Base class for all input handlers implemented in the `l1_adapter` package
        of the individual altimeter platform packages. Not to be called directly.

        :param cfg: Config dictionary
        :param raise_on_error: Flag how to handle errors
        :param cls_name: children class name (for error status)
        """
        super(Level1PInputHandlerBase, self).__init__(cls_name)
        self.error = ErrorStatus(caller_id=cls_name)
        self.cfg = cfg
        self.raise_on_error = raise_on_error

    def get_l1(self,
               filepath: Path,
               **kwargs,
               ) -> Union[None, Level1bData]:
        """
        This method needs to be overwritten by inherting class.

        :param filepath:
        :return:
        """
        raise NotImplementedError("Level1PInputHandlerBase is to be inherited only")


class Level1POutputHandler(DefaultLoggingClass):
    """
    The output handler for l1p product files
    NOTE: This is not a subclass of OutputHandlerbase due to the special nature of pysiral l1p products
    """

    def __init__(self, cfg: AttrDict) -> None:
        cls_name = self.__class__.__name__
        super(Level1POutputHandler, self).__init__(cls_name)
        self.error = ErrorStatus(caller_id=cls_name)
        self.cfg = cfg

        self.pysiral_cfg = psrlcfg

        # Init class properties
        self._path = None
        self._filename = None

    @staticmethod
    def remove_old_if_applicable(perid: DatePeriod) -> None:
        logger.warning("Not implemented: self.remove_old_if_applicable")
        return

    def export_to_netcdf(self, l1: "Level1bData") -> None:
        """
        Workflow to export a Level-1 object to l1p netCDF product. The workflow includes the generation of the
        output path (if applicable).

        :param l1: The Level-1 object to be exported

        :return: None
        """

        # Get filename and path
        self.set_output_filepath(l1)

        # Check if path exists
        Path(self.path).mkdir(exist_ok=True, parents=True)

        # Export the data object
        ncfile = L1bDataNC()
        ncfile.l1b = l1
        ncfile.output_folder = self.path
        ncfile.filename = self.filename
        ncfile.export()

    def set_output_filepath(self, l1: "Level1bData") -> None:
        """
        Sets the class properties required for the file export

        :param l1: The Level-1 object

        :return: None
        """

        local_machine_def_tag = self.cfg.get("local_machine_def_tag", None)
        if local_machine_def_tag is None:
            msg = "Missing mandatory option %s in l1p processor definition file -> aborting"
            msg %= "root.output_handler.options.local_machine_def_tag"
            msg = msg + "\nOptions: \n" + self.cfg.makeReport()
            self.error.add_error("missing-option", msg)
            self.error.raise_on_error()

        # TODO: Move to config file
        filename_template = "pysiral-l1p-{platform}-{source}-{timeliness}-{hemisphere}-{tcs}-{tce}-{file_version}.nc"
        time_fmt = "%Y%m%dT%H%M%S"
        values = {"platform": l1.info.mission,
                  "source": self.cfg.version.source_file_tag,
                  "timeliness": l1.info.timeliness,
                  "hemisphere": l1.info.hemisphere,
                  "tcs": l1.time_orbit.timestamp[0].strftime(time_fmt),
                  "tce": l1.time_orbit.timestamp[-1].strftime(time_fmt),
                  "file_version": self.cfg.version.version_file_tag}
        self._filename = filename_template.format(**values)

        local_repository = self.pysiral_cfg.local_machine.l1b_repository
        export_folder = Path(local_repository[l1.info.mission][local_machine_def_tag]["l1p"])
        yyyy = "%04g" % l1.time_orbit.timestamp[0].year
        mm = "%02g" % l1.time_orbit.timestamp[0].month
        self._path = export_folder / self.cfg.version["version_file_tag"] / l1.info.hemisphere / yyyy / mm

    @property
    def path(self) -> Path:
        return Path(self._path)

    @property
    def filename(self) -> str:
        return self._filename

    @property
    def last_written_file(self) -> Path:
        return self.path / self.filename


L1PInputCLS = TypeVar("L1PInputCLS", bound=Level1PInputHandlerBase)


class L1PreProcBase(DefaultLoggingClass):

    def __init__(self,
                 cls_name: str,
                 input_adapter: L1PInputCLS,
                 output_handler: Level1POutputHandler,
                 cfg: AttrDict
                 ) -> None:

        # Make sure the logger/error handler has the name of the parent class
        super(L1PreProcBase, self).__init__(cls_name)
        self.error = ErrorStatus(caller_id=cls_name)

        # The class that translates a given input file into an L1BData object
        self.input_adapter = input_adapter

        # Output data handler that creates l1p netCDF files from l1 data objects
        self.output_handler = output_handler

        # The configuration for the pre-processor
        self.cfg = cfg

        # Initialize the L1 processor items
        self.processor_item_dict = {}
        self._init_processor_items()

        # The stack of Level-1 objects is a simple list
        self.l1_stack = []

    def _init_processor_items(self) -> None:
        """
        Popuplate the processor item dictionary with initialized processor item classes
        for the different stages of the Level-1 pre-processor.

        This method will evaluate the configuration dictionary passed to this class,
        retrieve the corresponding classes, initialize them and store in a dictionary.

        Dictionary entries are the names of the processing stages and the content
        is a list of classes (cls, label) per processing stage.

        :return:
        """

        for processing_item_definition_dict in self.cfg.processing_items:
            def_ = L1PProcItemDef.from_l1procdef_dict(processing_item_definition_dict)
            cls = def_.get_initialized_processing_item_instance()
            if def_.stage in self.processor_item_dict:
                self.processor_item_dict[def_.stage].append((cls, def_.label,))
            else:
                self.processor_item_dict[def_.stage] = [(cls, def_.label)]

    def process_input_files(self, input_file_list: List[Union[Path, str]]) -> None:
        """
        Main entry point for the Level-1 pre-processor and start of the main
        loop of this class:

        1. Sequentially reads source file from input list of source files
        2. Applies processing items (stage: `post_source_file`)
        3. Extract polar ocean segments (can be multiple per source file)
        4. Applies processing items (stage: `post_polar_ocean_segment_extraction`)
        5. Creates a stack of spatially connecting segments
           a. if this includes all segments of the source file the loop
              progresses to the next source file
           b. if the stack contains a spatially unconnected segments
                1. merge all connected segments to a single l1 object
                2. Applies processing items (stage: `post_merge`)
                3. Export to l1p netCDF
        6. Export the last l1p segment at the end of the loop.

        :param input_file_list: A list full filepath for the pre-processor

        :return: None
        """

        # Validity Check
        n_input_files = len(input_file_list)
        if n_input_files == 0:
            logger.warning("Passed empty input file list to process_input_files()")
            return

        # Init helpers
        prgs = ProgressIndicator(n_input_files)

        # A class that is passed to the input adapter to check if the pre-processor wants the
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
            logger.info(f"+ Process input file {prgs.get_status_report(i)}")
            l1 = self.input_adapter.get_l1(input_file, polar_ocean_check=polar_ocean_check)
            if l1 is None:
                logger.info("- No polar ocean data for curent job -> skip file")
                continue

            # Step 2: Apply processor items on source data
            # for the sake of computational efficiency.
            self.l1_apply_processor_items(l1, "post_source_file")

            # Step 3: Extract and subset
            # The input files may contain unwanted data (low latitude/land segments).
            # It is the job of the L1PReProc children class to return only the relevant
            # segments over polar ocean as a list of l1 objects.
            l1_segments = self.extract_polar_ocean_segments(l1)

            self.l1_apply_processor_items(l1_segments, "post_ocean_segment_extraction")

            # Step 5: Merge orbit segments
            # Add the list of orbit segments to the l1 data stack and merge those that
            # are connected (e.g. two half orbits connected at the pole) into a single l1
            # object. Orbit segments that  are unconnected from other segments in the stack
            # will be exported to netCDF files.
            l1_merged = self.l1_stack_merge(l1_segments)
            if l1_merged is None:
                continue

            # Step 4: Processor items post
            # Computational expensive post-processing (e.g. computation of waveform shape parameters) can now be
            # executed as the Level-1 segments are cropped to the minimal length.
            self.l1_apply_processor_items(l1_merged, "post_merge")

            # Step 5: Export
            self.l1_export_to_netcdf(l1_merged)

        # Step : Export the last item in the stack
        l1_merged = self.l1_get_merged_stack()
        self.l1_export_to_netcdf(l1_merged)

    def extract_polar_ocean_segments(self, l1: "Level1bData") -> List["Level1bData"]:
        """
        Needs to be implemented by child classes, because the exact algorithm
        depends on th source data.

        :param l1:

        :return: None

        :raises: NotImplementedError
        """
        raise NotImplementedError("")

    def l1_apply_processor_items(self,
                                 l1: Union["Level1bData", List["Level1bData"]],
                                 stage_name: str
                                 ) -> None:
        """
        Apply the processor items defined in the l1 processor configuration file
        to either a l1 data object or a list of l2 data objects at a defined
        stage of the procesor.

        This method is a wrapper that deals with multiple input types.
        The functionality is implemented in `_l1_apply_proc_item`

        :param l1: Level-1 data object or list of Level-1 data objects
        :param stage_name: Name of the processing stage. Valid options are
            (`post_source_file`, `post_ocean_segment_extraction`, `post_merge`)

        :return: None, the l1_segments are changed in place
        """

        # Check if there is anything to do first
        if stage_name not in self.processor_item_dict:
            return
        (
            [self._l1_apply_proc_item(l1_item, stage_name) for l1_item in l1] if isinstance(l1, list)
            else self._l1_apply_proc_item(l1, stage_name)
        )

    def _l1_apply_proc_item(self,
                            l1_item: "Level1bData",
                            stage_name: str
                            ) -> None:
        """
        Apply a list of

        :param l1_item:
        :param stage_name:

        :return:
        """
        for procitem, label in self.processor_item_dict.get(stage_name, []):
            timer = StopWatch()
            timer.start()
            procitem.apply(l1_item)
            timer.stop()
            msg = f"- L1 processing item {stage_name}:{label} applied in {timer.get_seconds():.3f} seconds"
            logger.info(msg)
        # breakpoint()
        # # Get the post-processing options
        # pre_processing_items = self.cfg.get("pre_processing_items", None)
        # if pre_processing_items is None:
        #     logger.info("No pre processing items defined")
        #     return
        #
        # # Measure time for the different post processors
        # timer = StopWatch()
        #
        # # Get the list of post-processing items
        # for pp_item in pre_processing_items:
        #     timer.start()
        #     pp_class = get_cls(pp_item["module_name"], pp_item["class_name"], relaxed=False)
        #     post_processor = pp_class(**pp_item["options"])
        #     for l1 in l1_segments:
        #         post_processor.apply(l1)
        #     timer.stop()
        #     msg = "- L1 pre-processing item `%s` applied in %.3f seconds" % (pp_item["label"], timer.get_seconds())
        #     logger.info(msg)

    def l1_stack_merge(self, l1_segments: List["Level1bData"]) -> Union[None, Level1bData]:
        """
        Add the input Level-1 segments to the l1 stack and

        :param l1_segments: List of L1 data sets

        :return: None or merged stack if
        """

        # Loop over all input segments
        for l1 in l1_segments:

            # Test if l1 segment is connected to stack
            is_connected = self.l1_is_connected_to_stack(l1)

            # Case 1: Segment is connected
            # -> Add the l1 segment to the stack and check the next segment.
            if is_connected:
                logger.info("- L1 segment connected -> add to stack")
                self.l1_stack.append(l1)
                return None

            # Case 2: Segment is not connected
            # -> In this case all items in the l1 stack will be merged and the merged l1 object will be
            #    exported to a l1p netCDF product. The current l1 segment that was unconnected to the stack
            #    will become the next stack
            else:
                logger.info("- L1 segment unconnected -> exporting current stack")
                l1_merged = self.l1_get_merged_stack()
                self.l1_stack = [l1]
                return l1_merged

    def l1_is_connected_to_stack(self, l1: "Level1bData") -> bool:
        """
        Check if the start time of file i and the stop time if file i-1 indicate neighbouring orbit segments.
        (e.g. due to radar mode change, or two half-orbits)

        :param l1:

        :return: Flag if l1 is connected (True of False)
        """

        # Stack is empty (return True -> create a new stack)
        if self.stack_len == 0:
            return True

        # Test if segments are adjacent based on time gap between them
        tdelta = l1.info.start_time - self.last_stack_item.info.stop_time
        threshold = self.cfg.orbit_segment_connectivity.max_connected_segment_timedelta_seconds
        return tdelta.seconds <= threshold

    def l1_get_merged_stack(self) -> "Level1bData":
        """
        Concatenates all items in the l1 stack and returns the merged Level-1 data object.
        Note: This operation leaves the state of the Level-1 stack untouched

        :return: Level-1 data object
        """
        l1_merged = copy.deepcopy(self.l1_stack[0])
        for l1 in self.l1_stack[1:]:
            l1_merged.append(l1)
        return l1_merged

    def l1_export_to_netcdf(self, l1: "Level1bData") -> None:
        """
        Exports the Level-1 object as l1p netCDF

        :param l1: The Level-1 object to exported

        :return:
        """

        minimum_n_records = self.cfg.get("export_minimum_n_records", 0)
        if l1.n_records >= minimum_n_records:
            self.output_handler.export_to_netcdf(l1)
            logger.info(f"- Written l1p product: {self.output_handler.last_written_file}")
        else:
            logger.info("- Orbit segment below minimum size (%g), skipping" % l1.n_records)

    def trim_single_hemisphere_segment_to_polar_region(self, l1: "Level1bData") -> "Level1bData":
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

    def trim_two_hemisphere_segment_to_polar_regions(self, l1: "Level1bData"
                                                     ) -> Union[None, List["Level1bData"]]:
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
                is_polar = l1.time_orbit.latitude <= (-1.0 * polar_threshold)

            else:
                msg = f"Unknown hemisphere: {hemisphere} [north|south]"
                self.error.add_error("invalid-hemisphere", msg)
                self.error.raise_on_error()
                return None

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

        # Last step: Sort the list to maintain temporal order
        # (only if more than 1 segment)
        if len(l1_list) > 1:
            l1_list = sorted(l1_list, key=attrgetter("tcs"))

        return l1_list

    def trim_full_orbit_segment_to_polar_regions(self, l1: "Level1bData") -> Union[None, List["Level1bData"]]:
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

            # Compute full polar subset range
            if hemisphere == "north":
                is_polar = l1.time_orbit.latitude >= polar_threshold

            elif hemisphere == "south":
                is_polar = l1.time_orbit.latitude <= (-1.0 * polar_threshold)

            else:
                msg = f"Unknown hemisphere: {hemisphere} [north|south]"
                self.error.add_error("invalid-hemisphere", msg)
                self.error.raise_on_error()
                return None

            # Step: Extract the polar ocean segment for the given hemisphere
            polar_subset = np.where(is_polar)[0]
            n_records_subset = len(polar_subset)

            # Safety check
            if n_records_subset == 0:
                continue
            l1_segment = l1.extract_subset(polar_subset)

            # Step: Trim non-ocean segments
            l1_segment = self.trim_non_ocean_data(l1_segment)

            # Step: Split the polar subset to its marine regions
            l1_segment_list = self.split_at_large_non_ocean_segments(l1_segment)

            # Step: append the ocean segments
            l1_list.extend(l1_segment_list)

        # Last step: Sort the list to maintain temporal order
        # (only if more than 1 segment)
        if len(l1_list) > 1:
            l1_list = sorted(l1_list, key=attrgetter("tcs"))

        return l1_list

    def filter_small_ocean_segments(self, l1: "Level1bData") -> "Level1bData":
        """
        This method sets the surface type flag of very small ocean segments to land. This action should prevent
        large portions of land staying in the l1 segment is a small fjord et cetera is crossed. It should also filter
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

    @staticmethod
    def trim_non_ocean_data(l1: "Level1bData") -> Union[None, "Level1bData"]:
        """
        Remove leading and trailing data that is not if type ocean.

        :param l1: The input Level-1 objects

        :return: The subsetted Level-1 objects. (Segments with no ocean data are removed from the list)
        """

        ocean = l1.surface_type.get_by_name("ocean")
        first_ocean_index = get_first_array_index(ocean.flag, True)
        last_ocean_index = get_last_array_index(ocean.flag, True)
        if first_ocean_index is None or last_ocean_index is None:
            return None
        n = l1.info.n_records - 1
        is_full_ocean = first_ocean_index == 0 and last_ocean_index == n
        if not is_full_ocean:
            ocean_subset = np.arange(first_ocean_index, last_ocean_index + 1)
            l1.trim_to_subset(ocean_subset)
        return l1

    def split_at_large_non_ocean_segments(self, l1: "Level1bData") -> List["Level1bData"]:
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
            start_index = segments_start[index + 1]

        # Extract the last subset
        last_subset_list = np.arange(start_index, len(ocean.flag))
        l1_segments.append(l1.extract_subset(last_subset_list))

        # Return a list of segments
        return l1_segments

    def split_at_time_discontinuities(self, l1_list: List["Level1bData"]) -> List["Level1bData"]:
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

            # Get start/stop indices pairs
            segments_start = np.array([0])
            segments_start_indices = np.where(np.ediff1d(time) > dt_threshold)[0] + 1
            segments_start = np.append(segments_start, segments_start_indices)

            segments_stop = segments_start[1:] - 1
            segments_stop = np.append(segments_stop, len(time) - 1)

            # Check if only one segment found
            if len(segments_start) == 1:
                l1_segments.append(l1)
                continue

            # Extract subsets
            segment_indices = zip(segments_start, segments_stop)
            for start_index, stop_index in segment_indices:
                subset_indices = np.arange(start_index, stop_index + 1)
                l1_segment = l1.extract_subset(subset_indices)
                l1_segments.append(l1_segment)

        return l1_segments

    @property
    def polar_ocean_props(self) -> Union[Dict, AttrDict]:
        if "polar_ocean" not in self.cfg:
            msg = "Missing configuration key `polar_ocean` in Level-1 Pre-Processor Options"
            self.error.add_error("l1preproc-missing-option", msg)
            self.error.raise_on_error()
        return self.cfg.polar_ocean

    @property
    def orbit_segment_connectivity_props(self) -> Union[Dict, AttrDict]:
        if "orbit_segment_connectivity" not in self.cfg:
            msg = "Missing configuration key `orbit_segment_connectivity` in Level-1 Pre-Processor Options"
            self.error.add_error("l1preproc-missing-option", msg)
            self.error.raise_on_error()
        return self.cfg.orbit_segment_connectivity

    @property
    def stack_len(self) -> int:
        return len(self.l1_stack)

    @property
    def last_stack_item(self) -> "Level1bData":
        return self.l1_stack[-1]


class L1PreProcCustomOrbitSegment(L1PreProcBase):
    """ A Pre-Processor for input files with arbitrary segment lenght (e.g. CryoSat-2) """

    def __init__(self, *args):
        super(L1PreProcCustomOrbitSegment, self).__init__(self.__class__.__name__, *args)

    def extract_polar_ocean_segments(self, l1: "Level1bData") -> List["Level1bData"]:
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
        if "ocean_mininum_size_nrecords" in self.cfg.polar_ocean:
            logger.info("- filter ocean segments")
            l1 = self.filter_small_ocean_segments(l1)

        # Step: Trim the orbit segment to latitude range for the specific hemisphere
        # NOTE: There is currently no case of an input data set that is not of type half-orbit and that
        #       would have coverage in polar regions of both hemisphere. Therefore, `l1_subset` is assumed to
        #       be a single Level-1 object instance and not a list of instances.  This needs to be changed if
        #      `input_file_is_single_hemisphere=False`
        logger.info("- extracting polar region subset(s)")
        if self.cfg.polar_ocean.input_file_is_single_hemisphere:
            l1_list = [self.trim_single_hemisphere_segment_to_polar_region(l1)]
        else:
            l1_list = self.trim_two_hemisphere_segment_to_polar_regions(l1)

        # Step: Split the l1 segments at time discontinuities.
        # NOTE: This step is optional. It requires the presence of the options branch `timestamp_discontinuities`
        #       in the l1proc config file
        if "timestamp_discontinuities" in self.cfg:
            logger.info("- split at time discontinuities")
            l1_list = self.split_at_time_discontinuities(l1_list)

        # Step: Trim the non-ocean parts of the subset (e.g. land, land-ice, ...)
        # NOTE: Generally it can be assumed that the l1 object passed to this method contains polar ocean data.
        #       But there tests before only include if there is ocean data and data above the polar latitude
        #       threshold. It can therefore happen that trimming the non-ocean data leaves an empty Level-1 object.
        #       In this case an empty list is returned.
        logger.info("- trim outer non-ocean regions")
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


class L1PreProcHalfOrbit(L1PreProcBase):
    """ A Pre-Processor for input files with a full orbit around the earth (e.g. ERS-1/2) """

    def __init__(self, *args):
        super(L1PreProcHalfOrbit, self).__init__(self.__class__.__name__, *args)

    def extract_polar_ocean_segments(self, l1: "Level1bData") -> List["Level1bData"]:
        """
        Splits the input Level-1 object into the polar ocean segments (e.g. by trimming land at the edges
        or by splitting into several parts if there are land masses with the orbit segment). The returned
        polar ocean segments should be generally free of data over non-ocean parts of the orbit, except
        for smaller parts within the orbit.

        NOTE: This subclass of the Level-1 Pre-Processor is designed for input data type with coverage
              from pole to pole (e.g. Envisat SGDR)

        :param l1: A Level-1 data object

        :return: A list of Level-1 data objects (subsets of polar ocean segments from input l1)
        """

        # Step: Filter small ocean segments
        # NOTE: The objective is to remove any small marine regions (e.g. in fjords) that do not have any
        #       reasonable chance of freeboard/ssh retrieval early on in the pre-processing.
        if "ocean_mininum_size_nrecords" in self.cfg.polar_ocean:
            logger.info("- filter ocean segments")
            l1 = self.filter_small_ocean_segments(l1)

        # Step: Extract Polar ocean segments from full orbit respecting the selected target hemisphere
        logger.info("- extracting polar region subset(s)")
        l1_list = self.trim_two_hemisphere_segment_to_polar_regions(l1)

        # Step: Split the l1 segments at time discontinuities.
        # NOTE: This step is optional. It requires the presence of the options branch `timestamp_discontinuities`
        #       in the l1proc config file
        if "timestamp_discontinuities" in self.cfg:
            logger.info("- split at time discontinuities")
            l1_list = self.split_at_time_discontinuities(l1_list)

        # Step: Trim the non-ocean parts of the subset (e.g. land, land-ice, ...)
        # NOTE: Generally it can be assumed that the l1 object passed to this method contains polar ocean data.
        #       But there tests before only include if there is ocean data and data above the polar latitude
        #       threshold. It can therefore happen that trimming the non-ocean data leaves an empty Level-1 object.
        #       In this case an empty list is returned.
        logger.info("- trim outer non-ocean regions")
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


class L1PreProcFullOrbit(L1PreProcBase):
    """ A Pre-Processor for input files with a full orbit around the earth (e.g. ERS-1/2) """

    def __init__(self, *args):
        super(L1PreProcFullOrbit, self).__init__(self.__class__.__name__, *args)

    def extract_polar_ocean_segments(self, l1: "Level1bData") -> List["Level1bData"]:
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
        if "ocean_mininum_size_nrecords" in self.cfg.polar_ocean:
            logger.info("- filter ocean segments")
            l1 = self.filter_small_ocean_segments(l1)

        # Step: Extract Polar ocean segments from full orbit respecting the selected target hemisphere
        logger.info("- extracting polar region subset(s)")
        l1_list = self.trim_two_hemisphere_segment_to_polar_regions(l1)

        # Step: Split the l1 segments at time discontinuities.
        # NOTE: This step is optional. It requires the presence of the options branch `timestamp_discontinuities`
        #       in the l1proc config file
        if "timestamp_discontinuities" in self.cfg:
            logger.info("- split at time discontinuities")
            l1_list = self.split_at_time_discontinuities(l1_list)

        # Step: Trim the non-ocean parts of the subset (e.g. land, land-ice, ...)
        # NOTE: Generally it can be assumed that the l1 object passed to this method contains polar ocean data.
        #       But there tests before only include if there is ocean data and data above the polar latitude
        #       threshold. It can therefore happen that trimming the non-ocean data leaves an empty Level-1 object.
        #       In this case an empty list is returned.
        logger.info("- trim outer non-ocean regions")
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

    def __init__(self,
                 log_name: str,
                 cfg: AttrDict
                 ) -> None:
        cls_name = self.__class__.__name__
        super(L1PreProcPolarOceanCheck, self).__init__(log_name)
        self.error = ErrorStatus(caller_id=cls_name)

        # Save Parameter
        self.cfg = cfg

    def has_polar_ocean_segments(self, product_metadata: L1bMetaData) -> bool:
        """
        Checks if there are polar oceans segments based on the metadata of a L1 data object
        :param product_metadata: Metadata container of the l1b data product.
        :return: Boolean Flag (true: in region of interest, false: not in region of interest)
        """

        # 1 Check: Needs ocean data
        if product_metadata.open_ocean_percent <= 1e-6:
            logger.info("- No ocean data")
            return False

        # 2. Must be in target hemisphere
        # NOTE: the definition of hemisphere in l1 data is above or below the equator
        hemisphere = product_metadata.hemisphere
        target_hemisphere = self.cfg.get("target_hemisphere", None)
        if hemisphere != "global" and hemisphere not in target_hemisphere:
            logger.info(f'- No data in target hemishere: {"".join(self.cfg.target_hemispheres)}')
            return False

        # 3. Must be at higher latitude than the polar latitude threshold
        lat_range = np.abs([product_metadata.lat_min, product_metadata.lat_max])
        polar_latitude_threshold = self.cfg.get("polar_latitude_threshold", None)
        if np.amax(lat_range) < polar_latitude_threshold:
            msg = "- No data above polar latitude threshold (min:%.1f, max:%.1f) [req:+/-%.1f]"
            msg %= (product_metadata.lat_min, product_metadata.lat_max, polar_latitude_threshold)
            logger.info(msg)
            return False

        # 4. All tests passed
        return True


class Level1PreProcJobDef(DefaultLoggingClass):
    """
    A class that contains all information necessary for the Level-1 pre-processor.
    """

    def __init__(self,
                 l1p_settings_id_or_file: Union[str, Path],
                 tcs: List[int],
                 tce: List[int],
                 exclude_month: List[int] = None,
                 hemisphere: str = "global",
                 platform: str = None,
                 output_handler_cfg: Union[dict, AttrDict] = None,
                 source_repo_id: str = None):
        """
        The settings for the Level-1 pre-processor job

        :param l1p_settings_id_or_file: An id of a proc/l1 processor config file (filename excluding the .yaml
                                        extension) or a full filepath to a yaml config file
        :param tcs: [int list] Time coverage start (YYYY MM [DD])
        :param tce: [int list] Time coverage end (YYYY MM [DD]) [int list]
        :param exclude_month: [int list] A list of month that will be ignored
        :param hemisphere: [str] The target hemisphere (`north`, `south`, `global`:default).
        :param platform: [str] The target platform (pysiral id). Required if l1p settings files is valid for
                               multiple platforms (e.g. ERS-1/2, ...)
        :param output_handler_cfg: [dict] An optional dictionary with options of the output handler
                                   (`overwrite_protection`: [True, False], `remove_old`: [True, False])
        :param source_repo_id: [str] The tag in local_machine_def.yaml (l1b_repository.<platform>.<source_repo_id>)
                                  -> Overwrites the default source repo in the l1p settings
                                     (input_handler.options.local_machine_def_tag &
                                      output_handler.options.local_machine_def_tag)
        """

        super(Level1PreProcJobDef, self).__init__(self.__class__.__name__)
        self.error = ErrorStatus()

        # Get pysiral configuration
        # TODO: Move to global
        self._cfg = psrlcfg
        self._l1pprocdef = None

        # Store command line options
        self._hemisphere = hemisphere
        self._platform = platform
        self._source_repo_id = source_repo_id
        self._exclude_month = exclude_month

        # Parse the l1p settings file
        self._set_l1p_processor_def(l1p_settings_id_or_file)

        # Get full requested time range
        self._time_range = DatePeriod(tcs, tce)
        logger.info(f"Requested time range is {self.time_range.label}")

        # Store the data handler options
        if output_handler_cfg is None:
            output_handler_cfg = {}
        self._output_handler_cfg = output_handler_cfg

        # Measure execution time
        self.stopwatch = StopWatch()

    @classmethod
    def from_args(cls, args: AttrDict) -> "Level1PreProcJobDef":
        """
        Init the Processor Definition from the pysiral-l1preproc command line argument object

        :param args:
        :return:
        """

        # Optional Keywords
        kwargs = {}
        if args.exclude_month is not None:
            kwargs["exclude_month"] = args.exclude_month
        data_handler_cfg = {"overwrite_protection": args.overwrite_protection, "remove_old": args.remove_old}

        if args.source_repo_id is not None:
            data_handler_cfg["local_machine_def_tag"] = args.source_repo_id
        kwargs["output_handler_cfg"] = data_handler_cfg
        kwargs["hemisphere"] = args.hemisphere
        kwargs["platform"] = args.platform
        kwargs["source_repo_id"] = args.source_repo_id

        # Return the initialized class
        return cls(args.l1p_settings, args.start_date, args.stop_date, **kwargs)

    def _set_l1p_processor_def(self, l1p_settings_id_or_file: Union[str, Path]) -> None:
        """
        Parse the content of the processor definition file

        :param l1p_settings_id_or_file: A pysiral known id or list

        :return:
        """

        # 1. Resolve the absolute file path
        procdef_file_path = self._get_l1p_proc_def_filename(l1p_settings_id_or_file)

        # 2. Read the content
        logger.info(f"Parsing L1P processor definition file: {procdef_file_path}")
        self._l1pprocdef = get_yaml_config(procdef_file_path)
        self._check_if_unambiguous_platform()

        # 3. Expand info (input data lookup directories)
        self._get_local_input_directory()

        # 4. update hemisphere for input adapter
        self._l1pprocdef.level1_preprocessor.options.polar_ocean.target_hemisphere = self.target_hemisphere

    def _get_l1p_proc_def_filename(self, l1p_settings_id_or_file: Union[str, Path]) -> Union[str, Path]:
        """
        Query pysiral config to obtain filename for processor definition file

        :param l1p_settings_id_or_file:

        :return: The full file path
        """

        # A. Check if already filename
        if Path(l1p_settings_id_or_file).is_file():
            return l1p_settings_id_or_file

        # B. Not a file, try to resolve filename via pysiral config
        filename = self.pysiral_cfg.get_settings_file("proc", "l1", l1p_settings_id_or_file)
        if filename is None:
            msg = "Invalid Level-1 pre-processor definition filename or id: %s\n" % l1p_settings_id_or_file
            msg = msg + " \nRecognized Level-1 pre-processor definitions ids:\n"
            ids = self.pysiral_cfg.get_setting_ids("proc", "l1")
            for output_id in ids:
                msg = f'{msg}    - {output_id}' + "\n"
            self.error.add_error("invalid-l1p-outputdef", msg)
            self.error.raise_on_error()
        return filename

    def _get_local_input_directory(self) -> None:
        """
        Replace the tag for local machine def with the actual path info
        """

        input_handler_cfg = self.l1pprocdef.input_handler.options
        local_machine_def_tag = input_handler_cfg.local_machine_def_tag
        primary_input_def = self.pysiral_cfg.local_machine.l1b_repository
        platform, tag = self.platform, local_machine_def_tag

        # Overwrite the tag if specifically supplied
        if self._source_repo_id is not None:
            tag = self._source_repo_id

        # Get the value
        expected_branch_name = f"root.l1b_repository.{platform}.{tag}"
        branch = None
        try:
            branch = AttrDict(primary_input_def[platform][tag])
        except KeyError:
            msg = f"Missing definition {expected_branch_name} in {psrlcfg.local_machine_def_filepath}"
            self.error.add_error("local-machine-def-missing-tag", msg)
            self.error.raise_on_error()

        # Sanity Checks
        # TODO: Obsolete?
        if branch is None:
            msg = "Missing definition in `local_machine_def.yaml`. Expected branch: %s"
            msg %= expected_branch_name
            self.error.add_error("local-machine-def-missing-tag", msg)
            self.error.raise_on_error()

        # Validity checks
        # TODO: These checks are probably better located in a separate method?
        for key in ["source", "l1p"]:

            # 1. Branch must have specific keys for input and output
            if key not in branch:
                msg = "Missing definition in `local_machine_def.yaml`. Expected value: %s.%s"
                msg %= (expected_branch_name, key)
                self.error.add_error("local-machine-def-missing-tag", msg)
                self.error.raise_on_error()

            # 2. The value of each branch must be a valid directory or an
            #    attr (e.g. for different radar modes) with a list of directories
            directory_or_attrdict = branch[key]
            try:
                directories = directory_or_attrdict.values()
            except AttributeError:
                directories = [directory_or_attrdict]

            for directory in directories:
                if not Path(directory).is_dir():
                    logger.warning(f"l1p directory does not exist -> creating: {directory}")
                    Path(directory).mkdir(parents=True)

        # Update the lookup dir parameter
        self.l1pprocdef.input_handler["options"]["lookup_dir"] = branch.source

    def _check_if_unambiguous_platform(self) -> None:
        """
        Checks if the platform is unique, since some l1 processor definitions are valid for a series of
        platforms, such as ERS-1/2, Sentinel-3A/B, etc. The indicator is that the platform tag in the
        l1 preprocessor settings is comma separated list.

        For the location of the source data, it is however necessary that the exact platform is known.
        It must therefore be specified explicitly by the -platform argument

        :raises SysExit:

        """

        settings_is_ambigous = "," in self._l1pprocdef.platform
        platform_is_known = self.platform is not None

        if settings_is_ambigous:

            if not platform_is_known:
                msg = "Error: platform in l1p settings is ambiguous (%s), but no platform has been given (-platform)"
                msg %= self._l1pprocdef.platform
                sys.exit(msg)

            if platform_is_known and self.platform not in str(self._l1pprocdef.platform):
                msg = "Error: platform in l1p settings (%s) and given platform (%s) do not match"
                msg %= (self._l1pprocdef.platform, self.platform)
                sys.exit(msg)

        # If platform in settings is unambiguous, but not provided -> get platform from settings
        if not settings_is_ambigous and not platform_is_known:
            self._platform = self._l1pprocdef.platform
            logger.info(f"- get platform from l1p settings -> {self.platform}")

    @property
    def hemisphere(self) -> str:
        return self._hemisphere

    @property
    def target_hemisphere(self) -> str:
        values = {"north": ["north"], "south": ["south"], "global": ["north", "south"]}
        return values[self.hemisphere]

    @property
    def pysiral_cfg(self) -> AttrDict:
        return self._cfg

    @property
    def l1pprocdef(self) -> AttrDict:
        return self._l1pprocdef

    @property
    def time_range(self) -> DatePeriod:
        return self._time_range

    @property
    def period_segments(self) -> PeriodIterator:
        # TODO: Implement exclude months
        return self._time_range.get_segments("month", crop_to_period=True)

    @property
    def output_handler_cfg(self) -> AttrDict:
        return self._output_handler_cfg

    @property
    def platform(self) -> str:
        return self._platform


# Custom type hints
L1PPROC_CLS_TYPE = TypeVar("L1PPROC_CLS_TYPE", bound=L1PreProcBase)


def get_preproc(preproc_type: str,
                input_adapter: L1PInputCLS,
                output_handler: Level1POutputHandler,
                cfg: AttrDict
                ) -> L1PPROC_CLS_TYPE:
    """
    A function returning the pre-processor class corresponding the type definition.

    :param preproc_type: type of the pre-processor
    :param input_adapter: A class that return a L1bData object for a given input product file
    :param output_handler: A class that creates a pysiral l1p product from the merged L1bData object
    :param cfg: options for the pre-processor
    :return: Initialized pre-processor class
    """

    # A lookup dictionary for the appropriate class
    preproc_class_lookup_dict = {"custom_orbit_segment": L1PreProcCustomOrbitSegment,
                                 "half_orbit": L1PreProcHalfOrbit,
                                 "full_orbit": L1PreProcFullOrbit, }

    # Try the get the class
    cls = preproc_class_lookup_dict.get(preproc_type)

    # Error handling
    if cls is None:
        msg = f"Unrecognized Level-1 Pre-Processor class type: {preproc_type}"
        msg += "\nKnown types:"
        for key in preproc_class_lookup_dict:
            msg += "\n - %s" % key
        error = ErrorStatus(caller_id="Level1PreProcessor")
        error.add_error("invalid-l1preproc-class", msg)
        error.raise_on_error()

    # Return the initialized class
    return cls(input_adapter, output_handler, cfg)