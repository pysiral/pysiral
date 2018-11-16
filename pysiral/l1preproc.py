
import os
import numpy as np

from pysiral.clocks import StopWatch
from pysiral.config import ConfigInfo, TimeRangeRequest, get_yaml_config
from pysiral.helper import ProgressIndicator
from pysiral.errorhandler import ErrorStatus
from pysiral.logging import DefaultLoggingClass
from pysiral.output import Level1POutput, OutputHandlerBase


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

        # Init
        self.l1_stack = L1OrbitSegmentStack()

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
        polar_ocean_check = L1PreProcPolarOceanCheck(self.polar_ocean_props)

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
            # The input files may contain unwanted data (low latitude/land segments). It is the
            # job of the L1PReProc children class to return only the relevant segments over
            # polar ocean as a list of l1 objects.
            l1_segments = self.get_polar_ocean_segments(l1, **self.polar_ocean_props)

            # Step 3: Merge orbit segments
            # Add the list of orbit segments to the l1 data stack and merge those that are connected
            # (e.g. two half orbits connected at the pole) into a single l1 object.
            self.l1_stack.append_and_merge(l1_segments, **self.orbit_segment_connectivity_props)

            # Step 4: Post-Processing and Export
            # Orbit segments that are known to be unconnected to other part of the stack can be exported
            # to the l1p product after post-processing. The last segment in the stack is always assumed
            # to be connected as its state is only known when parsing the next input file. An empty list
            # indicates that all orbit segments in the stack are connected
            l1_separate_segments = self.l1_stack.extract_separate_segments()
            if len(l1_separate_segments) == 0:
                continue
            self.post_processing_and_export(l1_separate_segments)

        # Step 5: Export the final orbit segment
        # Per definition the last segment
        self.l1_stack.extract_all_segments()
        self.post_processing_and_export(l1_separate_segments)


    def post_processing_and_export(self, l1_separate_segments):
        pass

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


class L1OrbitSegmentStack(DefaultLoggingClass):
    """ A small helper class that handles the operations for the stacking of l1 orbit segments """

    def __init__(self):
        cls_name = self.__class__.__name__
        super(L1OrbitSegmentStack, self).__init__(cls_name)
        self.error = ErrorStatus(caller_id=cls_name)

        # Create an empty stack upon initialization
        self.clear_stack()


    def clear_stack(self):
        # The stack is a simple list, each entry is a tuple of the connected_flag and the l1 segment
        self._stack = []

    def append_and_merge(self, l1_segments, **cfg):
        pass

    def extract_separate_segments(self):
        pass

    def extract_all_segments(self):
        pass


class L1PreProcCustomOrbitSegment(L1PreProcBase):
    """ A Pre-Processor for input files with arbitrary segment lenght (e.g. CryoSat-2) """

    def __init__(self, *args):
        super(L1PreProcCustomOrbitSegment, self).__init__(self.__class__.__name__, *args)


class L1PreProcPolarOceanCheck(DefaultLoggingClass):
    """
    A small helper class that can be passed to input adapter to check whether the l1 segment is
    wanted or not
    """

    def __init__(self, cfg):
        cls_name = self.__class__.__name__
        super(L1PreProcPolarOceanCheck, self).__init__(cls_name)
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
            self.log.info("- No data above polar ocean latitude threshold (%.1f)" % polar_latitude_threshold)
            return False

        # 4. All tests passed
        return True









class Level1PreProcJobDef(DefaultLoggingClass):
    """ A class that contains the information for the Level-1 pre-processor JOB (not the pre-processor class!) """

    def __init__(self, l1p_settings_id_or_file, tcs, tce, exclude_month=[], output_handler_cfg={}):

        super(Level1PreProcJobDef, self).__init__(self.__class__.__name__)
        self.error = ErrorStatus()

        # Get pysiral configuration
        # TODO: Move to global
        self._cfg = ConfigInfo()

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
        data_handler_cfg = {}
        data_handler_cfg["overwrite_protection"] = args.overwrite_protection
        data_handler_cfg["remove_old"] = args.remove_old
        kwargs["output_handler_cfg"] = data_handler_cfg

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
        try:
            lookup_dir = primary_input_def[platform][tag]
        except KeyError:
            msg = "Missing local machine definition for: root.l1b_repository.%s.%s" % ((platform, tag))
            self.error.add_error("[local-machine-def-missing-tag", msg)
            self.error.raise_on_error()
        self.l1pprocdef.input_handler.options.lookup_dir = lookup_dir

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

    def remove_old_if_applicable(self, period):
        self.log.warning("Not implemented: self.remove_old_if_applicable")
        return




