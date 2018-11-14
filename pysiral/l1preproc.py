
import os

from pysiral.clocks import StopWatch
from pysiral.config import ConfigInfo, TimeRangeRequest, get_yaml_config
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
    preproc_class_lookup_dict = {"orbit_segment": L1PreProcOrbitSegment}

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


class L1PPreProcBase(DefaultLoggingClass):

    def __init__(self, cls_name, input_adapter, output_handler, cfg):

        # Make sure the logger/error handler has the name of the parent class
        super(L1PPreProcBase, self).__init__(cls_name)
        self.error = ErrorStatus(caller_id=cls_name)

        # The class that translates a given input file into an L1BData object
        self.input_adapter = input_adapter

        # The class that will create a pysiral pysiral l1p product file from the current stack
        self.output_handler = output_handler

        # The configuration for the pre-processor
        self.cfg = cfg


    def process_input_files(self, file_list):
        self.log.warning("Not implemented: self.process_input_files")


class L1PreProcOrbitSegment(L1PPreProcBase):
    """ A Pre-Processor for input files with arbitrary segment lenght (e.g. CryoSat-2) """

    def __init__(self, *args):
        super(L1PreProcOrbitSegment, self).__init__(self.__class__.__name__, *args)


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




