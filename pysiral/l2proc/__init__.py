# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 14:04:27 2015

@author: Stefan
"""

import sys

from loguru import logger
from pathlib import Path
from dateperiods import DatePeriod
from collections import deque, OrderedDict
from datetime import datetime

from pysiral import psrlcfg
from pysiral.l1bdata import L1bdataNCFile
from pysiral._class_template import DefaultLoggingClass
from pysiral.config import get_yaml_config
from pysiral.errorhandler import ErrorStatus, PYSIRAL_ERROR_CODES
from pysiral.datahandler import DefaultAuxdataClassHandler

from pysiral.l2data import Level2Data
from pysiral.l2proc.procsteps import Level2ProcessorStepOrder
from pysiral.output import (Level2Output, DefaultLevel2OutputHandler, get_output_class)


__all__ = ["Level2Processor", "Level2ProductDefinition", "L2ProcessorReport", "procsteps"]


class Level2Processor(DefaultLoggingClass):

    def __init__(self, product_def, auxclass_handler=None):
        """
        Init the Level-2 Processor
        :param product_def:
        :param auxclass_handler:
        """
        super(Level2Processor, self).__init__(self.__class__.__name__)

        # Level-2 Algorithm Definition
        # NOTE: This object should ony be called through the property self.l2def
        self._l2def = product_def

        # Auxiliary Data Handler
        # NOTE: retrieves and initializes the auxiliary data classes
        #       based on the l2 processor definition config file
        if auxclass_handler is None:
            auxclass_handler = DefaultAuxdataClassHandler()
        self._auxclass_handler = auxclass_handler

        # TODO: This may also be delegated to its own class (but works for now)
        # This variable will contain a list with the auxiliary data handlers
        self._registered_auxdata_handlers = []
        self._auxhandlers = {}

        # The processing step is a class that links to the different
        self.procsteps = None

        # Output_handler (can be one or many)
        self._output_handler = product_def.output_handler

        # List of Level-2 (processed) orbit segments
        self._orbit = deque()

        # List of Level-1b input files
        self._l1b_files = []

        # pysiral config
        self._config = psrlcfg

        # Processor Initialization Flag
        self._initialized = False

        # Processor summary report
        self.report = L2ProcessorReport()

        # Initialize the class
        self.initialize_processor()

    def set_l1b_files(self, l1b_files):
        self._l1b_files = l1b_files

    def process_l1b_files(self, l1b_files):
        self.set_l1b_files(l1b_files)
        self.run()

    def run(self):
        """ Run the processor """
        self._l2_processing_of_orbit_files()
        self._clean_up()

    def _clean_up(self):
        """ All procedures that need to be reset after a run """
        self.report.clean_up()

    def initialize_processor(self):
        """ Read required auxiliary data sets """

        # Instance can be reused
        if self._initialized:
            # Empty orbit list (or else orbits will accumulate)
            self._orbit.clear()
            return

        logger.info("Starting Initialization")

        # Initialize the auxiliary data handlers
        self.set_auxdata_handler()

        # Initialize the Level-2 processor step handler
        self.set_procstep_handler()

        # Report on output location
        self._report_output_location()

        # All done
        self._initialized = True
        logger.info("Initialization complete")

    def set_auxdata_handler(self):
        """
        Adds all auxdata types from the l2 config file to the Level-2 processor instance
        :return: None
        """

        # The auxdata definition is a list of dictionaries of form {auxdata_type: auxdata_def}
        for auxdata_dict in self.l2def.auxdata:

            # Extract the information
            auxdata_type = list(auxdata_dict.keys())[0]
            auxdata_def = auxdata_dict[auxdata_type]

            # Get options from l2 processor definition file
            # NOTE: These will intentionally override the default options defined in auxdata_def.yaml
            l2_procdef_opt = auxdata_def.get("options", None)

            # Retrieve the class (errors will be caught in method)
            auxhandler = self._auxclass_handler.get_pyclass(auxdata_type, auxdata_def["name"], l2_procdef_opt)

            # Get a unique id of the auxhandler
            auxhandler_id = "%s_%s" % (auxdata_type, auxhandler.pyclass.lower())

            # Add & register the class to the Level-2 processor instance
            self._auxhandlers[auxhandler_id] = auxhandler
            self._registered_auxdata_handlers.append((auxhandler_id, auxdata_type))

            # Report the auxdata class
            logger.info("Processor Settings - %s auxdata handler registered" % auxhandler_id.upper())

    def set_procstep_handler(self):
        """
        Add the proc
        :return:
        """

        logger.info("Init Processor Steps")
        cfg = self.l2def.procsteps
        self.procsteps = Level2ProcessorStepOrder(cfg)
        self.procsteps.validate()
        logger.info("Processor steps initialized and validated")

    def execute_l2_processor_steps(self, l1b, l2):
        """
        Executes all Level-2 processor steps.

        Note: This is where the algorithm are executed in their order defined
              in the Level-2 processor definition file
        :param l1b:
        :param l2:
        :return:
        """

        # Loop over all Level-2 processing steps.
        # Note: The property `class_instances` return freshly initialized
        #       class instances of the respective processor step
        for procstep_class in self.procsteps.class_instances:

            # Execute the processing step with the mandatory method.
            # Note: Each processing step is always supplied with both l1b and
            #       l2 data object, no matter if actually needed
            class_name = procstep_class.classname
            # logger.debug("Executing processing step {}".format(class_name))
            procstep_class.execute(l1b, l2)

    def _report_output_location(self):
        for output_handler in self._output_handler:
            msg = "Level-2 Output [%s]: %s" % (str(output_handler.id), output_handler.basedir)
            logger.info(msg)

    def _initialize_summary_report(self):
        """
        Only add report parameter that are not time range specific
        (e.g. the processor l2 settings)
        """
        self.report.l2_settings_file = self.l2def.l2_settings_file

    def _l2_processing_of_orbit_files(self):
        """ Orbit-wise level2 processing """

        # TODO: Evaluate parallelization
        logger.info("Start Orbit Processing")

        n_files = len(self._l1b_files)

        # loop over l1bdata preprocessed orbits
        for i, l1b_file in enumerate(self._l1b_files):

            # Log the current position in the file stack
            logger.info("+ [ %g of %g ] (%.2f%%)" % (i+1, n_files, float(i+1)/float(n_files)*100.))

            # Read the the level 1b file (l1bdata netCDF is required)
            l1b = self._read_l1b_file(l1b_file)
            source_primary_filename = Path(l1b_file).parts[-1]

            # Initialize the orbit level-2 data container
            # TODO: replace by proper product metadata transfer
            try:
                period = DatePeriod(l1b.info.start_time, l1b.info.stop_time)
            except SystemExit:
                msg = "Computation of data period caused exception"
                logger.warning("[invalid-l1b]", msg)
                continue

            # Init the Level-2 data object
            l2 = Level2Data(l1b.info, l1b.time_orbit, period=period)

            # Overwrite the timeliness value of the l1p input data
            # (requires settings of --force-l2def-record-type option in pysiral-l2proc)
            if self._l2def.force_l2def_record_type:
                l2.info.timeliness = self._l2def.record_type

            # Get auxiliary data from all registered auxdata handlers
            error_status, error_codes = self.get_auxiliary_data(l1b, l2)
            if True in error_status:
                self._discard_l1b_procedure(error_codes, l1b_file)
                continue

            # Execute all Level-2 processor steps
            self.execute_l2_processor_steps(l1b, l2)

            # Create output files
            l2.set_metadata(auxdata_source_dict=self.l2_auxdata_source_dict,
                            source_primary_filename=source_primary_filename,
                            l2_algorithm_id=self.l2def.label,
                            l2_version_tag=self.l2def.file_version_tag)
            self._create_l2_outputs(l2)

            # Add data to orbit stack
            self._add_to_orbit_collection(l2)

    def _read_l1b_file(self, l1b_file):
        """ Read a L1b data file (l1bdata netCDF) """
        filename = Path(l1b_file).name
        logger.info("- Parsing l1bdata file: %s" % filename)
        l1b = L1bdataNCFile(l1b_file)
        l1b.parse()
        l1b.info.subset_region_name = self.l2def.hemisphere
        return l1b

    def _discard_l1b_procedure(self, error_codes, l1b_file):
        """ Log and report discarded l1b orbit segment """
        logger.info("- skip file")
        for error_code in error_codes:
            self.report.add_orbit_discarded_event(error_code, l1b_file)

    def get_auxiliary_data(self, l1p: 'L1bdataNCFile', l2: 'Level2Data'):
        """ Transfer along-track data from all registered auxdata handler to the l2 data object """

        # Loop over all auxilary data types. Each type must have:
        # a) entry in the l2 processing definition under the root.auxdata
        # b) entry in .pysiral-cfg.auxdata.yaml
        # c) python class in pysiral.auxdata.$auxdata_type$

        auxdata_error_status = []
        auxdata_error_codes = []

        for (auxdata_id, auxdata_type) in self.registered_auxdata_handlers:

            # Get the class
            auxclass = self._auxhandlers[auxdata_id]

            # Pass optional l1p data container
            auxclass.receive_l1p_input(l1p)

            # Transfer variables (or at least attempt to)
            auxclass.add_variables_to_l2(l2)

            # Reporting
            for msg in auxclass.msgs:
                logger.info("- %s auxdata handler message: %s" % (auxdata_type.upper(), msg))

            # Check for errors
            if auxclass.error.status:
                error_messages = auxclass.error.get_all_messages()
                for error_message in error_messages:
                    logger.warning("! "+error_message)
                    # auxdata handler is persistent, therefore errors status
                    # needs to be reset before next orbit
                    auxclass.error.reset()

            auxdata_error_status.append(auxclass.error.status)
            auxdata_error_codes.extend(auxclass.error.codes)

            logger.info("- %s auxdata handler completed" % (auxdata_type.upper()))

        # Return error status list
        return auxdata_error_status, auxdata_error_codes

    def _create_l2_outputs(self, l2):
        for output_handler in self._output_handler:
            output = Level2Output(l2, output_handler)
            logger.info("- Write {} data file: {}".format(output_handler.id, output.export_filename))

    def _add_to_orbit_collection(self, l2):
        self._orbit.append(l2)

    @property
    def orbit(self):
        return self._orbit

    @property
    def has_empty_file_list(self):
        return len(self._l1b_files) == 0

    @property
    def l2def(self):
        return self._l2def

    @property
    def registered_auxdata_handlers(self):
        return list(self._registered_auxdata_handlers)

    @property
    def l2_auxdata_source_dict(self):
        """ A dictionary that contains the descriptions of the auxiliary data sources """
        auxdata_dict = {}
        for auxdata_id, auxdata_type in self.registered_auxdata_handlers:
            try:
                handler = self._auxhandlers[auxdata_id]
                if auxdata_type not in auxdata_dict:
                    auxdata_dict[auxdata_type] = handler.longname
                else:
                    auxdata_dict[auxdata_type] = auxdata_dict[auxdata_type]+", "+handler.longname
            except AttributeError:
                auxdata_dict[auxdata_type] = "unspecified"
        return auxdata_dict


class Level2ProductDefinition(DefaultLoggingClass):
    """ Main configuration class for the Level-2 Processor """

    def __init__(self,
                 run_tag: str,
                 l2_settings_file: str,
                 force_l2def_record_type: bool = False) -> None:

        super(Level2ProductDefinition, self).__init__(self.__class__.__name__)
        self.error = ErrorStatus(self.__class__.__name__)

        # Mandatory parameter
        self._l2_settings_file = l2_settings_file
        self._parse_l2_settings()
        self._run_tag = None
        self.force_l2def_record_type = force_l2def_record_type
        self._set_run_tag(run_tag)

        # Optional parameters (may be set to default values if not specified)
        self._output_handler = []

    def add_output_definition(self, output_def_file, period="default", overwrite_protection=True):

        # Set given or default output handler
        self._output_handler.append(DefaultLevel2OutputHandler(
            output_def=output_def_file, subdirectory=self.run_tag,
            period=period, overwrite_protection=overwrite_protection))

    def _parse_l2_settings(self):
        try:
            self._l2def = get_yaml_config(self._l2_settings_file)
        except Exception as ex:
            self.error.add_error("invalid-l2-settings", str(ex))
            self.error.raise_on_error()

    def _set_run_tag(self, run_tag):
        """
        Set the run tag (will be used for defining the output path)
        :param run_tag:
        :return:
        """
        # Take specified value or construct from metadata of config file if unspecified
        if run_tag is not None:
            value = run_tag
        else:
            value = "{}/{}/{}/{}/{}".format(self.product_line, self.record_type,
                                            self.platform, self.version, self.hemisphere_code)
        self._run_tag = value

    @property
    def run_tag(self):
        return self._run_tag

    @property
    def l2def(self):
        return self._l2def

    @property
    def auxdata(self):
        return self.l2def.auxdata

    @property
    def procsteps(self):
        return self.l2def.procsteps

    @property
    def product_line(self):
        return self.l2def.metadata.product_line

    @property
    def record_type(self):
        return self.l2def.metadata.record_type

    @property
    def platform(self):
        return self.l2def.metadata.platform

    @property
    def version(self):
        return self.l2def.metadata.version

    @property
    def file_version_tag(self):
        return self.l2def.metadata.file_version_tag

    @property
    def label(self):
        return self.l2def.metadata.label

    @property
    def hemisphere(self):
        return self.l2def.metadata.hemisphere

    @property
    def hemisphere_code(self):
        codes = dict(north="nh", south="sh")
        return codes.get(self.hemisphere, "global")

    @property
    def output_handler(self):
        # Revert to default output handler if non was specifically added
        if len(self._output_handler) == 0:
            self.add_output_definition("default")
        return self._output_handler


class L2ProcessorReport(DefaultLoggingClass):

    def __init__(self):

        super(L2ProcessorReport, self).__init__(self.__class__.__name__)

        self.n_files = 0
        self.data_period = None
        self.l2_settings_file = "none"
        self.l1b_repository = "none"

        # Counter for error codes
        # XXX: This is a first quick implementation of error codes
        #      (see pysiral.error_handler modules for more info) and
        #      the dev should make sure to use the correct names. A
        #      more formalized way of reporting errors will be added
        #      in future updates
        self._init_error_counters()

    def add_orbit_discarded_event(self, error_code, l1b_file):
        """ Add the l1b file to the list of files with a certain error code """

        # Only except defined error codes
        try:
            self.error_counter[error_code].append(l1b_file)
        except:
            logger.warning("Unknown error code (%s), ignoring" % error_code)

    def write_to_file(self, output_id, directory):
        """ Write a summary file to the defined export directory """

        # Create a simple filename
        filename = Path(directory) / "pysiral-l2proc-summary.txt"
        logger.info("Exporting summary report: %s" % filename)

        lfmt = "  %-16s : %s\n"
        current_time = str(datetime.now()).split(".")[0]
        with open(str(filename), "w") as fhandle:

            # Write infos on settings, host, os, ....
            fhandle.write("# pysiral Level2Processor Summary\n\n")
            fhandle.write(lfmt % ("created", current_time))

            # Brief statistics of files, errors, warnings
            fhandle.write("\n# Processor Statistics\n\n")
            fhandle.write(lfmt % ("l1b files", str(self.n_files)))
            fhandle.write(lfmt % ("errors", str(self.n_discarded_files)))
            fhandle.write(lfmt % ("warnings", str(self.n_warnings)))

            fhandle.write("\n# Processor & Local Machine Settings\n\n")
            fhandle.write(lfmt % ("pysiral version", psrlcfg.version))
            fhandle.write(lfmt % ("python version", sys.version))
            fhandle.write(lfmt % ("hostname", psrlcfg.hostname))

            # More info on this specific run
            fhandle.write(lfmt % ("data period", self.data_period_str))
            fhandle.write(lfmt % ("Level-2 settings", self.l2_settings_file))
            fhandle.write(lfmt % ("l1b repository", self.l1b_repository))

            # List discarded files and reason (error code & description)
            fhandle.write("\n# Detailed Error Breakdown\n\n")
            msg = "  No %s output generated for %g l1b files due " + \
                  "to following errors:\n"
            fhandle.write(msg % (output_id, self.n_discarded_files))

            for error_code in PYSIRAL_ERROR_CODES.keys():
                n_discarded_files = len(self.error_counter[error_code])
                if n_discarded_files == 0:
                    continue
                error_description = PYSIRAL_ERROR_CODES[error_code]
                msg = "\n  %g file(s): [error_code:%s] %s\n" % (
                    n_discarded_files, error_code, error_description)
                fhandle.write(msg)
                for discarded_file in self.error_counter[error_code]:
                    fn = Path(discarded_file).name
                    fhandle.write("  * %s\n" % fn)

    def clean_up(self):
        """ Remove all non-persistent parameter """
        self.data_period = None
        self.l1b_repository = "none"
        self._init_error_counters()

    def _init_error_counters(self):
        self.error_counter = OrderedDict([])
        for error_code in PYSIRAL_ERROR_CODES.keys():
            self.error_counter[error_code] = []

    @property
    def data_period_str(self):
        try:
            return self.time_range.label
        except:
            return "invalid/mission data period"

    @property
    def n_discarded_files(self):
        num_discarded_files = 0
        for error_code in self.error_counter.keys():
            num_discarded_files += len(self.error_counter[error_code])
        return num_discarded_files

    @property
    def n_warnings(self):
        return 0
