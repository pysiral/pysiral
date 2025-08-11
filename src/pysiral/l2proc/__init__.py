# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 14:04:27 2015

@author: Stefan
"""

import sys
from collections import OrderedDict, deque
from datetime import datetime
from pathlib import Path

from dateperiods import DatePeriod
from loguru import logger

from pysiral import psrlcfg
from pysiral.core.config import get_yaml_config
from pysiral.core.datahandler import DefaultAuxdataClassHandler
from pysiral.core.legacy_classes import ErrorStatus, DefaultLoggingClass
from pysiral.core.output import DefaultLevel2OutputHandler, Level2Output
from pysiral.l1data import L1bdataNCFile
from pysiral.l2data import Level2Data
from pysiral.l2proc.procsteps import Level2ProcessorStepOrder

__all__ = ["Level2Processor", "Level2ProductDefinition", "Level2ProcessorStepOrder", "procsteps"]


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
            auxhandler_id = f"{auxdata_type}_{auxhandler.pyclass.lower()}"

            # Add & register the class to the Level-2 processor instance
            self._auxhandlers[auxhandler_id] = auxhandler
            self._registered_auxdata_handlers.append((auxhandler_id, auxdata_type))

            # Report the auxdata class
            logger.info(f"Processor Settings - {auxhandler_id.upper()} auxdata handler registered")

    def set_procstep_handler(self):
        """
        Add and initialize Level-2 processor items

        :return: None
        """
        logger.info("Init Processor Steps")
        cfg = self.l2def.procsteps
        self.procsteps = Level2ProcessorStepOrder(cfg)
        is_valid = self.procsteps.validate()
        if not is_valid:
            raise IOError()
        logger.info("Processor steps initialized and validated")

    def execute_l2_processor_steps(self, l1b, l2):
        """
        Executes all Level-2 processor steps.

        Note: This is where the algorithm are executed in their order defined
              in the Level-2 processor definition file

        :param l1b:
        :param l2:

        :return: None
        """
        # Loop over all Level-2 processing steps.
        # Note: The property `class_instances` return freshly initialized
        #       class instances of the respective processor step
        for procstep_class in self.procsteps.class_instances:
            procstep_class.execute(l1b, l2)

    def _report_output_location(self):
        for output_handler in self._output_handler:
            msg = f"Level-2 Output [{str(output_handler.id)}]: {output_handler.basedir}"
            logger.info(msg)

    def _l2_processing_of_orbit_files(self):
        """ Orbit-wise level2 processing """

        # TODO: Evaluate parallelization
        logger.info("Start Orbit Processing")

        n_files = len(self._l1b_files)

        # loop over l1bdata preprocessed orbits
        for i, l1b_file in enumerate(self._l1b_files):

            # Log the current position in the file stack
            logger.info("+ [ %g of %g ] (%.2f%%)" % (i+1, n_files, float(i+1)/float(n_files)*100.))

            # Read the level-1p file
            try:
                l1b = self._read_l1b_file(l1b_file)
            except RuntimeError:
                logger.error(f"Cannot read {Path(l1b_file).name}, ... skipping")
                continue
            source_primary_filename = Path(l1b_file).parts[-1]

            # Initialize the orbit level-2 data container
            # TODO: replace by proper product metadata transfer
            try:
                period = DatePeriod(l1b.info.start_time, l1b.info.stop_time)
            except SystemExit:
                msg = "Computation of data period caused exception"
                logger.warning(f"[invalid-l1b] {msg}")
                continue

            # Init the Level-2 data object
            l2 = Level2Data(l1b.info, l1b.time_orbit, period=period)

            if l2.n_records < 2:
                msg = "Single-record L1b file"
                logger.warning(f"[invalid-l1b] {msg}")
                continue

            # Overwrite the timeliness value of the l1p input data
            # (requires settings of --force-l2def-record-type option in pysiral-l2proc)
            if self._l2def.force_l2def_record_type:
                l2.info.timeliness = self._l2def.record_type

            # Get auxiliary data from all registered auxdata handlers
            error_status, error_codes = self.get_auxiliary_data(l1b, l2)
            if True in error_status:
                logger.info("- skip file due to auxdata errors")
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
            # self._add_to_orbit_collection(l2)

    def _read_l1b_file(self, l1b_file):
        """ Read a L1b data file (l1bdata netCDF) """
        filename = Path(l1b_file).name
        logger.info(f"- Parsing l1bdata file: {filename}")
        l1b = L1bdataNCFile(l1b_file)
        l1b.parse()
        l1b.info.subset_region_name = self.l2def.hemisphere
        return l1b

    def get_auxiliary_data(self, l1p: 'L1bdataNCFile', l2: 'Level2Data'):
        """ Transfer along-track data from all registered auxdata handler to the l2 data object """

        # Loop over all auxiliary data types. Each type must have:
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
                logger.info(f"- {auxdata_type.upper()} auxdata handler message: {msg}")

            # Check for errors
            if auxclass.error.status:
                error_messages = auxclass.error.get_all_messages()
                for error_message in error_messages:
                    logger.warning(f"! {error_message}")
                    # auxdata handler is persistent, therefore errors status
                    # needs to be reset before next orbit
                    auxclass.error.reset()

            auxdata_error_status.append(auxclass.error.status)
            auxdata_error_codes.extend(auxclass.error.codes)

            logger.info(f"- {auxdata_type.upper()} auxdata handler completed")

        # Return error status list
        return auxdata_error_status, auxdata_error_codes

    def _create_l2_outputs(self, l2):
        for output_handler in self._output_handler:
            output = Level2Output(l2, output_handler)
            logger.info(f"- Write {output_handler.id} data file: {output.export_filename}")

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
            value = f"{self.product_line}/{self.record_type}/{self.platform}/{self.version}/{self.hemisphere_code}"

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
