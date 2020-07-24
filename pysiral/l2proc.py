# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 14:04:27 2015

@author: Stefan
"""

from collections import deque, OrderedDict
from datetime import datetime
import numpy as np
import time
import sys
from dateperiods import DatePeriod
from pathlib import Path

from pysiral import get_cls, psrlcfg
from pysiral.config import get_yaml_config
from pysiral.errorhandler import ErrorStatus, PYSIRAL_ERROR_CODES
from pysiral.datahandler import DefaultAuxdataClassHandler
from pysiral.l1bdata import L1bdataNCFile
from pysiral.l2data import Level2Data
from pysiral.logging import DefaultLoggingClass
from pysiral.ssh import get_l2_ssh_class
from pysiral.output import (Level2Output, DefaultLevel2OutputHandler, get_output_class)
from pysiral.surface_type import get_surface_type_class
from pysiral.retracker import get_retracker_class
from pysiral.filter import get_filter
from pysiral.validator import get_validator
from pysiral.frb import get_frb_algorithm
from pysiral.sit import get_sit_algorithm


class Level2Processor(DefaultLoggingClass):

    def __init__(self, product_def, auxclass_handler=None):
        """ Setup of the Level-2 Processor """

        super(Level2Processor, self).__init__(self.__class__.__name__)

        # Error Status Handler
        self.error = ErrorStatus(caller_id=self.__class__.__name__)

        # Level-2 Algorithm Definition
        # NOTE: This object should ony be called through the property self.l2def
        self._l2def = product_def.l2def

        # Auxiliary Data Handler
        # NOTE: retrieves and initializes the auxdata classes based on the l2 processor definition config file
        if auxclass_handler is None:
            auxclass_handler = DefaultAuxdataClassHandler()
        self._auxclass_handler = auxclass_handler

        # TODO: This may also be delegated to its own class (but works for now)
        # This variable will contain a list with the auxiliary data handlers
        self._registered_auxdata_handlers = []
        self._auxhandlers = {}

        # The processing step is a class that links to the different
        self._procsteps = []

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

    def set_l1b_files(self, l1b_files):
        self._l1b_files = l1b_files

    def process_l1b_files(self, l1b_files):
        self.set_l1b_files(l1b_files)
        self.run()

    def run(self):
        """ Run the processor """
        self._l2_processing_of_orbit_files()
        self._l2proc_summary_to_file()
        self._clean_up()

    def _l2proc_summary_to_file(self):
        if "output" not in self.l2def:
            return
        # TODO: This method is currently broken, as there is no output key in l2def
        for output_id, output_def in list(self.l2def.output.items()):
            output = get_output_class(output_def.pyclass)
            output.set_options(**output_def.options)
            output.set_base_export_path(output_def.path)
            time_range = self.report.time_range
            export_folder = output.get_full_export_path(time_range.start)
            self.report.write_to_file(output_id, export_folder)

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

        self.log.info("Starting Initialization")

        # Initialize the auxiliary data handlers
        self.set_auxdata_handler()

        # Initialize the Level-2 processor step handler
        self.set_procstep_handler()

        # Report on output location
        self._report_output_location()

        # All done
        self._initialized = True
        self.log.info("Initialization complete")

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
            self.log.info("Processor Settings - %s auxdata handler registered" % auxhandler_id.upper())

    def set_procstep_handler(self):
        """
        Add the proc
        :return:
        """

        self.log.info("Init Processor Steps")
        cfg = self.l2def.procsteps
        self._procsteps = Level2ProcessorStepOrder(cfg)
        self._procsteps.validate()
        self.log.info("Processor steps initialized and validated")

    def _report_output_location(self):
        for output_handler in self._output_handler:
            msg = "Level-2 Output [%s]: %s" % (str(output_handler.id), output_handler.basedir)
            self.log.info(msg)

    def _initialize_summary_report(self):
        """
        Only add report parameter that are not time range specific
        (e.g. the processor l2 settings)
        """
        self.report.l2_settings_file = self.l2def.l2_settings_file

# %% Level2Processor: orbit processing

    def _l2_processing_of_orbit_files(self):
        """ Orbit-wise level2 processing """

        # TODO: Evaluate parallelization
        self.log.info("Start Orbit Processing")

        n_files = len(self._l1b_files)

        # loop over l1bdata preprocessed orbits
        for i, l1b_file in enumerate(self._l1b_files):

            # Log the current position in the file stack
            self.log.info("+ [ %g of %g ] (%.2f%%)" % (i+1, n_files, float(i+1)/float(n_files)*100.))

            # Read the the level 1b file (l1bdata netCDF is required)
            l1b = self._read_l1b_file(l1b_file)
            source_primary_filename = Path(l1b_file).parts[-1]

            # Apply the geophysical range corrections on the waveform range
            # bins in the l1b data container
            # TODO: move to level1bData class
            self._apply_range_corrections(l1b)

            # Apply a pre-filter of the l1b data (can be none)
            self._apply_l1b_prefilter(l1b)

            # Initialize the orbit level-2 data container
            # TODO: replace by proper product metadata transfer
            try:
                period = DatePeriod(l1b.info.start_time, l1b.info.stop_time)
            except SystemExit:
                msg = "Computation of data period caused exception"
                self.log.warning("[invalid-l1b]", msg)
                continue
            l2 = Level2Data(l1b.info, l1b.time_orbit, period=period)

            # Transfer l1p parameter to the l2 data object (if applicable)
            # NOTE: This is only necessary, if parameters from the l1p files (classifiers) should
            #       be present in the l2i product
            self._transfer_l1p_vars(l1b, l2)

            # Get auxiliary data from all registered auxdata handlers
            error_status, error_codes = self._get_auxiliary_data(l2)
            if True in error_status:
                self._discard_l1b_procedure(error_codes, l1b_file)
                continue

            # Surface type classification (ocean, ice, lead, ...)
            # (ice type classification comes later)
            self._classify_surface_types(l1b, l2)

            # Validate surface type classification
            # yes/no decision on continuing with orbit
            error_status, error_codes = self._validate_surface_types(l2)
            if error_status:
                self._discard_l1b_procedure(error_codes, l1b_file)
                continue

            # Get elevation by retracking of different surface types
            # adds parameter elevation to l2
            error_status, error_codes = self._waveform_retracking(l1b, l2)
            if error_status:
                self._discard_l1b_procedure(error_codes, l1b_file)
                continue

            # Compute the sea surface anomaly (from mss and lead tie points)
            # adds parameter ssh, ssa, afrb to l2
            self._estimate_sea_surface_height(l2)

            # Compute the radar freeboard and its uncertainty
            self._get_altimeter_freeboard(l1b, l2)

            # get radar(-derived) from altimeter freeboard
            self._get_freeboard_from_radar_freeboard(l1b, l2)

            # Apply freeboard filter
            self._apply_freeboard_filter(l2)

            # Convert to thickness
            self._convert_freeboard_to_thickness(l2)

            # Filter thickness
            self._apply_thickness_filter(l2)

            # Post processing
            self._post_processing_items(l2)

            # Create output files
            l2.set_metadata(auxdata_source_dict=self.l2_auxdata_source_dict,
                            source_primary_filename=source_primary_filename,
                            l2_algorithm_id=self.l2def.id,
                            l2_version_tag=self.l2def.version_tag)
            self._create_l2_outputs(l2)

            # Add data to orbit stack
            self._add_to_orbit_collection(l2)

    def _read_l1b_file(self, l1b_file):
        """ Read a L1b data file (l1bdata netCDF) """
        filename = Path(l1b_file).name
        self.log.info("- Parsing l1bdata file: %s" % filename)
        l1b = L1bdataNCFile(l1b_file)
        l1b.parse()
        l1b.info.subset_region_name = self.l2def.hemisphere
        return l1b

    def _discard_l1b_procedure(self, error_codes, l1b_file):
        """ Log and report discarded l1b orbit segment """
        self.log.info("- skip file")
        for error_code in error_codes:
            self.report.add_orbit_discarded_event(error_code, l1b_file)

    # TODO: Marked as obsolete
    def _apply_range_corrections(self, l1b):
        """ Apply the range corrections """
        # XXX: This should be applied to the L2 data not l1b
        for correction in self.l2def.corrections:
            l1b.apply_range_correction(correction)

    # TODO: Marked as obsolete
    def _apply_l1b_prefilter(self, l1b):
        """ Apply filtering of l1b variables """
        # Backward compatibility with older l2 setting files
        if "l1b_pre_filtering" not in self.l2def:
            return
        # Apply filters
        names, filters = self.l2def.l1b_pre_filtering.items()
        for name, filter_def in zip(names, filters):
            self.log.info("- Apply l1b pre-filter: %s" % filter_def.pyclass)
            l1bfilter = get_filter(filter_def.pyclass)
            l1bfilter.set_options(**filter_def.options)
            l1bfilter.apply_filter(l1b)

    # TODO: Marked as obsolete
    def _transfer_l1p_vars(self, l1b, l2):
        """ Transfer variables from l1p to l2 object"""

        # Make this a backward compatible feature (should work without tag in l2 processor definition file)
        if "transfer_from_l1p" not in self.l2def:
            return

        # Don't spam the log
        try:
            verbose = self.l2def.transfer_from_l1p.options.get("verbose")
        except:
            verbose = False

        # Get and loop over data groups
        l1p_items = self.l2def.transfer_from_l1p.items()
        for data_group, varlist in list(l1p_items):

            if data_group == "options":
                continue

            # Get and loop over variables per data group
            l1p_variables = varlist.items()
            for var_name, vardef in list(l1p_variables):

                # Get variable via standard getter method
                # NOTE: Will return None if not found -> create an empty array
                var = l1b.get_parameter_by_name(data_group, var_name)
                if var is None:
                    var = np.full((l2.n_records), np.nan)

                # Add variable to l2 object as auxiliary variable
                l2.set_auxiliary_parameter(vardef["aux_id"], vardef["aux_name"], var, None)

                if verbose:
                    self.log.info("- Transfered l1p variable: %s.%s" % (data_group, var_name))

    def _get_auxiliary_data(self, l2):
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

            # Transfer variables (or at least attempt to)
            auxclass.add_variables_to_l2(l2)

            # Reporting
            for msg in auxclass.msgs:
                self.log.info("- %s auxdata handler message: %s" % (auxdata_type.upper(), msg))

            # Check for errors
            if auxclass.error.status:
                error_messages = auxclass.error.get_all_messages()
                for error_message in error_messages:
                    self.log.warning("! "+error_message)
                    # auxdata handler is persistent, therefore errors status
                    # needs to be reset before next orbit
                    auxclass.error.reset()

            auxdata_error_status.append(auxclass.error.status)
            auxdata_error_codes.extend(auxclass.error.codes)

            self.log.info("- %s auxdata handler completed" % (auxdata_type.upper()))

        # Return error status list
        return auxdata_error_status, auxdata_error_codes

    # TODO: Marked as obsolete
    def _classify_surface_types(self, l1b, l2):
        """ Run the surface type classification """
        pyclass = self.l2def.surface_type.pyclass
        surface_type = get_surface_type_class(pyclass)
        surface_type.set_options(**self.l2def.surface_type.options)
        surface_type.classify(l1b, l2)
        l2.set_surface_type(surface_type.result)

    # TODO: Marked as obsolete
    def _validate_surface_types(self, l2):
        """ Loop over stack of surface type validators """
        surface_type_validators = self.l2def.validator.surface_type
        error_codes = ["l2proc_surface_type_discarded"]
        error_states = []
        error_messages = []
        for name, validator_def in list(surface_type_validators.items()):
            validator = get_validator(validator_def["pyclass"])
            validator.set_options(**validator_def["options"])
            state, message = validator.validate(l2)
            error_states.append(state)
            error_messages.append(message)
            if state:
                self.log.info("- Validator message: "+message)
        error_status = True in error_states
        return error_status, error_codes

    # TODO: Marked as obsolete
    def _waveform_retracking(self, l1b, l2):
        """ Retracking: Obtain surface elevation from l1b waveforms """
        # loop over retrackers for each surface type

        for surface_type, retracker_def in list(self.l2def.retracker.items()):

            # Check if any waveforms need to be retracked for given
            # surface type
            surface_type_flag = l2.surface_type.get_by_name(surface_type)
            if surface_type_flag.num == 0:
                self.log.info("- no waveforms of type %s" % surface_type)
                continue

            # Benchmark retracker performance
            # XXX: is currently the bottleneck of level2 processing
            timestamp = time.time()

            # Retrieve the retracker assiciated with surface type
            # from the l2 settings
            retracker = get_retracker_class(retracker_def["pyclass"])

            # Set options (if any)
            if retracker_def["options"] is not None:
                retracker.set_options(**retracker_def["options"])

            # set subset of waveforms
            retracker.set_indices(surface_type_flag.indices)

            # Add classifier data (some retracker need that)
            retracker.set_classifier(l1b.classifier)

            # Start the retracking
            retracker.retrack(l1b, l2)

            # Retrieve the range after retracking
            l2.update_retracked_range(retracker)

            # XXX: Let the retracker return other parameters?
            l2.set_radar_mode(l1b.waveform.radar_mode)

            # retrieve potential error status and update surface type flag
            if retracker.error_flag.num > 0:
                l2.surface_type.add_flag(retracker.error_flag.flag, "invalid")
            self.log.info("- Retrack class %s with %s in %.3f seconds" % (
                surface_type, retracker_def["pyclass"],
                time.time()-timestamp))

        # Error handling not yet implemented, return dummy values
        return False, None

    # TODO: Marked as obsolete
    def _estimate_sea_surface_height(self, l2):

        # 2. get get sea surface anomaly
        ssa = get_l2_ssh_class(self.l2def.ssa.pyclass)
        ssa.set_options(**self.l2def.ssa.options)
        ssa.interpolate(l2)

        # dedicated setters, else the uncertainty, bias attributes are broken
        l2.ssa.set_value(ssa.value)
        l2.ssa.set_uncertainty(ssa.uncertainty)

    # TODO: Marked as obsolete
    def _get_altimeter_freeboard(self, l1b, l2):
        """ Compute radar freeboard and its uncertainty """

        afrbalg = get_frb_algorithm(self.l2def.afrb.pyclass)
        if self.l2def.afrb.options is not None:
            afrbalg.set_options(**self.l2def.afrb.options)
        afrb, afrb_unc = afrbalg.get_radar_freeboard(l1b, l2)

        # Check and return error status and codes
        # (unlikely in this case)
        error_status = afrbalg.error.status
        error_codes = afrbalg.error.codes

        if not error_status:
            # Add to l2data
            l2.afrb.set_value(afrb)
            l2.afrb.set_uncertainty(afrb_unc)

        # on error: display error messages as warning and return status flag
        # (this will cause the processor to report and skip this orbit segment)
        else:
            error_messages = self._snow.error.get_all_messages()
            for error_message in error_messages:
                self.log.warning("! "+error_message)

        return error_status, error_codes

    # TODO: Marked as obsolete
    def _get_freeboard_from_radar_freeboard(self, l1b, l2):
        """ Convert the altimeter freeboard in radar freeboard """

        frbgeocorr = get_frb_algorithm(self.l2def.frb.pyclass)
        frbgeocorr.set_options(**self.l2def.frb.options)
        frb, frb_unc = frbgeocorr.get_freeboard(l1b, l2)

        # Check and return error status and codes (e.g. missing file)
        error_status = frbgeocorr.error.status
        error_codes = frbgeocorr.error.codes

        # Add to l2data
        if not error_status:
            # Add to l2data
            l2.frb.set_value(frb)
            l2.frb.set_uncertainty(frb_unc)

        # on error: display error messages as warning and return status flag
        # (this will cause the processor to report and skip this orbit segment)
        else:
            error_messages = frbgeocorr.get_all_messages()
            for error_message in error_messages:
                self.log.warning("! "+error_message)

    # TODO: Marked as obsolete
    def _apply_freeboard_filter(self, l2):
        """ Apply freeboard filters as defined in the level-2 settings file
        under `root.filter.freeboard`

        Filtering means:
        - setting the freeboard value to nan
        - setting the surface type classification to invalid
        """

        #TODO: Transform this method into optional processing item

        # Loop over freeboard filters
        for name, filter_def in list(self.l2def.filter.freeboard.items()):

            # Get corresponding class name in pysiral.filter and transfer options
            # XXX: This should be rewritten as (e.g.)
            #   `frbfilter = VariableFilter(filter_def.pyclass, **filter_def.options)`
            frbfilter = get_filter(filter_def["pyclass"])
            frbfilter.set_options(**filter_def["options"])

            # XXX: This is a temporary fix of an error in the algorithm
            #
            # Explanation: The filter target was wrongly set to radar freeboard,
            # meaning that whether a freeboard value was filtered was determined on
            # the wrong parameter. Both values differ by the geometric snow propagation
            # correction (22% of snow depth). While the impact on the high freeboard end
            # is negligible, at the lower (negative) end more freeboard where filtered
            # than necessary since radar freeboard is always lower.
            #
            # The `afrb` filter target was hard coded, thus an option is added to replace
            # the filter target (`root.filter.freeboard.frb_valid_range.filter_target`).
            # The default option is the wrong one only for consistency reasons.
            filter_target = "afrb"
            if "filter_target" in filter_def["options"]:
                filter_target = filter_def["options"]["filter_target"]

            # Check if action is required
            frbfilter.apply_filter(l2, filter_target)
            if frbfilter.flag.num == 0:
                continue

            # Logging
            self.log.info("- Filter message: %s has flagged %g waveforms" % (
                filter_def["pyclass"], frbfilter.flag.num))

            # Set surface type flag (contains invalid)
            l2.surface_type.add_flag(frbfilter.flag.flag, "invalid")

            # Remove invalid elevations / freeboards
            l2.frb.set_nan_indices(frbfilter.flag.indices)

    # TODO: Marked as obsolete
    def _convert_freeboard_to_thickness(self, l2):
        """
        Convert Freeboard to Thickness
        Note: This step requires the definition of sea ice density
              (usually in the l2 settings)
        """

        frb2sit = get_sit_algorithm(self.l2def.sit.pyclass)
        frb2sit.set_options(**self.l2def.sit.options)

        sit, sit_unc, ice_dens, ice_dens_unc = frb2sit.get_thickness(l2)

        # Check and return error status and codes (e.g. missing file)
        error_status = frb2sit.error.status

        # Add to l2data
        if not error_status:
            # Add to l2data
            l2.sit.set_value(sit)
            l2.sit.set_uncertainty(sit_unc)
            l2.set_auxiliary_parameter("idens", "sea_ice_density", ice_dens, ice_dens_unc)

    # TODO: Marked as obsolete
    def _apply_thickness_filter(self, l2):
        for name, filter_def in list(self.l2def.filter.thickness.items()):
            sitfilter = get_filter(filter_def["pyclass"])
            sitfilter.set_options(**filter_def["options"])
            sitfilter.apply_filter(l2, "sit")
            if sitfilter.flag.num == 0:
                continue
            self.log.info("- Filter message: %s has flagged %g waveforms" % (
                filter_def["pyclass"], sitfilter.flag.num))
            # Set surface type flag (contains invalid)
            l2.surface_type.add_flag(sitfilter.flag.flag, "invalid")
            # Remove invalid thickness values
            l2.sit.set_nan_indices(sitfilter.flag.indices)

    # TODO: Marked as obsolete
    def _post_processing_items(self, l2):
        """

        :param l2:
        :return:
        """
        # Get the post processing options
        post_processing_items = self._l2def.get("post_processing", None)
        if post_processing_items is None:
            self.log.info("No post-processing items defined")
            return

        # Get the list of post-processing items
        for pp_item in post_processing_items:
            pp_class = get_cls(pp_item["module_name"], pp_item["class_name"], relaxed=False)
            post_processor = pp_class(**pp_item["options"])
            post_processor.apply(l2)
            msg = "- Level-2 post-processing item `%s` applied" % (pp_item["label"])
            self.log.info(msg)

    def _create_l2_outputs(self, l2):
        for output_handler in self._output_handler:
            output = Level2Output(l2, output_handler)
            self.log.info("- Write {} data file: {}".format(output_handler.id, output.export_filename))

    def _add_to_orbit_collection(self, l2):
        self._orbit.append(l2)


class Level2ProductDefinition(DefaultLoggingClass):
    """ Main configuration class for the Level-2 Processor """

    def __init__(self, run_tag, l2_settings_file):

        super(Level2ProductDefinition, self).__init__(self.__class__.__name__)
        self.error = ErrorStatus(self.__class__.__name__)

        # Mandatory parameter
        self._run_tag = run_tag
        self._l2_settings_file = l2_settings_file
        self._parse_l2_settings()

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

    @property
    def run_tag(self):
        return self._run_tag

    @property
    def l2def(self):
        return self._l2def

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
            self.log.warning("Unknown error code (%s), ignoring" % error_code)

    def write_to_file(self, output_id, directory):
        """ Write a summary file to the defined export directory """

        # Create a simple filename
        filename = Path(directory) / "pysiral-l2proc-summary.txt"
        self.log.info("Exporting summary report: %s" % filename)

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


class Level2ProcessorStep(DefaultLoggingClass):
    """
    Parent class for any Level-2 processor step class, which may be distributed over the
    different pysiral modules.

    This class also serves as a template for all sub-classes. Mandatory methods and properties
    in this class which raise a NotImplementedException must be overwritten by the subclass
    """

    def __init__(self, cfg):
        """
        Init the
        :param cfg:
        """
        super(Level2ProcessorStep, self).__init__(self.__class__.__name__)

        # -- Properties --

        # Class configuration
        self.cfg = cfg

        # Log messages
        self.msgs = []

        # Error Status
        self.error = ErrorStatus()

        # Error flag dict {code: component}
        self.error_flag_bit_dict = {
            "l1b": 0,
            "l2proc": 0,
            "auxdata": 2,
            "surface_type": 3,
            "retracker": 4,
            "range_correction": 5,
            "frb": 6,
            "sit": 7,
            "filter": 8,
            "other": 17}

    def execute(self, l1b, l2):
        """
        The main entry point for the
        :param l1b:
        :param l2:
        :return:
        """

        # Execute the method of the subclass. The class needs to
        error_flag = self.execute_step(l1b, l2)

        # Update the status flag
        self.update_error_flag(l2, error_flag)

    def update_error_flag(self, l2, error_flag):
        """
        Add the error_flag of the the processing step to the
        :param: l2: The Level-2 data container
        :param: error_flag: An array with the shape of l2.records containing the error flag
            (False: nominal, True: error)
        :return:
        """
        raise NotImplementedError("Deactivated for test purposes")

    def execute_step(self, l1b, l2):
        raise NotImplementedError("This method needs to implemented in {}".format(self.classname))

    @property
    def classname(self):
        return self.__class__.__name__


class Level2ProcessorStepOrder(DefaultLoggingClass):
    """
    A container providing the ordered list of processing steps
    as initialized classes for each trajectory
    """

    def __init__(self, cfg):
        """
        Initialize this class
        :param cfg: the procsteps tag from the Level-2 processor definitions file
        """
        super(Level2ProcessorStepOrder, self).__init__(self.__class__.__name__)

        # Properties
        self.cfg = cfg
        self.error = ErrorStatus()

        # A list of the class object (not initialized!)
        self._classes = []
        self.get_procstep_classes()

    def get_procstep_classes(self):
        """
        Retrieves the required classes from the processor definition files and stores them in a list
        without initializing them. This way a freshly initialized version can be supplied to each
        l2 data object without risk of interference of class properties
        :return:
        """

        # Loop
        for procstep_item in self.cfg:

            # Extract the information
            procstep_type = list(procstep_item.keys())[0]
            procstep_def = procstep_item[procstep_type]

            # Get the pyclass
            obj = get_cls(procstep_type, procstep_def["pyclass"])

            # This object should not be None
            if obj is None:
                msg = "Could not find L2 processing step class: {}.{}".format(procstep_type, procstep_def["pyclass"])
                self.error.add_error("missing-class", msg)
                self.error.raise_on_error()

            breakpoint()

    def validate(self):
        """
        Checkout the difference processing steps and validate input/output variables in
        the order of the steps
        :return:
        """

        breakpoint()
