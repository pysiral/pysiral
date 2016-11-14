# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 14:04:27 2015

@author: Stefan
"""

from pysiral.config import (td_branches, ConfigInfo, TimeRangeRequest,
                            get_yaml_config, PYSIRAL_VERSION, HOSTNAME)
from pysiral.errorhandler import ErrorStatus, PYSIRAL_ERROR_CODES
from pysiral.l1bdata import L1bdataNCFile
from pysiral.iotools import get_local_l1bdata_files
from pysiral.l2data import Level2Data
from pysiral.logging import DefaultLoggingClass
from pysiral.sic import get_l2_sic_handler
from pysiral.sitype import get_l2_sitype_handler
from pysiral.snow import get_l2_snow_handler
from pysiral.mss import get_l2_ssh_class
from pysiral.roi import get_roi_class
from pysiral.surface_type import get_surface_type_class
from pysiral.retracker import get_retracker_class
from pysiral.filter import get_filter
from pysiral.validator import get_validator
from pysiral.frb import get_frb_algorithm
from pysiral.sit import get_sit_algorithm
from pysiral.output import get_output_class
from pysiral.path import filename_from_path

from collections import deque, OrderedDict
from datetime import datetime
import numpy as np
import time
import glob
import sys
import os


class Level2Processor(DefaultLoggingClass):

    def __init__(self, job):

        super(Level2Processor, self).__init__(self.__class__.__name__)

        # Error Status Handler
        self.error = ErrorStatus(caller_id=self.__class__.__name__)

        # job definition
        self._job = job

        # List of Level-2 (processed) orbit segments
        self._orbit = deque()

        # List of Level-1b input files
        self._l1b_files = []

        # pysiral config
        self._config = ConfigInfo()

        # Processor Initialization Flag
        self._initialized = False

        # Processor summary report
        self.report = L2ProcessorReport()


# %% Level2Processor: class properties

    @property
    def orbit(self):
        return self._orbit

    @property
    def has_empty_file_list(self):
        return len(self._l1b_files) == 0


# %% Level2Processor: public methods

    def initialize(self):
        self._initialize_processor()
        self._initialize_summary_report()

    def get_input_files_local_machine_def(self, time_range, version="default"):
        mission_id = self._job.mission_id
        hemisphere = self._job.hemisphere
        l1b_files = get_local_l1bdata_files(
            mission_id, time_range, hemisphere, version=version)
        self.set_l1b_files(l1b_files)

        # Update the report
        self.report.n_files = len(l1b_files)
        self.report.time_range = time_range
        self.report.l1b_repository = os.path.split(l1b_files[0])[0]

    def set_l1b_files(self, l1b_files):
        self._l1b_files = l1b_files

    def remove_old_l2data(self, time_range):
        """ Clean up old l2 output data """
        # can be several oututs
        output_ids, output_defs = td_branches(self._job.config.output)
        for output_id, output_def in zip(output_ids, output_defs):
            output = get_output_class(output_def.pyclass)
            output.set_options(**output_def.options)
            output.set_base_export_path(output_def.path)
            export_folder = output.get_full_export_path(time_range.start)

            # Get list of output files
            search_pattern = os.path.join(export_folder, "*.*")
            l2output_files = glob.glob(search_pattern)

            # Delete files
            self.log.info("Removing %g output files [ %s ] in %s" % (
                len(l2output_files), output_id, export_folder))
            for l2output_file in l2output_files:
                os.remove(l2output_file)

    def run(self):
        """ Run the processor """
        self._initialize_processor()
        self._l2_processing_of_orbit_files()
        self._l2proc_summary_to_file()
        self._clean_up()

    def purge(self):
        """ Clean the orbit collection """
        pass

# %% Level2Processor: house keeping methods

    def _l2proc_summary_to_file(self):
        output_ids, output_defs = td_branches(self._job.config.output)
        for output_id, output_def in zip(output_ids, output_defs):
            output = get_output_class(output_def.pyclass)
            output.set_options(**output_def.options)
            output.set_base_export_path(output_def.path)
            time_range = self.report.time_range
            export_folder = output.get_full_export_path(time_range.start)
            self.report.write_to_file(output_id, export_folder)

    def _clean_up(self):
        """ All procedures that need to be reset after a run """
        self.report.clean_up()

# %% Level2Processor: initialization

    def _initialize_processor(self):
        """ Read required auxiliary data sets """

        # Instance can be reused
        if self._initialized:
            # Empty orbit list (or else orbits will acculumate)
            self._orbit.clear()
            return

        self.log.info("Initializing processor")

        self.log.info("Processor Settings - range correction list:")
        for correction in self._job.config.corrections:
            self.log.info("- %s" % correction)
        self.log.info("Processor Settings - surface type classificator: %s" % (
            self._job.config.surface_type.pyclass))
        self.log.info("Processor Settings - lead interpolator: %s" % (
            self._job.config.ssa.pyclass))

        # Set the region of interest option
        # (required for MSS subsetting)
        self._set_roi()

        # Load static background field

        # Read the mean surface height auxiliary file
        self._set_mss()

        # Handler for dynamic data sets (sea ice concentration, ...)
        # need to be called with timestamps and positions

        # Sea ice concentration data handler
        self._set_sic_handler()

        # sea ice type data handler (needs to be before snow)
        self._set_sitype_handler()

        # snow data handler (needs to provide snow depth and density)
        self._set_snow_handler()

        # All done
        self._initialized = True
        self.log.info("Initializing done")

    def _set_roi(self):
        self.log.info("Processor Settings - ROI: %s" % self._job.roi.pyclass)
        self._roi = get_roi_class(self._job.roi.pyclass)
        self._roi.set_options(**self._job.roi.options)

    def _set_mss(self):
        """ Loading the mss product file from a static file """
        settings = self._job.config.auxdata.mss
        filename = os.path.join(settings.local_repository, settings.file)
        self.log.info("Processor Settings - MSS: %s" % settings.pyclass)
        self.log.info("- loading roi subset from: %s" % filename)
        self._mss = get_l2_ssh_class(settings.pyclass)
        self._mss.set_filename(filename)
        self._mss.set_roi(self._roi)
        self._mss.parse()

    def _set_sic_handler(self):
        """ Set the sea ice concentration handler """
        settings = self._job.config.auxdata.sic
        self._sic = get_l2_sic_handler(settings.pyclass)
        if settings.options is not None:
            self._sic.set_options(**settings.options)
        self._sic.set_local_repository(settings.local_repository)
        self._sic.set_filenaming(settings.filenaming)
        self._sic.set_subfolders(settings.subfolders)
        self._sic.initialize()
        self.log.info("Processor Settings - SIC handler: %s" % (
            settings.pyclass))

    def _set_sitype_handler(self):
        """ Set the sea ice type handler """
        settings = self._job.config.auxdata.sitype
        self._sitype = get_l2_sitype_handler(settings.pyclass)
        if settings.options is not None:
            self._sitype.set_options(**settings.options)
        self._sitype.set_local_repository(settings.local_repository)
        self._sitype.set_filenaming(settings.filenaming)
        self._sitype.set_subfolders(settings.subfolders)
        self.log.info("Processor Settings - SIType handler: %s" % (
            settings.pyclass))

    def _set_snow_handler(self):
        """ Set the snow (depth and density) handler """
        settings = self._job.config.auxdata.snow
        self._snow = get_l2_snow_handler(settings.pyclass)
        if settings.options is not None:
            self._snow.set_options(**settings.options)
        self._snow.set_local_repository(settings.local_repository)
        self._snow.set_filenaming(settings.filenaming)
        self._snow.set_subfolders(settings.subfolders)
        self.log.info("Processor Settings - Snow handler: %s" % (
            settings.pyclass))

    def _initialize_summary_report(self):
        """
        Only add report parameter that are not time range specific
        (e.g. the processor l2 settings)
        """
        self.report.l2_settings_file = self._job.l2_settings_file

# %% Level2Processor: orbit processing

    def _l2_processing_of_orbit_files(self):
        """ Orbit-wise level2 processing """
        # TODO: Evaluate parallelization
        self.log.info("Start Orbit Processing")

        # loop over l1bdata preprocessed orbits
        for i, l1b_file in enumerate(self._l1b_files):

            # Log the current position in the file stack
            self.log.info("+ [ %g of %g ] (%.2f%%)" % (
                i+1, len(self._l1b_files),
                float(i+1)/float(len(self._l1b_files))*100.))

            # Read the the level 1b file (l1bdata netCDF is required)
            l1b = self._read_l1b_file(l1b_file)

            # Apply the geophysical range corrections on the waveform range
            # bins in the l1b data container
            # TODO: move to level1bData class
            self._apply_range_corrections(l1b)

            # Initialize the orbit level-2 data container
            l2 = Level2Data(l1b)

            # Add sea ice concentration (can be used as classifier)
            error_status, error_codes = self._get_sea_ice_concentration(l2)
            if error_status:
                self._discard_l1b_procedure(error_codes, l1b_file)
                continue

            # Get sea ice type (may be required for geometrical corrcetion)
            error_status, error_codes = self._get_sea_ice_type(l2)
            if error_status:
                self._discard_l1b_procedure(error_codes, l1b_file)
                continue

            # Surface type classification (ocean, ice, lead, ...)
            # (ice type classification comes later)
            # TODO: Add L2 classifiers (ice concentration, ice type)
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
            self._estimate_ssh_and_radar_freeboard(l2)

            # Get snow depth & density
            error_status, error_codes = self._get_snow_parameters(l2)
            if error_status:
                self.report.add_orbit_discarded_event(error_codes, l1b_file)
                continue

            # get radar(-derived) from altimeter freeboard
            self._get_freeboard_from_radar_freeboard(l2)

            # Apply freeboard filter
            self._apply_freeboard_filter(l2)

            # Convert to thickness
            self._convert_freeboard_to_thickness(l2)

            # Filter thickness
            self._apply_thickness_filter(l2)

            # Create output files
            self._create_l2_outputs(l2)

            # Add data to orbit stack
            self._add_to_orbit_collection(l2)

    def _read_l1b_file(self, l1b_file):
        """ Read a L1b data file (l1bdata netCDF) """
        filename = filename_from_path(l1b_file)
        self.log.info("- Parsing l1bdata file: %s" % filename)
        l1b = L1bdataNCFile(l1b_file)
        l1b.parse()
        l1b.info.subset_region_name = self._job.roi.hemisphere
        return l1b

    def _discard_l1b_procedure(self, error_codes, l1b_file):
        """ Log and report discarded l1b orbit segment """
        self.log.info("- skip file")
        for error_code in error_codes:
            self.report.add_orbit_discarded_event(error_code, l1b_file)

    def _apply_range_corrections(self, l1b):
        """ Apply the range corrections """
        for correction in self._job.config.corrections:
            l1b.apply_range_correction(correction)

    def _get_sea_ice_concentration(self, l2):
        """ Get sea ice concentration along track from auxdata """

        # Get along-track sea ice concentrations via the SIC handler class
        # (see self._set_sic_handler)
        sic, msg = self._sic.get_along_track_sic(l2)

        # Report any messages from the SIC handler
        if not msg == "":
            self.log.info("- "+msg)

        # Check and return error status and codes (e.g. missing file)
        error_status = self._sic.error.status
        error_codes = self._sic.error.codes

        # No error: Set sea ice concentration data to the l2 data container
        if not error_status:
            l2.sic.set_value(sic)

        # on error: display error messages as warning and return status flag
        # (this will cause the processor to report and skip this orbit segment)
        else:
            error_messages = self._sic.error.get_all_messages()
            for error_message in error_messages:
                self.log.warning("! "+error_message)
                # SIC Handler is persistent, therefore errors status
                # needs to be reset before next orbit
                self._sic.error.reset()

        return error_status, error_codes

    def _get_sea_ice_type(self, l2):
        """ Get sea ice type (myi fraction) along track from auxdata """

        # Call the sitype handler
        sitype, msg = self._sitype.get_along_track_sitype(l2)

        # Report any messages from the sitype handler
        if not msg == "":
            self.log.info("- "+msg)

        # Check and return error status and codes (e.g. missing file)
        error_status = self._sitype.error.status
        error_codes = self._sitype.error.codes

        # Add to l2data
        if not error_status:
            l2.sitype.set_value(sitype)

        # on error: display error messages as warning and return status flag
        # (this will cause the processor to report and skip this orbit segment)
        else:
            error_messages = self._sitype.error.get_all_messages()
            for error_message in error_messages:
                self.log.warning("! "+error_message)
                # SIC Handler is persistent, therefore errors status
                # needs to be reset before next orbit
            self._sitype.error.reset()

        return error_status, error_codes

    def _classify_surface_types(self, l1b, l2):
        """ Run the surface type classification """
        pyclass = self._job.config.surface_type.pyclass
        surface_type = get_surface_type_class(pyclass)
        surface_type.set_options(**self._job.config.surface_type.options)
        surface_type.classify(l1b, l2)
        l2.set_surface_type(surface_type.result)

    def _validate_surface_types(self, l2):
        """ Loop over stack of surface type validators """
        surface_type_validators = self._job.config.validator.surface_type
        names, validators = td_branches(surface_type_validators)
        error_codes = ["l2proc_surface_type_discarded"]
        error_states = []
        error_messages = []
        for name, validator_def in zip(names, validators):
            validator = get_validator(validator_def.pyclass)
            validator.set_options(**validator_def.options)
            state, message = validator.validate(l2)
            error_states.append(state)
            error_messages.append(message)
            if state:
                self.log.info("- Validator message: "+message)
        error_status = True in error_states
        return error_status, error_codes

    def _waveform_retracking(self, l1b, l2):
        """ Retracking: Obtain surface elevation from l1b waveforms """
        # loop over retrackers for each surface type
        surface_types, retracker_def = td_branches(self._job.config.retracker)

        for i, surface_type in enumerate(surface_types):

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
            retracker = get_retracker_class(retracker_def[i].pyclass)

            # Set options (if any)
            if retracker_def[i].options is not None:
                retracker.set_options(**retracker_def[i].options)

            # set subset of waveforms
            retracker.set_indices(surface_type_flag.indices)

            # Add classifier data (some retracker need that)
            retracker.set_classifier(l1b.classifier)

            # Start the retracking
            retracker.retrack(l1b, l2)

            # Retrieve the range after retracking
            l2.update_retracked_range(retracker)

            # XXX: Let the retracker return other parameters?

            # retrieve potential error status and update surface type flag
            if retracker.error_flag.num > 0:
                l2.surface_type.add_flag(retracker.error_flag.flag, "invalid")
            self.log.info("- Retrack class %s with %s in %.3f seconds" % (
                surface_type, retracker_def[i].pyclass,
                time.time()-timestamp))

        # Error handling not yet implemented, return dummy values
        return False, None

    def _estimate_ssh_and_radar_freeboard(self, l2):
        # 1. get mss for orbit
        l2.mss = self._mss.get_track(l2.track.longitude, l2.track.latitude)
        # 2. get get sea surface anomaly
        ssa = get_l2_ssh_class(self._job.config.ssa.pyclass)
        ssa.set_options(**self._job.config.ssa.options)
        ssa.interpolate(l2)
        # dedicated setters, else the uncertainty, bias attributes are broken
        l2.ssa.set_value(ssa.value)
        l2.ssa.set_uncertainty(ssa.uncertainty)
        # altimeter freeboard (from radar altimeter w/o corrections)
        l2.afrb = l2.elev - l2.mss - l2.ssa

    def _get_snow_parameters(self, l2):
        """ Get snow depth and density """
        snow_depth, snow_dens, msg = self._snow.get_along_track_snow(l2)
        if not msg == "":
            self.log.info("- "+msg)
        # Add to l2data
        l2.snow_depth.set_value(snow_depth)
        l2.snow_dens.set_value(snow_dens)
        # XXX: Error Handling not yet implemted, return dummies
        return False, None

    def _get_freeboard_from_radar_freeboard(self, l2):
        """ Convert the altimeter freeboard in radar freeboard """
        frbgeocorr = get_frb_algorithm(self._job.config.frb.pyclass)
        frbgeocorr.set_options(**self._job.config.frb.options)
        frb, msg = frbgeocorr.get_freeboard(l2)
        if not msg == "":
            self.log.info("- "+msg)
        l2.frb.set_value(frb)

    def _apply_freeboard_filter(self, l2):
        freeboard_filters = self._job.config.filter.freeboard
        names, filters = td_branches(freeboard_filters)
        for name, filter_def in zip(names, filters):
            frbfilter = get_filter(filter_def.pyclass)
            frbfilter.set_options(**filter_def.options)
            frbfilter.apply_filter(l2, "afrb")
            if frbfilter.flag.num == 0:
                continue
            self.log.info("- Filter message: %s has flagged %g waveforms" % (
                filter_def.pyclass, frbfilter.flag.num))
            # Set surface type flag (contains invalid)
            l2.surface_type.add_flag(frbfilter.flag.flag, "invalid")
            # Remove invalid elevations / freeboards
            l2.frb[frbfilter.flag.indices] = np.nan

    def _convert_freeboard_to_thickness(self, l2):
        frb2sit = get_sit_algorithm(self._job.config.sit.pyclass)
        frb2sit.set_options(**self._job.config.sit.options)
        sit, ice_dens, msg = frb2sit.get_thickness(l2)
        if not msg == "":
            self.log.info("- "+msg)
        l2.sit.set_value(sit)
        l2.ice_dens.set_value(ice_dens)

    def _apply_thickness_filter(self, l2):
        thickness_filters = self._job.config.filter.thickness
        names, filters = td_branches(thickness_filters)
        for name, filter_def in zip(names, filters):
            sitfilter = get_filter(filter_def.pyclass)
            sitfilter.set_options(**filter_def.options)
            sitfilter.apply_filter(l2, "sit")
            if sitfilter.flag.num == 0:
                continue
            self.log.info("- Filter message: %s has flagged %g waveforms" % (
                filter_def.pyclass, sitfilter.flag.num))
            # Set surface type flag (contains invalid)
            l2.surface_type.add_flag(sitfilter.flag.flag, "invalid")
            # Remove invalid elevations / freeboards
            l2.frb[sitfilter.flag.indices] = np.nan

    def _create_l2_outputs(self, l2):
        output_ids, output_defs = td_branches(self._job.config.output)
        for output_id, output_def in zip(output_ids, output_defs):
            output = get_output_class(output_def.pyclass)
            output.set_options(**output_def.options)
            output.set_processor_settings(self._job.settings.level2)
            output.set_base_export_path(output_def.path)
            output.write_to_file(l2)
            self.log.info("- Write %s data file: %s" % (
                output_id, output.filename))

    def _add_to_orbit_collection(self, l2):
        self._orbit.append(l2)


class L2ProcJob(DefaultLoggingClass):
    """ Container for the definition and handling of a pre-processor job """

    def __init__(self):

        super(L2ProcJob, self).__init__(self.__class__.__name__)

        # Error Status
        self.error = ErrorStatus(caller_id=self.__class__.__name__)

        # Save pointer to pysiral configuration
        self.pysiral_config = ConfigInfo()

        # Initialize the time range and set to monthly per default
        self.time_range = TimeRangeRequest()

        # Initialize job parameter
        self.options = L2ProcJobOptions()

        # Level-2 processor settings from external file
        self.settings = None

        # List for iterations (currently only month-wise)
        self.iterations = []

    def parse_l2_settings(self):
        """ Read the Level-2 settings yaml file """

        # Check if l2_settings (id or filename is already known)
        # a) self.options.from_dict must have been called
        # b) self.options.l2_settings must have been manually set
        if self.options.l2_settings is None:
            self.error.add_error(
                "no-l2-settings", "L2 settings file unspecified (is None)")
            return

        # Check
        if os.path.isfile(self.options.l2_settings):
            l2_settings_filename = self.options.l2_settings
        # if not filename, than it need to be id of settings file in
        # pysiral\config\l2
        else:
            settings_path = os.path.join(
                self.pysiral_config.pysiral_local_path, "settings", "l2")
            l2_settings_filename = os.path.join(
                settings_path, self.options.l2_settings+".yaml")
            if not os.path.isfile(l2_settings_filename):
                self.error.add_error(
                    "l2-settings-not-found",
                    "Level-2 yaml settings file not found: %s" % (
                        l2_settings_filename))
                return

        # All clear, read the settings
        self.settings = get_yaml_config(l2_settings_filename)
        self.options.l2_settings_filename = l2_settings_filename
        self.log.info("Level-2 settings file: %s" % l2_settings_filename)

    def generate_preprocessor_iterations(self):
        """ Break the requested time range into monthly iterations """

        # The input data is organized in folder per month, therefore
        # the processing period is set accordingly
        self.time_range.set_period("monthly")
        self.log.info("Level-2 processor base period is monthly")
        self.time_range.set_exclude_month(self.options.exclude_month)
        self.log.info("Excluding month: %s" % str(self.options.exclude_month))
        self.iterations = self.time_range.get_iterations()
        self.log.info("Number of iterations: %g" % len(self.iterations))

    def process_requested_time_range(self):
        """
        Verify time range with mission data availability and create
        datetime objects
        """

        mission_info = self.pysiral_config.get_mission_info(self.mission_id)

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

    def validate(self):
        self._validate_and_expand_auxdata()
        self._validate_and_create_output_directory()

    def _validate_and_expand_auxdata(self):
        """
        - Verifies auxdata information in config.auxdata with
          content of config/auxdata_def.yaml
        - Transfer relevant content of auxdata_def.yaml into configuration
          structure
        """

        # Loop over all auxiliary data types
        auxtypes, auxinfos = td_branches(self.settings.level2.auxdata)
        for auxtype, auxinfo in zip(auxtypes, auxinfos):

            self.log.info("Validating local auxiliary repository [%s:%s]" % (
                auxtype, auxinfo.name))

            # Locate auxdata setting in config/auxdata_def.yaml
            try:
                pysiral_def = self.pysiral_config.auxdata[auxtype][auxinfo.name]
            except:
                msg = "id %s for type %s not in config/auxdata_def.yaml" % (
                    auxinfo.name, auxtype)
                self.error.add_error("invalid-auxdata-def", msg)

            # Repace local machine directory placeholder with actual directory
            auxdata_def = self.pysiral_config.local_machine.auxdata_repository
            auxdata_id = pysiral_def.local_repository

            # Specific type might not need a local repository, but update
            # anyway
            if auxdata_id is None:
                self.settings.level2.auxdata[auxtype].update(pysiral_def)
                if self.settings.level2.auxdata[auxtype].has_key("source"):
                    del self.settings.level2.auxdata[auxtype].source
                continue

            # Check if entry is in local_machine_def.yaml
            try:
                local_repository = auxdata_def[auxtype][auxdata_id]
                pysiral_def.local_repository = local_repository
            except:
                msg = "No auxdata definition in local_machine_def.yaml" + \
                      " for %s:%s" % (auxtype, auxdata_id)
                self.error.add_error("missing-auxdata-def", msg)
                continue

            # Check if local directory (only main) does exist
            # XXX: Check for each iteration?
            if os.path.isdir(local_repository):
                # Expand the settings with actual path
                self.settings.level2.auxdata[auxtype].update(pysiral_def)
                if self.settings.level2.auxdata[auxtype].has_key("source"):
                    del self.settings.level2.auxdata[auxtype].source
            else:
                msg = "Missing local auxiliary directory (%s:%s): %s " % (
                          auxtype, auxdata_id, local_repository)
                self.error.add_error("missing-auxdata-dir", msg)

    def _validate_and_create_output_directory(self):

        output_ids, output_defs = td_branches(self.settings.level2.output)
        export_path = self.pysiral_config.local_machine.product_repository
        export_path = os.path.join(export_path, self.options.run_tag)
        time = datetime.now()
        tstamp = time.strftime("%Y%m%dT%H%M%S")
        for output_id, output_def in zip(output_ids, output_defs):
            if self.overwrite_protection:
                product_export_path = os.path.join(
                   export_path, tstamp, output_id)
            else:
                product_export_path = os.path.join(export_path, output_id)
            self.config.output[output_id].path = product_export_path
            self.log.info("Exporting %s data in directory %s" % (
                output_id, product_export_path))
        self.settings.export_path = product_export_path

    @property
    def mission_id(self):
        return self.settings.mission.id

    @property
    def hemisphere(self):
        return self.settings.roi.hemisphere

    @property
    def input_version(self):
        return self.options.input_version

    @property
    def remove_old(self):
        return self.options.remove_old

    @property
    def overwrite_protection(self):
        return self.options.overwrite_protection

    @property
    def config(self):
        return self.settings.level2

    @property
    def roi(self):
        return self.settings.roi

    @property
    def l2_settings_file(self):
        return self.options.l2_settings_filename


class L2ProcJobOptions(object):
    """ Simple container for Level-2 processor options """

    def __init__(self):
        self.l2_settings_filename = "none"
        self.l2_settings = None
        self.run_tag = None
        self.input_version = "default"
        self.remove_old = False
        self.overwrite_protection = True
        self.start_date = None
        self.stop_date = None
        self.exclude_month = []
        self.export_path = None

    def from_dict(self, options_dict):
        for parameter in options_dict.keys():
            if hasattr(self, parameter):
                setattr(self, parameter, options_dict[parameter])


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
        filename = os.path.join(directory, "pysiral-l2proc-summary.txt")
        self.log.info("Exporting summary report: %s" % filename)

        lfmt = "  %-16s : %s\n"
        current_time = str(datetime.now()).split(".")[0]
        with open(filename, "w") as fhandle:

            # Write infos on settings, host, os, ....
            fhandle.write("# pysiral Level2Processor Summary\n\n")
            fhandle.write(lfmt % ("created", current_time))

            # Brief statistics of files, errors, warnings
            fhandle.write("\n# Processor Statistics\n\n")
            fhandle.write(lfmt % ("l1b files", str(self.n_files)))
            fhandle.write(lfmt % ("errors", str(self.n_discarded_files)))
            fhandle.write(lfmt % ("warnings", str(self.n_warnings)))

            fhandle.write("\n# Processor & Local Machine Settings\n\n")
            fhandle.write(lfmt % ("pysiral version", PYSIRAL_VERSION))
            fhandle.write(lfmt % ("python version", sys.version))
            fhandle.write(lfmt % ("hostname", HOSTNAME))

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
                    fn = filename_from_path(discarded_file)
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
