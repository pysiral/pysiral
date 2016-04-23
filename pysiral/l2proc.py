# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 14:04:27 2015

@author: Stefan
"""
from pysiral.config import td_branches
from pysiral.l1bdata import L1bdataNCFile
from pysiral.l2data import Level2Data
from pysiral.logging import DefaultLoggingClass
from pysiral.mss import *
from pysiral.roi import *
from pysiral.surface_type import *
from pysiral.retracker import *
from pysiral.filter import get_filter
from pysiral.validator import get_validator

from collections import deque
import numpy as np
import time
import os


class Level2Processor(DefaultLoggingClass):

    def __init__(self, job):

        super(Level2Processor, self).__init__("Level2Processor")
        self._job = job
        self._orbit = deque()
        self._l1b_files = []
        self._config = None
        self._initialized = False
        self._error_handler = {"raise_on_error": True}

# %% Level2Processor: class properties

    @property
    def orbit(self):
        return self._orbit

# %% Level2Processor: public methods

    def error_handling(self, **keyw):
        self._error_handler.update(keyw)

    def set_config(self, config):
        self._config = config

    def set_l1b_files(self, l1b_files):
        self._l1b_files = l1b_files

    def run(self):
        """ Run the processor """
        self._initialize_processor()
        self._run_processor()
        self._clean_up()

    def purge(self):
        """ Clean the orbit collection """
        pass

    def grid_orbit_collection(self, grid_definitions):
        pass

# %% Level2Processor: house keeping methods

    def _clean_up(self):
        """ Make sure to deallocate memory """
        pass

# %% Level2Processor: initialization

    def _initialize_processor(self):
        """ Read required auxiliary data sets """

        # Instance can be reused
        if self._initialized:
            return

        self.log.info("Initializing processor")

        self.log.info("Processor Options - range corrections:")
        for correction in self._job.config.corrections:
            self.log.info("- %s" % correction)
        self.log.info("Processor Options - surface type classificator: %s" % (
            self._job.config.surface_type.pyclass))
        self.log.info("Processor Options - lead interpolator: %s" % (
            self._job.config.ssh.ssa.pyclass))

        # Set the region of interest option
        # (required for MSS subsetting)
        self._get_roi()

        # Read the mean surface height auxiliary file
        self._get_mss()

        # TODO: Read the ice concentration data
        #       getter function, only reloads if necessary

        # TODO: snow depth information
        #       requires getter (only reload data if necessary)

        # TODO: Ice Type Information?
        #       same as SIC -> getter function

        self.initialized = True
        self.log.info("Initializing done")

    def _get_roi(self):
        self.log.info("Setting ROI type: %s" % self._job.roi.pyclass)
        self._roi = globals()[self._job.roi.pyclass]()
        self._roi.set_options(**self._job.roi.options)

    def _get_mss(self):
        settings = self._job.config.ssh.mss
        local_repository = self._job.local_machine.auxdata_repository.static
        directory = local_repository[settings.local_machine_directory]
        filename = os.path.join(directory, settings.file)
        self.log.info("Loading mss (%s) file: %s" % (
            settings.pyclass, filename))
        self._mss = globals()[settings.pyclass]()
        self._mss.set_filename(filename)
        self._mss.set_roi(self._roi)
        self._mss.parse()

# %% Level2Processor: orbit processing

    def _run_processor(self):
        """ Orbit-wise level2 processing """
        # TODO: Evaluate parallelization
        self.log.info("Start Orbit Processing")

        # loop over l1bdata preprocessed orbits
        for i, l1b_file in enumerate(self._l1b_files):

            # Log the current position in the file stack
            self.log.info("l1b orbit file %g of %g (%.2f%%)" % (
                i+1, len(self._l1b_files),
                float(i+1)/float(len(self._l1b_files))*100.))

            # Read the the level 1b file (l1bdata netCDF is required)
            l1b = self._read_l1b_file(l1b_file)

            # File subsetting
            # XXX: This is obsolete due to pre-processing, keep for later?
            # in_roi = self._trim_to_roi(l1b)
            # over_ocean = self._trim_land_margins(l1b)
            # if not in_roi or not over_ocean:
            #     continue

            # Apply the geophysical range corrections on the waveform range
            # bins in the l1b data container
            self._apply_range_corrections(l1b)

            # Initialize the orbit level-2 data container
            l2 = Level2Data(l1b)

            # Surface type classification (ocean, ice, lead, ...)
            # (ice type classification comes later)
            # TODO: Add L2 classifiers (ice concentration, ice type)
            self._classify_surface_types(l1b, l2)

            # Validate surface type classification
            # yes/no decision on continuing with orbit
            error_status, error_messages = self._validate_surface_types(l2)
            if error_status:
                for error_message in error_messages:
                    self.log.info(". validator message: "+error_message)
                self.log.info(". skip file")
                continue

            # Get elevation by retracking of different surface types
            # adds parameter elevation to l2
            self._retrack_waveforms(l1b, l2)

            # Compute the sea surface anomaly (from mss and lead tie points)
            # adds parameter ssh, ssa, afrb to l2
            self._reference_to_ssh(l2)

            # Apply freeboard filter
            self._apply_freeboard_filter(l2)

            # self._apply_data_quality_filter(l2)
            # self._post_processing(l2)
            # self._create_outputs(l2)
            self._add_to_orbit_collection(l2)

    def _read_l1b_file(self, l1b_file):
        """ Read a L1b data file (l1bdata netCDF) """
        self.log.info(". Parsing l1bdata file: %s" % l1b_file)
        l1b = L1bdataNCFile(l1b_file)
        l1b.parse()
        return l1b

    def _trim_to_roi(self, l1b):
        """
        Trims orbit for ice/ocean areas (excluding land, land ice)
        and returns true if data in ROI and false otherwise
        """
        roi_list = self._roi.get_roi_list(
            l1b.time_orbit.longitude, l1b.time_orbit.latitude)
        if len(roi_list) == 0:  # No match
            return False
        if len(roi_list) == l1b.n_records:  # Full match (no trimming)
            return True
        l1b.trim_to_subset(roi_list)  # Partial match (trimming)
        return True

    def _trim_land_margins(self, l1b):
        """
        Trim land areas at the margins of the orbit data

        points over land surrounded by ocean measurements (e.g. inside the
        Canadian Archipelago) will be excluded later in the processing
        """
        if not l1b.surface_type.has_flag("ocean"):
            # TODO: Think of method to classify lon/lat point as land/ocean
            #       (e.g. basemap functionality?) if land flag is not part
            #       of l1b data set
            pass
        # TODO: This assumes that in the initial classification all waveforms
        #       of relevance are classified as ocean
        flag = l1b.surface_type.ocean.flag
        is_ocean_list = np.nonzero(flag)
        if len(is_ocean_list) == l1b.n_records:  # Nothing to do here
            return True
        if len(is_ocean_list) == 0:              # Nothing left to do
            return False
        trimmed_list = np.arange(np.amin(is_ocean_list),
                                 np.amax(is_ocean_list)+1)
        l1b.trim_to_subset(trimmed_list)
        return True

    def _apply_range_corrections(self, l1b):
        """ Apply the range corrections """
        for correction in self._job.config.corrections:
            l1b.apply_range_correction(correction)

    def _classify_surface_types(self, l1b, l2):
        """ Run the surface type classificator """
        surface_type = globals()[self._job.config.surface_type.pyclass]()
        surface_type.set_options(**self._job.config.surface_type.options)
        surface_type.set_classifiers(l1b.classifier)
        surface_type.set_l1b_surface_type(l1b.surface_type)
        surface_type.classify()
        l2.set_surface_type(surface_type.result)

    def _validate_surface_types(self, l2):
        """ Loop over stack of surface type validators """
        surface_type_validators = self._job.config.validator.surface_type
        names, validators = td_branches(surface_type_validators)
        error_states = []
        error_messages = []
        for name, validator_def in zip(names, validators):
            validator = get_validator(validator_def.pyclass)
            validator.set_options(**validator_def.options)
            state, message = validator.validate(l2)
            error_states.append(state)
            error_messages.append(message)
        error_status = True in error_states
        return error_status, error_messages

    def _retrack_waveforms(self, l1b, l2):
        """ Retracking: Obtain surface elevation from l1b waveforms """
        # loop over retrackers for each surface type
        surface_types, retracker_def = td_branches(self._job.config.retracker)
        for i, surface_type in enumerate(surface_types):
            surface_type_flag = l2.surface_type.get_by_name(surface_type)
            if surface_type_flag.num == 0:
                self.log.info(". no waveforms of type %s" % surface_type)
                continue
            timestamp = time.time()
            retracker = globals()[retracker_def[i].pyclass]()
            retracker.set_options(**retracker_def[i].options)
            retracker.set_indices(surface_type_flag.indices)
            retracker.retrack(l1b, l2)
            l2.update_retracked_range(retracker)
            self.log.info(". Retrack class %s with %s in %.3f seconds" % (
                surface_type, retracker_def[i].pyclass,
                time.time()-timestamp))

    def _reference_to_ssh(self, l2):
        # 1. get mss for orbit
        l2.mss = self._mss.get_track(l2.track.longitude, l2.track.latitude)
        # 2. get get sea surface anomaly
        ssa = globals()[self._job.config.ssh.ssa.pyclass]()
        ssa.set_options(**self._job.config.ssh.ssa.options)
        ssa.interpolate(l2)
        # dedicated setters, else the uncertainty, bias attributes are broken
        l2.ssa.set_value(ssa.value)
        l2.ssa.set_uncertainty(ssa.uncertainty)
        # get apparent freeboard
        l2.afrb = l2.elev - l2.mss - l2.ssa

    def _apply_freeboard_filter(self, l2):
        freeboard_filters = self._job.config.filter.freeboard
        names, filters = td_branches(freeboard_filters)
        for name, filter_def in zip(names, filters):
            frbfilter = get_filter(filter_def.pyclass)
            frbfilter.set_options(**filter_def.options)
            frbfilter.apply_filter(l2, "afrb")
            if frbfilter.flag.num == 0:
                continue
            self.log.info(". filter message: %s has flagged %g waveforms" % (
                filter_def.pyclass, frbfilter.flag.num))
            # Set surface type flag (contains invalid)
            l2.surface_type.add_flag(frbfilter.flag.flag, "invalid")
            # Remove invalid elevations / freeboards
            l2.range[frbfilter.flag.indices] = np.nan
            l2.elev[frbfilter.flag.indices] = np.nan
            l2.afrb[frbfilter.flag.indices] = np.nan

    def _post_processing(self, l2):
        pass

    def _create_outputs(self, l2):
        pass

    def _add_to_orbit_collection(self, l2):
        self._orbit.append(l2)
