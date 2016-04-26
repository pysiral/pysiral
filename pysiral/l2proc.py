# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 14:04:27 2015

@author: Stefan
"""
from pysiral.config import td_branches
from pysiral.l1bdata import L1bdataNCFile
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
        self._l2_processing_of_orbit_files()
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
        self.initialized = True
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

# %% Level2Processor: orbit processing

    def _l2_processing_of_orbit_files(self):
        """ Orbit-wise level2 processing """
        # TODO: Evaluate parallelization
        self.log.info("Start Orbit Processing")

        # loop over l1bdata preprocessed orbits
        for i, l1b_file in enumerate(self._l1b_files):

            # Log the current position in the file stack
            self.log.info("+ l1b orbit file %g of %g (%.2f%%)" % (
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
            # TODO: move to level1bData class
            self._apply_range_corrections(l1b)

            # Initialize the orbit level-2 data container
            l2 = Level2Data(l1b)

            # Add sea ice concentration (can be used as classifier)
            self._get_sea_ice_concentration(l2)

            # Surface type classification (ocean, ice, lead, ...)
            # (ice type classification comes later)
            # TODO: Add L2 classifiers (ice concentration, ice type)
            self._classify_surface_types(l1b, l2)

            # Validate surface type classification
            # yes/no decision on continuing with orbit
            error_status, error_messages = self._validate_surface_types(l2)
            if error_status:
                for error_message in error_messages:
                    self.log.info("- validator message: "+error_message)
                self.log.info("- skip file")
                continue

            # Get elevation by retracking of different surface types
            # adds parameter elevation to l2
            self._waveform_range_retracking(l1b, l2)

            # Compute the sea surface anomaly (from mss and lead tie points)
            # adds parameter ssh, ssa, afrb to l2
            self._estimate_ssh_and_radar_freeboard(l2)

            # Get sea ice type (may be required for geometrical corrcetion)
            self._get_sea_ice_type(l2)

            # Get snow depth & density
            self._get_snow_parameters(l2)

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
        self.log.info("- Parsing l1bdata file: %s" % l1b_file)
        l1b = L1bdataNCFile(l1b_file)
        l1b.parse()
        l1b.info.subset_region_name = self._job.roi.hemisphere
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

    def _get_sea_ice_concentration(self, l2):
        """ Get sea ice concentration along track from auxdata """
        sic, msg = self._sic.get_along_track_sic(l2)
        if not msg == "":
            self.log.info("- "+msg)
        # Add to l2data
        l2.sic.set_value(sic)

    def _classify_surface_types(self, l1b, l2):
        """ Run the surface type classification """
        pyclass = self._job.config.surface_type.pyclass
        surface_type = get_surface_type_class(pyclass)
        surface_type.set_options(**self._job.config.surface_type.options)
        # Add all classifiers from l1bdata
        for classifier_name in l1b.classifier.parameter_list:
            classifier = getattr(l1b.classifier, classifier_name)
            surface_type.add_classifiers(classifier, classifier_name)
        # add sea ice concentration
        surface_type.add_classifiers(l2.sic, "sic")
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

    def _waveform_range_retracking(self, l1b, l2):
        """ Retracking: Obtain surface elevation from l1b waveforms """
        # loop over retrackers for each surface type
        surface_types, retracker_def = td_branches(self._job.config.retracker)
        for i, surface_type in enumerate(surface_types):
            surface_type_flag = l2.surface_type.get_by_name(surface_type)
            if surface_type_flag.num == 0:
                self.log.info("- no waveforms of type %s" % surface_type)
                continue
            timestamp = time.time()
            retracker = get_retracker_class(retracker_def[i].pyclass)
            retracker.set_options(**retracker_def[i].options)
            retracker.set_indices(surface_type_flag.indices)
            retracker.retrack(l1b, l2)
            l2.update_retracked_range(retracker)
            self.log.info("- Retrack class %s with %s in %.3f seconds" % (
                surface_type, retracker_def[i].pyclass,
                time.time()-timestamp))

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

    def _get_sea_ice_type(self, l2):
        """ Get sea ice concentration along track from auxdata """
        sitype, msg = self._sitype.get_along_track_sitype(l2)
        if not msg == "":
            self.log.info("- "+msg)
        # Add to l2data
        l2.sitype.set_value(sitype)

    def _get_snow_parameters(self, l2):
        """ Get snow depth and density """
        snow_depth, snow_dens, msg = self._snow.get_along_track_snow(l2)
        if not msg == "":
            self.log.info("- "+msg)
        # Add to l2data
        l2.snow_depth.set_value(snow_depth)
        l2.snow_dens.set_value(snow_dens)

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

    def _post_processing(self, l2):
        pass

    def _create_outputs(self, l2):
        pass

    def _add_to_orbit_collection(self, l2):
        self._orbit.append(l2)
