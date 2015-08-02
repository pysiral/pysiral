# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 14:04:27 2015

@author: Stefan
"""
from pysiral.config import td_branches
from pysiral.l1bdata import L1bConstructor
from pysiral.l2data import Level2Data
from pysiral.mss import *
from pysiral.roi import *
from pysiral.surface_type import *
from pysiral.retracker import *

import os


class Level2Processor(object):

    def __init__(self, job):

        self._job = job
        self._orbit = []
        self._l1b_files = []
        self._config = None
        self._initialized = False
        self._error_handler = {"raise_on_error": True}

    @property
    def orbit(self):
        return self._orbit

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

    def _initialize_processor(self):
        """ Read required auxiliary data sets """
        if self._initialized:
            return
        self._get_roi()
        # Read the mean surface height auxiliary file
        self._get_mss()
        # TODO: Read the ice concentration data
        # TODO: snow depth information
        # XXX: Ice Type Information?
        self.initialized = True

    def _get_roi(self):
        self._roi = globals()[self._job.roi.pyclass]()
        self._roi.set_options(**self._job.roi.options)

    def _get_mss(self):
        settings = self._job.config.ssh.mss
        local_repository = self._job.local_machine.auxdata_repository.static
        directory = local_repository[settings.local_machine_directory]
        filename = os.path.join(directory, settings.file)
        self._mss = globals()[settings.pyclass]()
        self._mss.set_filename(filename)
        self._mss.set_roi(self._roi)
        self._mss.parse()

    def _run_processor(self):
        """ Orbit-wise level2 processing """
        # TODO: Evaluate parallelization
        for l1b_file in self._l1b_files:
            l1b = self._read_l1b_file(l1b_file)        # Read the l1b file
            in_roi = self._trim_to_roi(l1b)            # region of interest
            over_ocean = self._trim_land_margins(l1b)  # trim edges (land)
            if not in_roi or not over_ocean:
                continue
            self._range_corrections(l1b)               # geophys. corrections
            l2 = Level2Data(l1b)
            # TODO: Get L2 classifiers (ice concentration, ice type)
            self._classify_surface_types(l1b, l2)      # (sea ice, lead, ...)
            self._retrack_waveforms(l1b, l2)           # (elevation)
            self._reference_to_ssh(l2)                 # (ssh, ssa, app frb)
            self._apply_data_quality_filter(l2)
            self._post_processing(l2)                  # (radar freeboard)
            self._create_outputs(l2)                   # (plots, files)
            self._add_to_orbit_collection(l2)

    def _clean_up(self):
        """ Make sure to deallocate memory """
        pass

    def _read_l1b_file(self, l1b_file):
        """ Read a L1b data file """
        l1b = L1bConstructor(self._config)
        l1b.mission = self._job.mission.id
        l1b.set_mission_options(**self._job.mission.options)
        l1b.filename = l1b_file
        l1b.construct()
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

    def _range_corrections(self, l1b):
        """ Apply the range corrections """
        for correction in self._job.config.corrections:
            l1b.apply_range_correction(correction)

    def _classify_surface_types(self, l1b, l2):
        """ Run the surface type classificator """
        surface_type = globals()[self._job.config.surface_type.pyclass]()
        surface_type.set_options(**self._job.config.surface_type.options)
        surface_type.set_classifiers(l1b.classifier)
        # XXX: Not sure if this is useful:
        # if hasattr(l1b, "surface_type"):
        #    classifier.set_initial_classification(l1b.surface_type)
        surface_type.classify()
        l2.set_surface_type(surface_type.result)

    def _retrack_waveforms(self, l1b, l2):
        """ Do the retracking """
        # loop over retrackers for each surface type
        surface_types, retracker_def = td_branches(self._job.config.retracker)
        for i, surface_type in enumerate(surface_types):
            retracker = globals()[retracker_def[i].pyclass]()
            retracker.set_options(**retracker_def[i].options)
            retracker.set_surface_type(surface_type)
            retracker.retrack(l1b, l2)
            l2.update_retracked_range(retracker)

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

    def _apply_data_quality_filter(self, l2):
        pass

    def _post_processing(self, l2):
        pass

    def _create_outputs(self, l2):
        pass

    def _add_to_orbit_collection(self, l2):
        self._orbit.append(l2)
