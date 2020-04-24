# -*- coding: utf-8 -*-

"""
"""

__author__ = "Stefan Hendricks"

import numpy as np
from pathlib import Path
from cftime import num2pydate

from pysiral import psrlcfg
from pysiral.clocks import StopWatch
from pysiral.errorhandler import ErrorStatus
from pysiral.ers.sgdrfile import ERSSGDR
from pysiral.flag import ORCondition
from pysiral.l1bdata import Level1bData
from pysiral.logging import DefaultLoggingClass
from pysiral.surface_type import ESA_SURFACE_TYPE_DICT


class ERSReaperSGDR(DefaultLoggingClass):
    """ Converts a Envisat SGDR object into a L1bData object """

    def __init__(self, cfg, raise_on_error=False):
        """
        Input handler for Sentinel-3 L2WAT netCDF files from the CODA.
        :param cfg: A treedict object (root.input_handler.options) from the corresponding Level-1 pre-processor
                    config file
        :param raise_on_error: Boolean value if the class should raise an exception upon an error (default: False)
        """

        cls_name = self.__class__.__name__
        super(ERSReaperSGDR, self).__init__(cls_name)
        self.error = ErrorStatus(caller_id=cls_name)

        # Store arguments
        self.raise_on_error = raise_on_error
        self.cfg = cfg

        # Debug variables
        self.timer = None

    def get_l1(self, filepath, polar_ocean_check=None):
        """
        Read the Envisat SGDR file and transfers its content to a Level1Data instance
        :param filepath: The full file path to the netCDF file
        :param polar_ocean_check: Mandatory parameter but will be ignored as ERS Data is full orbit
        :return: The parsed (or empty) Level-1 data container
        """

        # Store arguments
        self.filepath = filepath

        # Create an empty Level-1 data object
        self.l1 = Level1bData()

        #  for debug purposes
        self.timer = StopWatch()
        self.timer.start()

        # Read the file
        # NOTE: This will create the variable `self.sgdr`
        self._read_sgdr()

        # Get metadata
        self._set_input_file_metadata()

        # Polar ocean check passed, now fill the rest of the l1 data groups
        self._set_l1_data_groups()

        self.timer.stop()
        self.log.info("- Created L1 object in %.3f seconds" % self.timer.get_seconds())

        return self.l1

    def _read_sgdr(self):
        """ Read the L1b file and create a ERS native L1b object """
        self.sgdr = ERSSGDR(self.cfg)
        self.sgdr.filename = self.filepath
        self.sgdr.parse()
        error_status = self.sgdr.get_status()
        if error_status:
            # TODO: Needs ErrorHandler
            raise IOError()
        self.sgdr.post_processing()

    def _set_input_file_metadata(self):
        """ Extract essential metadata information from SGDR file """
        info = self.l1.info
        sgdr = self.sgdr
        info.set_attribute("pysiral_version", psrlcfg.version)
        try:
            info.set_attribute("mission", self.cfg.platform_name_dict[str(sgdr.nc.mission)])
        except KeyError:
            mission_id = self.sgdr.guess_mission_from_filename()
            info.set_attribute("mission", self.cfg.platform_name_dict[str(mission_id)])

        info.set_attribute("mission_data_version", sgdr.nc.software_ver)
        info.set_attribute("orbit", sgdr.nc.abs_orbit)
        info.set_attribute("cycle", sgdr.nc.cycle)
        mission_data_source = Path(sgdr.nc.filename).name
        info.set_attribute("mission_data_source", mission_data_source)
        info.set_attribute("timeliness", self.cfg.timeliness)

    def _set_l1_data_groups(self):
        self._transfer_timeorbit()            # (lon, lat, alt, time)
        self._transfer_waveform_collection()  # (power, range)
        self._transfer_range_corrections()    # (range corrections)
        self._transfer_surface_type_data()    # (land flag, ocean flag, ...)
        self._transfer_classifiers()          # (beam parameters, flags, ...)

    def _transfer_timeorbit(self):
        """ Extracts the time/orbit data group from the SGDR data """

        # Transfer the orbit position
        self.l1.time_orbit.set_position(
            self.sgdr.nc.lon_20hz.flatten(),
            self.sgdr.nc.lat_20hz.flatten(),
            self.sgdr.nc.alt_20hz.flatten())

        # Transfer the timestamp
        sgdr_timestamp = self.sgdr.nc.time_20hz.flatten()
        units = self.cfg.sgdr_timestamp_units
        calendar = self.cfg.sgdr_timestamp_calendar
        timestamp = num2pydate(sgdr_timestamp, units, calendar)
        self.l1.time_orbit.timestamp = timestamp

        # Mandatory antenna pointing parameter (but not available for ERS)
        dummy_angle = np.full(timestamp.shape, 0.0)
        self.l1.time_orbit.set_antenna_attitude(dummy_angle, dummy_angle, dummy_angle)

        # Update meta data container
        self.l1.update_data_limit_attributes()

    def _transfer_waveform_collection(self):
        """ Transfers the waveform data (power & range for each range bin) """

        # Transfer the reformed 18Hz waveforms
        self.l1.waveform.set_waveform_data(
            self.sgdr.wfm_power,
            self.sgdr.wfm_range,
            self.sgdr.radar_mode)

        # Set valid flag to exclude calibration data
        # (see section 3.5 of Reaper handbook)
        tracking_state = self.sgdr.nc.alt_state_flag_20hz.flatten()
        valid = ORCondition()
        valid.add(tracking_state == 2)
        valid.add(tracking_state == 3)
        self.l1.waveform.set_valid_flag(valid.flag)

    def _transfer_range_corrections(self):
        """
        Transfer range correction data from the SGDR netCDF to the
        l1bdata object. The parameter are defined in
        config/mission_def.yaml for ers1/ers2
        -> ersX.settings.sgdr_range_correction_targets

        For a description of the parameter see section 3.10 in the
        REAPER handbook
        """
        grc_dict = self.cfg.range_correction_targets
        for name in grc_dict.keys():
            target_parameter = grc_dict[name]
            if target_parameter is None:
                continue
            correction = getattr(self.sgdr.nc, target_parameter)
            correction = np.repeat(correction, self.cfg.sgdr_n_blocks)
            self.l1.correction.set_parameter(name, correction)

    def _transfer_classifiers(self):
        """
        Transfer classifier parameter from the SGDR netCDF to the
        l1bdata object. Most parameter are defined in
        config/mission_def.yaml for ers1/ers2
        -> ersX.settings.sgdr_range_correction_targets
        """
        target_dict = self.cfg.classifier_targets
        for parameter_name in target_dict.keys():
            nc_parameter_name = target_dict[parameter_name]
            nc_parameter = getattr(self.sgdr.nc, nc_parameter_name)
            self.l1.classifier.add(nc_parameter.flatten(), parameter_name)

    def _transfer_surface_type_data(self):
        surface_type = self.sgdr.nc.surface_type
        surface_type = np.repeat(surface_type, self.cfg.sgdr_n_blocks)
        for key in ESA_SURFACE_TYPE_DICT.keys():
            flag = surface_type == ESA_SURFACE_TYPE_DICT[key]
            self.l1.surface_type.add_flag(flag, key)

    @property
    def empty(self):
        """
        Default return object, if nodata should be returned
        :return: Representation of an empty object (None)
        """
        return None