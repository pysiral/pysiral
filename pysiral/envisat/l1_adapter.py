# -*- coding: utf-8 -*-

"""
"""

__author__ = "Stefan Hendricks"

import re
import numpy as np
from attrdict import AttrDict
from loguru import logger
from typing import Union, Dict
from scipy import interpolate
from cftime import num2pydate
from pathlib import Path

from pysiral import psrlcfg
from pysiral.clocks import StopWatch
from pysiral.l1preproc import Level1PInputHandlerBase
from pysiral.envisat.functions import (get_envisat_window_delay, get_envisat_wfm_range)
from pysiral.errorhandler import ErrorStatus
from pysiral.iotools import ReadNC
from pysiral.l1bdata import Level1bData
from pysiral.core import DefaultLoggingClass
from pysiral.core.flags import ESA_SURFACE_TYPE_DICT


class EnvisatSGDRNC(Level1PInputHandlerBase):
    """ Converts a Envisat SGDR object into a L1bData object """

    def __init__(self,
                 cfg: Union[Dict, AttrDict],
                 raise_on_error: bool = False
                 ) -> None:
        """
        Input handler for Sentinel-3 L2WAT netCDF files from the CODA.

        :param cfg: Options from the corresponding Level-1 pre-processor config file
        :param raise_on_error: Boolean value if the class should raise an exception upon an error (default: False)
        """

        cls_name = self.__class__.__name__
        super(EnvisatSGDRNC, self).__init__(cfg, raise_on_error, cls_name)

        # Debug variables
        self.timer = None

        # Properties
        self.sgdr = None
        self.filepath = None
        self.l1 = None

    def get_l1(self,
               filepath: Union[str, Path]
               ) -> "Level1bData":
        """
        Read the Envisat SGDR file and transfers its content to a Level1Data instance
        :param filepath: The full file path to the netCDF file
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
        self.sgdr = ReadNC(self.filepath, nan_fill_value=True)

        # Get metadata
        self._set_input_file_metadata()

        # Polar ocean check passed, now create all l1 data groups
        self._set_l1_data_groups()

        self.timer.stop()
        logger.info("- Created L1 object in %.3f seconds" % self.timer.get_seconds())

        return self.l1

    def _set_input_file_metadata(self) -> None:
        """ Extract essential metadata information from SGDR file """
        info = self.l1.info
        sgdr = self.sgdr
        info.set_attribute("pysiral_version", psrlcfg.version)
        info.set_attribute("mission", "envisat")
        info.set_attribute("mission_data_version", sgdr.software_version)
        info.set_attribute("orbit", sgdr.absolute_orbit_number)
        info.set_attribute("cycle", sgdr.cycle_number)
        info.set_attribute("mission_data_source", sgdr.product_name)
        info.set_attribute("timeliness", self.cfg.timeliness)

    def _set_l1_data_groups(self) -> None:
        """
        Transfer data from the necessary data groups
        :return:
        """
        self._transfer_timeorbit()            # (lon, lat, alt, time)
        self._transfer_waveform_collection()  # (power, range)
        self._transfer_range_corrections()    # (range corrections)
        self._transfer_surface_type_data()    # (land flag, ocean flag, ...)
        self._transfer_classifiers()          # (beam parameters, flags, ...)

    def _transfer_timeorbit(self) -> None:
        """
        Extracts the time/orbit data group from the SGDR data
        :return:
        """

        # Transfer the orbit position
        self.l1.time_orbit.set_position(self.sgdr.lon_20, self.sgdr.lat_20, self.sgdr.alt_20)

        # Transfer the timestamp
        sgdr_timestamp = self.sgdr.time_20
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
        """
        Transfers the waveform data (power & range for each range bin)
        :return:
        """

        # Transfer the reformed 18Hz waveforms
        # "waveform samples (I2+Q2, 1/2048 FFT power unit): 18 Hz Ku band";
        # "the echo is corrected for the intermediate frequency filter effect";
        wfm_power = self.sgdr.waveform_fft_20_ku
        n_records, n_range_bins = wfm_power.shape

        # In Envisat data v3.0 the tracker range is already corrected for all internal effects,
        # but the nominal tracking bin is variable.
        nominal_tracking_bin = 63 + self.sgdr.offset_tracking_20.values/256
        window_delay_m = self.sgdr.tracker_range_20_ku - nominal_tracking_bin * self.cfg.bin_width_meter

        # Compute the range value for each range bin of the 18hz waveform
        wfm_range = get_envisat_wfm_range(
            window_delay_m,
            n_range_bins,
            bin_width_meter=self.cfg.bin_width_meter
        )

        # Transfer data to the waveform group
        self.l1.waveform.set_waveform_data(wfm_power, wfm_range, self.cfg.radar_mode)

        # Set valid flag to exclude calibration data
        # (see section 3.5 of Reaper handbook)
        valid_flag = np.logical_not(self.sgdr.waveform_fault_id_20.astype(bool))
        self.l1.waveform.set_valid_flag(valid_flag)

    def _transfer_range_corrections(self):
        """
        Transfer range correction data from the SGDR netCDF to the
        l1bdata object. The parameters are defined in
        config/mission_def.yaml for ers1/ers2
        -> settings.sgdr_range_correction_targets

        For a description of the parameter see section 3.10 in the
        REAPER handbook
        """

        # Get the reference times for interpolating the range corrections from 1Hz -> 20Hz
        time_1Hz = np.array(self.sgdr.time_01)
        time_20Hz = np.array(self.sgdr.time_20)

        # Loop over all range correction in config file
        grc_dict = self.cfg.range_correction_targets
        for name in grc_dict.keys():

            # Get the variable
            target_parameter = grc_dict[name]
            if target_parameter is None:
                continue
            correction = np.array(getattr(self.sgdr, target_parameter))

            # Debug code
            # -> in this case discard the variable
            n_nans = len(np.where(np.isnan(correction))[0])
            # TODO: 500 is hardcoded
            if 500 < n_nans < len(correction):
                msg = "Significant number of NaNs (%g) in range correction variable: %s"
                msg %= (n_nans, target_parameter)
                logger.warning(msg)
            elif n_nans == len(correction):
                msg = "All-NaN array encountered in range correction variable: %s"
                msg %= target_parameter
                logger.warning(msg)

            # Some Envisat range corrections are 1Hz others 20Hz
            # -> Those with "_01" in the variable name need to be
            # extrapolated to 20 Hz
            error = False
            if re.search(self.cfg.variable_identifier_1Hz, target_parameter):
                correction, error = self.interp_1Hz_to_20Hz(correction, time_1Hz, time_20Hz,
                                                            fill_on_error_value=0.0)
            if error:
                msg = f"Failing to create 20Hz range correction variable for {target_parameter}"
                logger.warning(msg)

            # Interpolate NaN's or return a zero-filled array for all-nan input variables
            correction_filtered = self.find_and_interpolate_nans(correction, fill_on_error_value=0.0)

            # Debug code
            # -> in this case discard the variable
            n_nans = len(np.where(np.isnan(correction_filtered))[0])
            if n_nans > 0:
                msg = "Remaining NaN's after filtering in %s" % target_parameter
                logger.warning(msg)

            # Set the parameter
            self.l1.correction.set_parameter(name, correction_filtered)

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
            nc_parameter = getattr(self.sgdr, nc_parameter_name)
            self.l1.classifier.add(nc_parameter.flatten(), parameter_name)

    def _transfer_surface_type_data(self):
        surface_type = self.sgdr.surf_class_20
        for key in ESA_SURFACE_TYPE_DICT.keys():
            flag = surface_type == ESA_SURFACE_TYPE_DICT[key]
            self.l1.surface_type.add_flag(flag, key)

    @staticmethod
    def find_and_interpolate_nans(variable, fill_on_error_value=np.nan):
        """
        Replace NaN's in variable with linear interpolated values
        :param variable:
        :return: interpolated variable
        """
        is_nan = np.isnan(variable)
        n_nans = len(np.where(is_nan)[0])
        if n_nans == 0:
            variable_filtered = variable
        else:
            x = np.arange(len(variable))
            valid = np.where(np.logical_not(is_nan))[0]
            try:
                f = interpolate.interp1d(x[valid], variable[valid], fill_value="extrapolate",
                                         bounds_error=False)
                variable_filtered = f(x)
            except ValueError:
                variable_filtered = np.full(variable.shape, fill_on_error_value)
        return variable_filtered

    @staticmethod
    def interp_1Hz_to_20Hz(variable_1Hz, time_1Hz, time_20Hz, fill_on_error_value=np.nan, **kwargs):
        """
        Computes a simple linear interpolation to transform a 1Hz into a 20Hz variable
        :param variable_1Hz: an 1Hz variable array
        :param time_1Hz: 1Hz reference time
        :param time_20Hz: 20 Hz reference time
        :return: the interpolated 20Hz variable
        """
        error_status = False
        try:
            is_valid = np.logical_and(np.isfinite(time_1Hz), np.isfinite(variable_1Hz))
            valid_indices = np.where(is_valid)[0]
            f = interpolate.interp1d(time_1Hz[valid_indices], variable_1Hz[valid_indices],
                                     fill_value="extrapolate", bounds_error=False, **kwargs)
            variable_20Hz = f(time_20Hz)
        except ValueError:
            variable_20Hz = np.full(time_20Hz.shape, fill_on_error_value)
            error_status = True
        return variable_20Hz, error_status

    @property
    def empty(self):
        """
        Default return object, if nodata should be returned
        :return: Representation of an empty object (None)
        """
        return None