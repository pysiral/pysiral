# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 15:10:04 2015

@author: Stefan
"""

from pysiral.cryosat2.functions import (
    tai2utc, get_tai_datetime_from_timestamp,
    get_cryosat2_wfm_power, get_cryosat2_wfm_range)
from pysiral.esa.functions import get_structarr_attr
from pysiral.cryosat2.l1bfile import CryoSatL1B
from pysiral.envisat.sgdrfile import EnvisatSGDR
from pysiral.helper import parse_datetime_str
from pysiral.classifier import CS2OCOGParameter, CS2PulsePeakiness

import numpy as np



ESA_SURFACE_TYPE_DICT = {
    "ocean": 0,
    "closed_sea": 1,
    "land_ice": 2,
    "land": 3}


class L1bAdapterCryoSat(object):
    """ Converts a CryoSat2 L1b object into a L1bData object """

    def __init__(self, config):
        self.filename = None
        self._config = config
        self._mission = "cryosat2"

    def construct_l1b(self, l1b):
        self.l1b = l1b                        # pointer to L1bData object
        self._read_cryosat2l1b()              # Read CryoSat-2 L1b data file
        self._transfer_metadata()             # (orbit, radar mode, ..)
        self._transfer_timeorbit()            # (lon, lat, alt, time)
        self._transfer_waveform_collection()  # (power, range)
        self._transfer_range_corrections()    # (range corrections)
        self._transfer_surface_type_data()    # (land flag, ocean flag, ...)
        self._transfer_classifiers()          # (beam parameters, flags, ...)

    def _read_cryosat2l1b(self):
        """ Read the L1b file and create a CryoSat-2 native L1b object """
        self.cs2l1b = CryoSatL1B()
        self.cs2l1b.filename = self.filename
        self.cs2l1b.parse()
        error_status = self.cs2l1b.get_status()
        if error_status:
            # TODO: Needs ErrorHandler
            raise IOError()
        self.cs2l1b.post_processing()

    def _transfer_metadata(self):
        self.l1b.info.mission = self._mission
        self.l1b.info.mission_data_version = self.cs2l1b.baseline
        self.l1b.info.radar_mode = self.cs2l1b.radar_mode
        self.l1b.info.orbit = self.cs2l1b.sph.abs_orbit_start
        self.l1b.info.start_time = parse_datetime_str(
            self.cs2l1b.sph.start_record_tai_time)
        self.l1b.info.stop_time = parse_datetime_str(
            self.cs2l1b.sph.stop_record_tai_time)

    def _transfer_timeorbit(self):
        # Transfer the orbit position
        longitude = get_structarr_attr(self.cs2l1b.time_orbit, "longitude")
        latitude = get_structarr_attr(self.cs2l1b.time_orbit, "latitude")
        altitude = get_structarr_attr(self.cs2l1b.time_orbit, "altitude")
        self.l1b.time_orbit.set_position(longitude, latitude, altitude)
        # Transfer the timestamp
        tai_objects = get_structarr_attr(
            self.cs2l1b.time_orbit, "tai_timestamp")
        tai_timestamp = get_tai_datetime_from_timestamp(tai_objects)
        utc_timestamp = tai2utc(tai_timestamp)
        self.l1b.time_orbit.timestamp = utc_timestamp

    def _transfer_waveform_collection(self):
        # Create the numpy arrays for power & range
        dtype = np.float32
        n_records = len(self.cs2l1b.waveform)
        n_range_bins = len(self.cs2l1b.waveform[0].wfm)
        echo_power = np.ndarray(shape=(n_records, n_range_bins), dtype=dtype)
        echo_range = np.ndarray(shape=(n_records, n_range_bins), dtype=dtype)
        # Set the echo power in dB and calculate range
        for i, record in enumerate(self.cs2l1b.waveform):
            echo_power[i, :] = get_cryosat2_wfm_power(
                np.array(record.wfm).astype(np.float32),
                record.linear_scale, record.power_scale)
            echo_range[i, :] = get_cryosat2_wfm_range(
                self.cs2l1b.measurement[i].window_delay, n_range_bins)
        # Transfer to L1bData
        self.l1b.waveform.add_waveforms(echo_power, echo_range)

    def _transfer_range_corrections(self):
        # Transfer all the correction in the list
        for key in self.cs2l1b.corrections[0].keys():
            if key in self._config.parameter.correction_list:
                self.l1b.correction.set_parameter(
                    key, get_structarr_attr(self.cs2l1b.corrections, key))
        # CryoSat-2 specific: Two different sources of ionospheric corrections
        options = self._config.get_mission_defaults(self._mission)
        key = options["ionospheric_correction_source"]
        ionospheric = get_structarr_attr(self.cs2l1b.corrections, key)
        self.l1b.correction.set_parameter("ionospheric", ionospheric)

    def _transfer_surface_type_data(self):
        # L1b surface type flag word
        surface_type = get_structarr_attr(
            self.cs2l1b.corrections, "surface_type")
        for key in ESA_SURFACE_TYPE_DICT.keys():
            flag = surface_type == self._SURFACE_TYPE_DICT[key]
            self.l1b.surface_type.add_flag(flag, key)

    def _transfer_classifiers(self):
        # Add L1b beam parameter group
        beam_parameter_list = [
            "stack_standard_deviation", "stack_centre",
            "stack_scaled_amplitude", "stack_skewness", "stack_kurtosis"]
        for beam_parameter_name in beam_parameter_list:
            recs = get_structarr_attr(self.cs2l1b.waveform, "beam")
            beam_parameter = [rec[beam_parameter_name] for rec in recs]
            self.l1b.classifier.add(beam_parameter, beam_parameter_name)
        # Calculate Parameters from waveform counts
        # XXX: This is a legacy of the CS2AWI IDL processor
        #      Threshold defined for waveform counts not power in dB
        wfm = get_structarr_attr(self.cs2l1b.waveform, "wfm")
        # Calculate the OCOG Parameter (CryoSat-2 notation)
        ocog = CS2OCOGParameter(wfm)
        self.l1b.classifier.add(ocog.width, "ocog_width")
        self.l1b.classifier.add(ocog.amplitude, "ocog_amplitude")
        # Calculate the Peakiness (CryoSat-2 notation)
        pulse = CS2PulsePeakiness(wfm)
        self.l1b.classifier.add(pulse.peakiness, "peakiness")
        self.l1b.classifier.add(pulse.peakiness_r, "peakiness_r")
        self.l1b.classifier.add(pulse.peakiness_l, "peakiness_l")


class L1bAdapterEnvisat(object):
    """ Converts a Envisat SGDR object into a L1bData object """

    def __init__(self, config):
        self.filename = None
        self._config = config
        self._mission = "envisat"

    def construct_l1b(self, l1b):
        self.l1b = l1b                        # pointer to L1bData object
        self._read_envisat_sgdr()             # Read Envisat SGDR data file
        self._transfer_metadata()             # (orbit, radar mode, ..)
        self._transfer_timeorbit()            # (lon, lat, alt, time)
        self._transfer_waveform_collection()  # (power, range)
        self._transfer_range_corrections()    # (range corrections)
        self._transfer_surface_type_data()    # (land flag, ocean flag, ...)
        self._transfer_classifiers()          # (beam parameters, flags, ...)

    def _read_envisat_sgdr(self):
        """ Read the L1b file and create a CryoSat-2 native L1b object """
        self.sgdr = EnvisatSGDR()
        self.sgdr.filename = self.filename
        self.sgdr.parse()
        error_status = self.sgdr.get_status()
        if error_status:
            # TODO: Needs ErrorHandler
            raise IOError()
        self.sgdr.post_processing()

    def _transfer_metadata(self):
        """ Extract essential metadata information from SGDR file """
        info = self.l1b.info
        sgdr = self.sgdr
        info.set_attribute("mission", self._mission)
        info.set_attribute("mission_data_version", "final v9.3p5")
        info.set_attribute("orbit", sgdr.mph.abs_orbit)
        info.set_attribute("cycle", sgdr.mph.cycle)
        info.set_attribute("cycle", sgdr.mph.cycle)
        info.set_attribute("mission_data_source", sgdr.mph.product)

    def _transfer_timeorbit(self):
        """ Extracts the time/orbit data group from the SGDR data """
        # Transfer the orbit position
        self.l1b.time_orbit.set_position(
            self.sgdr.mds_18hz.longitude,
            self.sgdr.mds_18hz.latitude,
            self.sgdr.mds_18hz.altitude)
        # Transfer the timestamp
        self.l1b.time_orbit.timestamp = self.sgdr.mds_18hz.timestamp
        # Update meta data container
        self.l1b.update_data_limit_attributes()

    def _transfer_waveform_collection(self):
        # Transfer to L1bData
        echo_power = self.sgdr.mds_18hz.power
        # XXX: Debug purposes only
        echo_range = self.sgdr.mds_18hz.range
        self.l1b.waveform.add_waveforms(echo_power, echo_range)

    def _transfer_range_corrections(self):
        pass
#        # Transfer all the correction in the list
#        for key in self.cs2l1b.corrections[0].keys():
#            if key in self._config.parameter.correction_list:
#                self.l1b.correction.set_parameter(
#                    key, get_structarr_attr(self.cs2l1b.corrections, key))
#        # CryoSat-2 specific: Two different sources of ionospheric corrections
#        options = self._config.get_mission_defaults(self._mission)
#        key = options["ionospheric_correction_source"]
#        ionospheric = get_structarr_attr(self.cs2l1b.corrections, key)
#        self.l1b.correction.set_parameter("ionospheric", ionospheric)

    def _transfer_surface_type_data(self):
        surface_type = self.sgdr.mds_18hz.surface_type
        for key in ESA_SURFACE_TYPE_DICT.keys():
            flag = surface_type == ESA_SURFACE_TYPE_DICT[key]
            self.l1b.surface_type.add_flag(flag, key)

    def _transfer_classifiers(self):
        pass
#        # Add L1b beam parameter group
#        beam_parameter_list = [
#            "stack_standard_deviation", "stack_centre",
#            "stack_scaled_amplitude", "stack_skewness", "stack_kurtosis"]
#        for beam_parameter_name in beam_parameter_list:
#            recs = get_structarr_attr(self.cs2l1b.waveform, "beam")
#            beam_parameter = [rec[beam_parameter_name] for rec in recs]
#            self.l1b.classifier.add(beam_parameter, beam_parameter_name)
#        # Calculate Parameters from waveform counts
#        # XXX: This is a legacy of the CS2AWI IDL processor
#        #      Threshold defined for waveform counts not power in dB
#        wfm = get_structarr_attr(self.cs2l1b.waveform, "wfm")
#        # Calculate the OCOG Parameter (CryoSat-2 notation)
#        ocog = CS2OCOGParameter(wfm)
#        self.l1b.classifier.add(ocog.width, "ocog_width")
#        self.l1b.classifier.add(ocog.amplitude, "ocog_amplitude")
#        # Calculate the Peakiness (CryoSat-2 notation)
#        pulse = CS2PulsePeakiness(wfm)
#        self.l1b.classifier.add(pulse.peakiness, "peakiness")
#        self.l1b.classifier.add(pulse.peakiness_r, "peakiness_r")
#        self.l1b.classifier.add(pulse.peakiness_l, "peakiness_l")
