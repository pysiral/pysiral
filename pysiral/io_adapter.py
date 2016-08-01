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
from pysiral.ers.sgdrfile import ERSSGDR
from pysiral.sentinel3.sral_l1b import Sentinel3SRALL1b
from pysiral.helper import parse_datetime_str
from pysiral.path import filename_from_path
from pysiral.classifier import (CS2OCOGParameter, CS2PulsePeakiness,
                                EnvisatWaveformParameter)

from scipy import interpolate
import numpy as np
# import time

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

    def construct_l1b(self, l1b, header_only=False):
        self.l1b = l1b                        # pointer to L1bData object
        self._read_cryosat2l1b_header()       # Read CryoSat-2 L1b header
        if not header_only:
            self.read_msd()

    def read_msd(self):

        # Read the content of the .DLB & .HDR files and store content in
        # this class for
        self._read_cryosat2l1b_data()

        # Transfer relevant metdata
        self._transfer_metadata()

        # Extract time orbit data group (lon, lat, alt, time)
        self._transfer_timeorbit()

        # Extract waveform data (power, range, mode, flags)
        self._transfer_waveform_collection()

        # Extract range correction and tidal information
        self._transfer_range_corrections()

        # Get the surface type data from the L1b file
        # (will serve mainly as source for land mask at this stage)
        self._transfer_surface_type_data()

        # Extract waveform parameters that will be used for
        # surface type classification in the L2 processing
        self._transfer_classifiers()

        # now that all data is transfered to the l1bdata object,
        # complete the attributes in the metadata container
        self.l1b.update_l1b_metadata()

    def _read_cryosat2l1b_header(self):
        """
        Read the header and L1b file and stores information
        on open ocean coverage and geographical location in the metadata
        file
        """

        self.cs2l1b = CryoSatL1B()
        self.cs2l1b.filename = self.filename
        self.cs2l1b.parse_header()

        # Populate metadata object with key information to decide if
        # l1b file contains sea ice information

        # open ocean percent from speficic product header in *.DBL
        self.l1b.info.set_attribute(
            "open_ocean_percent", self.cs2l1b.sph.open_ocean_percent*(10.**-2))

        # Geographical coverage of data from start/stop positions in
        # specific product header in *.DBL
        start_lat = self.cs2l1b.sph.start_lat*(10.**-6)
        stop_lat = self.cs2l1b.sph.stop_lat*(10.**-6)
        start_lon = self.cs2l1b.sph.start_long*(10.**-6)
        stop_lon = self.cs2l1b.sph.stop_long*(10.**-6)
        self.l1b.info.set_attribute("lat_min", np.amin([start_lat, stop_lat]))
        self.l1b.info.set_attribute("lat_max", np.amax([start_lat, stop_lat]))
        self.l1b.info.set_attribute("lon_min", np.amin([start_lon, stop_lon]))
        self.l1b.info.set_attribute("lon_max", np.amax([start_lon, stop_lon]))
        self.l1b.update_region_name()
        error_status = self.cs2l1b.get_status()
        if error_status:
            # TODO: Needs ErrorHandler
            raise IOError()

    def _read_cryosat2l1b_data(self):
        """ Read the L1b file and create a CryoSat-2 native L1b object """
        self.cs2l1b.parse_mds()
        error_status = self.cs2l1b.get_status()
        if error_status:
            # TODO: Needs ErrorHandler
            raise IOError()
        self.cs2l1b.post_processing()

    def _transfer_metadata(self):

        info = self.l1b.info

        # Processing System info
        info.set_attribute("pysiral_version", self._config.PYSIRAL_VERSION)

        # General CryoSat-2 metadata
        info.set_attribute("mission", self._mission)
        info.set_attribute("mission_data_version", self.cs2l1b.baseline)
        info.set_attribute("orbit", self.cs2l1b.sph.abs_orbit_start)
        info.set_attribute("cycle", self.cs2l1b.mph.cycle)
        mission_data_source = filename_from_path(self.cs2l1b.filename)
        info.set_attribute("mission_data_source", mission_data_source)

        # Time-Orbit Metadata
        start_time = parse_datetime_str(self.cs2l1b.sph.start_record_tai_time)
        stop_time = parse_datetime_str(self.cs2l1b.sph.stop_record_tai_time)
        info.set_attribute("start_time", start_time)
        info.set_attribute("stop_time", stop_time)

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
        # XXX: This might need to be switchable
        for i, record in enumerate(self.cs2l1b.waveform):
            echo_power[i, :] = get_cryosat2_wfm_power(
                np.array(record.wfm).astype(np.float32),
                record.linear_scale, record.power_scale)
            echo_range[i, :] = get_cryosat2_wfm_range(
                self.cs2l1b.measurement[i].window_delay, n_range_bins)

            # Transfer to L1bData
        self.l1b.waveform.set_waveform_data(
            echo_power, echo_range, self.cs2l1b.radar_mode)

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
            flag = surface_type == ESA_SURFACE_TYPE_DICT[key]
            self.l1b.surface_type.add_flag(flag, key)

    def _transfer_classifiers(self):

        from pysiral.waveform import get_waveforms_peak_power

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
        # Add the peak power (in Watts)
        # (use l1b waveform power array that is already in physical units)
        peak_power = get_waveforms_peak_power(self.l1b.waveform.power, dB=True)
        self.l1b.classifier.add(peak_power, "peak_power_db")


class L1bAdapterEnvisat(object):
    """ Converts a Envisat SGDR object into a L1bData object """

    def __init__(self, config):
        self.filename = None
        self._config = config
        self._mission = "envisat"

    def construct_l1b(self, l1b):
        """
        Read the Envisat SGDR file and transfers its content to a
        Level1bData instance
        """
        self.l1b = l1b                        # pointer to L1bData object
        # t0 = time.time()
        self._read_envisat_sgdr()             # Read Envisat SGDR data file
        # t1 = time.time()
        # print "Parsing Envisat SGDR file in %.3g seconds" % (t1 - t0)
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
        """ Transfers the waveform data (power & range for each range bin) """
        from pysiral.flag import ANDCondition
        # Transfer the reformed 18Hz waveforms
        self.l1b.waveform.set_waveform_data(
            self.sgdr.mds_18hz.power,
            self.sgdr.mds_18hz.range,
            self.sgdr.radar_mode)
        # This is from SICCI-1 processor, check of mcd flags
        valid = ANDCondition()
        valid.add(np.logical_not(self.sgdr.mds_18hz.flag_packet_length_error))
        valid.add(np.logical_not(self.sgdr.mds_18hz.flag_obdh_invalid))
        valid.add(np.logical_not(self.sgdr.mds_18hz.flag_agc_fault))
        valid.add(np.logical_not(self.sgdr.mds_18hz.flag_rx_delay_fault))
        valid.add(np.logical_not(self.sgdr.mds_18hz.flag_waveform_fault))
        # ku_chirp_band_id (0 -> 320Mhz)
        valid.add(self.sgdr.mds_18hz.ku_chirp_band_id == 0)
        self.l1b.waveform.set_valid_flag(valid.flag)

    def _transfer_range_corrections(self):
        # Transfer all the correction in the list
        mds = self.sgdr.mds_18hz
        for correction_name in mds.sgdr_geophysical_correction_list:
            if correction_name not in self._config.parameter.correction_list:
                continue
            self.l1b.correction.set_parameter(
                    correction_name, getattr(mds, correction_name))
        # Envisat specific: There are several options for sources
        #   of geophysical correction in the SGDR files. Selct those
        #   specified in the mission defaults
        #   (see config/mission_def.yaml)
        mission_defaults = self._config.get_mission_options(self._mission)
        correction_options = mission_defaults.geophysical_corrections
        for option in correction_options.iterbranches():
            self.l1b.correction.set_parameter(
                    option.target, getattr(mds, option.selection))

    def _transfer_surface_type_data(self):
        surface_type = self.sgdr.mds_18hz.surface_type
        for key in ESA_SURFACE_TYPE_DICT.keys():
            flag = surface_type == ESA_SURFACE_TYPE_DICT[key]
            self.l1b.surface_type.add_flag(flag, key)

    def _transfer_classifiers(self):
        """
        OCOC retracker parameters are needed for surface type classification
        in Envisat L2 processing
        """
        wfm = self.sgdr.mds_18hz.power
        parameter = EnvisatWaveformParameter(wfm)
        self.l1b.classifier.add(parameter.pulse_peakiness, "pulse_peakiness")
        sea_ice_backscatter = self.sgdr.mds_18hz.sea_ice_backscatter
        self.l1b.classifier.add(sea_ice_backscatter, "sea_ice_backscatter")


class L1bAdapterERS(object):
    """ Converts a Envisat SGDR object into a L1bData object """

    def __init__(self, config, mission):
        self.filename = None
        self._config = config
        self._mission = mission

    def construct_l1b(self, l1b):
        """
        Read the Envisat SGDR file and transfers its content to a
        Level1bData instance
        """
        self.l1b = l1b                        # pointer to L1bData object
        # t0 = time.time()
        self._read_ers_sgdr()                 # Read Envisat SGDR data file
        # t1 = time.time()
        # print "Parsing Envisat SGDR file in %.3g seconds" % (t1 - t0)
        self._transfer_metadata()             # (orbit, radar mode, ..)
        self._transfer_timeorbit()            # (lon, lat, alt, time)
        self._transfer_waveform_collection()  # (power, range)
        self._transfer_range_corrections()    # (range corrections)
        self._transfer_surface_type_data()    # (land flag, ocean flag, ...)
        self._transfer_classifiers()          # (beam parameters, flags, ...)

    def _read_ers_sgdr(self):
        """ Read the L1b file and create a ERS native L1b object """
        self.sgdr = ERSSGDR()
        self.sgdr.filename = self.filename
        self.sgdr.parse()
        error_status = self.sgdr.get_status()
        if error_status:
            # TODO: Needs ErrorHandler
            raise IOError()
        self.sgdr.post_processing()

    def _transfer_metadata(self):
        pass
        """ Extract essential metadata information from SGDR file """
        info = self.l1b.info
        sgdr = self.sgdr
        info.set_attribute("mission", self._mission)
        info.set_attribute("mission_data_version", sgdr.nc.software_ver)
        info.set_attribute("orbit", sgdr.nc.abs_orbit)
        info.set_attribute("cycle", sgdr.nc.cycle)
        info.set_attribute("mission_data_source", sgdr.nc.proc_centre)

    def _transfer_timeorbit(self):
        """ Extracts the time/orbit data group from the SGDR data """
        from netCDF4 import num2date
        # Transfer the orbit position
        self.l1b.time_orbit.set_position(
            self.sgdr.nc.lon_20hz.flatten(),
            self.sgdr.nc.lat_20hz.flatten(),
            self.sgdr.nc.alt_20hz.flatten())
        # Transfer the timestamp
        units = "seconds since 1990-01-01 00:00:00"
        calendar = "gregorian"
        timestamp = num2date(
            self.sgdr.nc.time_20hz.flatten(), units, calendar)
        self.l1b.time_orbit.timestamp = timestamp
        # Update meta data container
        self.l1b.update_data_limit_attributes()

    def _transfer_waveform_collection(self):
        """ Transfers the waveform data (power & range for each range bin) """

        from pysiral.flag import ORCondition

        # Transfer the reformed 18Hz waveforms
        self.l1b.waveform.set_waveform_data(
            self.sgdr.wfm_power,
            self.sgdr.wfm_range,
            self.sgdr.radar_mode)

        # Set valid flag to exclude calibration data
        # (see section 3.5 of Reaper handbook)
        tracking_state = self.sgdr.nc.alt_state_flag_20hz.flatten()
        valid = ORCondition()
        valid.add(tracking_state == 2)
        valid.add(tracking_state == 3)
        self.l1b.waveform.set_valid_flag(valid.flag)

    def _transfer_range_corrections(self):

        # (see section 3.10 in REAPER handbook)
        # TODO: move selection dict to configuration files
        range_correction_target_dict = {
            "dry_troposphere": "model_dry_tropo_corr",
            "wet_troposphere": "model_wet_tropo_corr",
            "inverse_barometric": "inv_bar_corr",
            "dynamic_atmosphere": "hf_fluctuations_corr",
            "ionosphere": "iono_corr_model",
            "ocean_tide_elastic": "ocean_tide_sol1",
            "ocean_tide_long_period": "ocean_tide_equil",
            "ocean_loading_tide": "load_tide_sol1",
            "solid_earth_tide": "solid_earth_tide",
            "geocentric_polar_tide": "pole_tide",
            "total_geocentric_ocean_tide": None}

        for name in range_correction_target_dict.keys():
            target_parameter = range_correction_target_dict[name]
            if target_parameter is None:
                continue
            correction = getattr(self.sgdr.nc, target_parameter)
            correction = np.repeat(correction, self.sgdr.n_blocks)
            self.l1b.correction.set_parameter(name, correction)

    def _transfer_classifiers(self):
        chirp_type = self.sgdr.nc.alt_state_flag_chirp_type_20hz.flatten()
        self.l1b.classifier.add(chirp_type, "chirp_type")
        parameter = EnvisatWaveformParameter(self.l1b.waveform.power)
        self.l1b.classifier.add(parameter.pulse_peakiness, "pulse_peakiness")

    def _transfer_surface_type_data(self):
        surface_type = self.sgdr.nc.surface_type
        surface_type = np.repeat(surface_type, self.sgdr.n_blocks)
        for key in ESA_SURFACE_TYPE_DICT.keys():
            flag = surface_type == ESA_SURFACE_TYPE_DICT[key]
            self.l1b.surface_type.add_flag(flag, key)


class L1bAdapterERS1(L1bAdapterERS):
    """ Class for ERS-1 """

    def __init__(self, config):
        super(L1bAdapterERS1, self).__init__(config, "ers1")


class L1bAdapterERS2(L1bAdapterERS):
    """ Class for ERS-1 """

    def __init__(self, config):
        super(L1bAdapterERS2, self).__init__(config, "ers2")


class L1bAdapterSentinel3(object):
    """ Class for Sentinel3 """
    """ Converts a Envisat SGDR object into a L1bData object """

    def __init__(self, config, mission):
        self.filename = None
        self._mission = mission
        self._config = config
        self.error_status = False

    def construct_l1b(self, l1b, header_only=False):
        """
        Read the Envisat SGDR file and transfers its content to a
        Level1bData instance
        """
        self.l1b = l1b
        self._read_sentinel3_sral_l1b()
        self._transfer_metadata()
        self._test_ku_data_present()
        if not header_only and not self.l1b.error_status:
            self._transfer_timeorbit()
            self._transfer_waveform_collection()
            self._transfer_range_corrections()
            self._transfer_surface_type_data()
            self._transfer_classifiers()

    def _read_sentinel3_sral_l1b(self):
        """ Read the L1b file and create a ERS native L1b object """
        self.sral = Sentinel3SRALL1b()
        self.sral.filename = self.filename
        self.sral.parse()
        self.error_status = self.sral.get_status()
        if not self.error_status:
            # TODO: Needs ErrorHandler
            self.sral.post_processing()

    def _transfer_metadata(self):
        pass
        """ Extract essential metadata information from SGDR file """
        info = self.l1b.info
        product = self.sral.product_info
        info.set_attribute("mission", self._mission)
        info.set_attribute("mission_data_source", self.sral.nc.institution)
        info.set_attribute("sar_mode_percent", product.sar_mode_percentage)
        info.set_attribute("open_ocean_percent", product.open_ocean_percentage)

    def _test_ku_data_present(self):
        if len(self.sral.nc.time_l1b_echo_sar_ku) == 0:
            self.error_status = True
            self.l1b.error_status = True

    def _transfer_timeorbit(self):
        """ Extracts the time/orbit data group from the SGDR data """
        from netCDF4 import num2date
        # Transfer the orbit position
#        self.l1b.time_orbit.set_position(
#            self.sral.nc.lon_20_ku,
#            self.sral.nc.lat_20_ku,
#            self.sral.nc.alt_20_ku)
        self.l1b.time_orbit.set_position(
            self.sral.nc.lon_l1b_echo_sar_ku,
            self.sral.nc.lat_l1b_echo_sar_ku,
            self.sral.nc.alt_l1b_echo_sar_ku)
        # Transfer the timestamp
        units = "seconds since 2000-01-01 00:00:00.0"
        calendar = "gregorian"
        timestamp = num2date(
                self.sral.nc.time_l1b_echo_sar_ku, units, calendar)

        self.l1b.time_orbit.timestamp = timestamp
        # Update meta data container
        self.l1b.update_data_limit_attributes()

    def _transfer_waveform_collection(self):
        """ Transfers the waveform data (power & range for each range bin) """

        # Transfer the 20Hz waveforms
        self.l1b.waveform.set_waveform_data(
            self.sral.wfm_power,
            self.sral.wfm_range,
            self.sral.radar_mode)

    def _transfer_range_corrections(self):

        # (see section 3.10 in REAPER handbook)
        # TODO: move selection dict to configuration files
#        range_correction_target_dict = {
#            "dry_troposphere": "mod_dry_tropo_cor_meas_altitude_01",
#            "wet_troposphere": "mod_wet_tropo_cor_meas_altitude_01",
#            "inverse_barometric": "iono_cor_alt_20_ku",
#            "dynamic_atmosphere": "hf_fluct_cor_01",
#            "ionosphere": "inv_bar_cor_01",
#            "ocean_tide_elastic": None,
#            "ocean_tide_long_period": "ocean_tide_eq_01",
#            "ocean_loading_tide": "load_tide_sol1_01",
#            "solid_earth_tide": "solid_earth_tide_01",
#            "geocentric_polar_tide": "pole_tide_01",
#            "total_geocentric_ocean_tide": None}

        range_correction_target_dict = {
            "dry_troposphere": None,
            "wet_troposphere": None,
            "inverse_barometric": None,
            "dynamic_atmosphere": None,
            "ionosphere": None,
            "ocean_tide_elastic": None,
            "ocean_tide_long_period": None,
            "ocean_loading_tide": None,
            "solid_earth_tide": None,
            "geocentric_polar_tide": None,
            "total_geocentric_ocean_tide": None}

        dummy_val = np.zeros(shape=(self.l1b.n_records), dtype=np.float32)
        for name in range_correction_target_dict.keys():
            target_parameter = range_correction_target_dict[name]
            if target_parameter is None:
                correction = dummy_val
            else:
                correction = getattr(self.sral.nc, target_parameter)

#            # Some corrections are only available in 1Hz, while 20Hz is needed
#            if len(correction) != self.l1b.n_records:
#                corr_interp = interpolate.interp1d(
#                    self.sral.nc.time_01, correction, bounds_error=False)
#                correction = corr_interp(self.sral.nc.time_20_ku)

            self.l1b.correction.set_parameter(name, correction)

    def _transfer_classifiers(self):
        """
        Transfer the classifiers from L2, L1 (if available) and from
        additional waveform shape analysis
        XXX: This is a makeshift implementation for the expert user
             assessment of early access data
        """

        # Try to get l1 stack parameters (might not be present)
        # TODO: Move this into the mission config structure
        l1_classifier_target_dict = {
            "stack_standard_deviation": "stdev_stack_l1b_echo_sar_ku",
            "stack_skewness": "skew_stack_l1b_echo_sar_ku",
            "stack_kurtosis": "kurt_stack_l1b_echo_sar_ku",
            "stack_maximum_power": "max_stack_l1b_echo_sar_ku",
            "num_stack_waveforms": "nb_stack_l1b_echo_sar_ku"}

        for name in l1_classifier_target_dict.keys():
            target_parameter = l1_classifier_target_dict[name]
            classifier = getattr(self.sral.nc, target_parameter)
            self.l1b.classifier.add(classifier, name)

#        if hasattr(self.sral, "l1nc"):
#            time_l2 = self.sral.nc.time_20_ku
#            time_l1 = self.sral.l1nc.time_l1b_echo_sar_ku
#            for name in l1_classifier_target_dict.keys():
#                target_parameter = l1_classifier_target_dict[name]
#                classifier = getattr(self.sral.l1nc, target_parameter)
#                # interpolation is required since L1 and L2 records are
#                # not identical
#                interp_l1_l2 = interpolate.interp1d(
#                    time_l1, classifier, bounds_error=False)
#                classifier = interp_l1_l2(time_l2)
#                self.l1b.classifier.add(classifier, name)

        # Transfer L2 classifier parameter
        # XXX: Currently None
        l2_classifier_target_dict = {}
        for name in l2_classifier_target_dict.keys():
            target_parameter = l2_classifier_target_dict[name]
            classifier = getattr(self.sral.nc, target_parameter)
            self.l1b.classifier.add(classifier, name)

        # Compute sar specific waveform classifiers after Ricker et al. 2014
        wfm = self.l1b.waveform.power
        ocog = CS2OCOGParameter(wfm)
        self.l1b.classifier.add(ocog.width, "ocog_width")
        self.l1b.classifier.add(ocog.amplitude, "ocog_amplitude")
        # Calculate the Peakiness (CryoSat-2 notation)
        pulse = CS2PulsePeakiness(wfm)
        self.l1b.classifier.add(pulse.peakiness, "peakiness")
        self.l1b.classifier.add(pulse.peakiness_r, "peakiness_r")
        self.l1b.classifier.add(pulse.peakiness_l, "peakiness_l")

    def _transfer_surface_type_data(self):
        # surface_type = self.sral.nc.surf_type_20_ku
        surface_type = self.sral.nc.surf_type_l1b_echo_sar_ku
        # surface_type = np.repeat(surface_type, self.sgdr.n_blocks)
        for key in ESA_SURFACE_TYPE_DICT.keys():
            flag = surface_type == ESA_SURFACE_TYPE_DICT[key]
            self.l1b.surface_type.add_flag(flag, key)


class L1bAdapterSentinel3A(L1bAdapterSentinel3):
    """ Class for ERS-1 """

    def __init__(self, config):
        super(L1bAdapterSentinel3A, self).__init__(config, "sentinel3a")
