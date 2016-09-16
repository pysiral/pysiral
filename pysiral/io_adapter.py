# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 15:10:04 2015

@author: Stefan
"""

from pysiral.config import PYSIRAL_VERSION
from pysiral.cryosat2.l1bfile import CryoSatL1B
from pysiral.envisat.sgdrfile import EnvisatSGDR
from pysiral.ers.sgdrfile import ERSSGDR
from pysiral.sentinel3.sral_l1b import Sentinel3SRALL1b

from pysiral.cryosat2.functions import (
    get_tai_datetime_from_timestamp, get_cryosat2_wfm_power,
    get_cryosat2_wfm_range)
from pysiral.classifier import (CS2OCOGParameter, CS2PulsePeakiness,
                                EnvisatWaveformParameter)

from pysiral.clocks import UTCTAIConverter
from pysiral.esa.functions import get_structarr_attr
from pysiral.flag import ORCondition
from pysiral.helper import parse_datetime_str
from pysiral.path import filename_from_path
from pysiral.surface_type import ESA_SURFACE_TYPE_DICT
from pysiral.waveform import (get_waveforms_peak_power, TFMRALeadingEdgeWidth,
                              get_sar_sigma0)

from netCDF4 import num2date
import numpy as np


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
        info.set_attribute("pysiral_version", PYSIRAL_VERSION)

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

        # Convert the TAI timestamp to UTC
        # XXX: Note, the leap seconds are only corrected based on the
        # first timestamp to avoid having identical timestamps.
        # In the unlikely case this will cause problems in orbits that
        # span over a leap seconds change, set check_all=True
        converter = UTCTAIConverter()
        utc_timestamp = converter.tai2utc(tai_timestamp, check_all=False)
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

        # Compute the leading edge width (requires TFMRA retracking)
        wfm = self.l1b.waveform.power
        rng = self.l1b.waveform.range
        radar_mode = self.l1b.waveform.radar_mode
        is_ocean = self.l1b.surface_type.get_by_name("ocean").flag
        lew = TFMRALeadingEdgeWidth(rng, wfm, radar_mode, is_ocean)
        lew1 = lew.get_width_from_thresholds(0.05, 0.5)
        lew2 = lew.get_width_from_thresholds(0.5, 0.95)
        self.l1b.classifier.add(lew1, "leading_edge_width_first_half")
        self.l1b.classifier.add(lew2, "leading_edge_width_second_half")

        # Compute sigma nought
        peak_power = get_waveforms_peak_power(self.l1b.waveform.power)
        tx_power = get_structarr_attr(self.cs2l1b.measurement, "tx_power")
        altitude = self.l1b.time_orbit.altitude
        v = get_structarr_attr(self.cs2l1b.time_orbit, "satellite_velocity")
        vx2, vy2, vz2 = v[:, 0]**2., v[:, 1]**2., v[:, 2]**2
        vx2, vy2, vz2 = vx2.astype(float), vy2.astype(float), vz2.astype(float)
        velocity = np.sqrt(vx2+vy2+vz2)
        sigma0 = get_sar_sigma0(peak_power, tx_power, altitude, velocity)
        self.l1b.classifier.add(sigma0, "sigma0")



class L1bAdapterEnvisat(object):
    """ Converts a Envisat SGDR object into a L1bData object """

    def __init__(self, config):
        self.filename = None
        self._config = config
        self._mission = "envisat"
        self.settings = config.get_mission_settings(self._mission)

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
        self.sgdr = EnvisatSGDR(self.settings)
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
        info.set_attribute("pysiral_version", PYSIRAL_VERSION)
        info.set_attribute("mission", self._mission)
        info.set_attribute("mission_data_version", "final v9.3p5")
        info.set_attribute("orbit", sgdr.mph.abs_orbit)
        info.set_attribute("cycle", sgdr.mph.cycle)
        info.set_attribute("cycle", sgdr.mph.cycle)
        mission_data_source = filename_from_path(sgdr.filename)
        info.set_attribute("mission_data_source", mission_data_source)

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
            self.l1b.correction.set_parameter(
                    correction_name, getattr(mds, correction_name))


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
        self.settings = config.get_mission_settings(mission)

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
        self.sgdr = ERSSGDR(self.settings)
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
        info.set_attribute("pysiral_version", PYSIRAL_VERSION)
        info.set_attribute("mission", self._mission)
        info.set_attribute("mission_data_version", sgdr.nc.software_ver)
        info.set_attribute("orbit", sgdr.nc.abs_orbit)
        info.set_attribute("cycle", sgdr.nc.cycle)
        mission_data_source = filename_from_path(sgdr.nc.filename)
        info.set_attribute("mission_data_source", mission_data_source)

    def _transfer_timeorbit(self):
        """ Extracts the time/orbit data group from the SGDR data """

        # Transfer the orbit position
        self.l1b.time_orbit.set_position(
            self.sgdr.nc.lon_20hz.flatten(),
            self.sgdr.nc.lat_20hz.flatten(),
            self.sgdr.nc.alt_20hz.flatten())

        # Transfer the timestamp
        sgdr_timestamp = self.sgdr.nc.time_20hz.flatten()
        units = self.settings.sgdr_timestamp_units
        calendar = self.settings.sgdr_timestamp_calendar
        timestamp = num2date(sgdr_timestamp, units, calendar)
        self.l1b.time_orbit.timestamp = timestamp

        # Update meta data container
        self.l1b.update_data_limit_attributes()

    def _transfer_waveform_collection(self):
        """ Transfers the waveform data (power & range for each range bin) """

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
        """
        Transfer range correction data from the SGDR netCDF to the
        l1bdata object. The parameter are defined in
        config/mission_def.yaml for ers1/ers2
        -> ersX.settings.sgdr_range_correction_targets

        For a description of the parameter see section 3.10 in the
        REAPER handbook
        """
        grc_dict = self.settings.sgdr_range_correction_targets
        for name in grc_dict.keys():
            target_parameter = grc_dict[name]
            if target_parameter is None:
                continue
            correction = getattr(self.sgdr.nc, target_parameter)
            correction = np.repeat(correction, self.settings.sgdr_n_blocks)
            self.l1b.correction.set_parameter(name, correction)

    def _transfer_classifiers(self):
        """
        Transfer classifier parameter from the SGDR netCDF to the
        l1bdata object. Most parameter are defined in
        config/mission_def.yaml for ers1/ers2
        -> ersX.settings.sgdr_range_correction_targets
        """
        target_dict = self.settings.sgdr_classifier_targets
        for parameter_name in target_dict.keys():
            nc_parameter_name = target_dict[parameter_name]
            nc_parameter = getattr(self.sgdr.nc, nc_parameter_name)
            self.l1b.classifier.add(nc_parameter.flatten(), parameter_name)

        # Add consistent definition of pulse peakiness
        parameter = EnvisatWaveformParameter(self.l1b.waveform.power)
        self.l1b.classifier.add(parameter.pulse_peakiness, "pulse_peakiness")

    def _transfer_surface_type_data(self):
        surface_type = self.sgdr.nc.surface_type
        surface_type = np.repeat(surface_type, self.settings.sgdr_n_blocks)
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
    """ Class for Sentinel3X """

    def __init__(self, config, mission):
        self.filename = None
        self._mission = mission
        self._config = config
        self.settings = config.get_mission_settings(mission)
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
            self.l1b.update_l1b_metadata()

    def _read_sentinel3_sral_l1b(self):
        """ Read the L1b file and create a ERS native L1b object """
        self.sral = Sentinel3SRALL1b()
        self.sral.filename = self.filename
        self.sral.parse_xml_header(self.settings)
        self.sral.parse()
        self.error_status = self.sral.get_status()
        if not self.error_status:
            # TODO: Needs ErrorHandler
            self.sral.post_processing()

    def _transfer_metadata(self):
        """ Extract essential metadata information from SRAL L1 nc file """
        info = self.l1b.info
        product = self.sral.product_info
        info.set_attribute("mission", self._mission)
        info.set_attribute("pysiral_version", PYSIRAL_VERSION)
        mission_data_source = filename_from_path(self.sral.filename)
        info.set_attribute("mission_data_source", mission_data_source)
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
        self.l1b.time_orbit.set_position(
            self.sral.nc.lon_l1b_echo_sar_ku,
            self.sral.nc.lat_l1b_echo_sar_ku,
            self.sral.nc.alt_l1b_echo_sar_ku)

        # Transfer the timestamp
        units = self.settings.time_units
        calendar = self.settings.time_calendar
        seconds = self.sral.nc.time_l1b_echo_sar_ku
        timestamp = num2date(seconds, units, calendar)
        self.l1b.time_orbit.timestamp = timestamp

    def _transfer_waveform_collection(self):
        """ Transfers the waveform data (power & range for each range bin) """

        # self.wfm_power = self.nc.waveform_20_ku
        wfm_power = self.sral.nc.i2q2_meas_ku_l1b_echo_sar_ku
        n_records, n_range_bins = shape = wfm_power.shape

        # Get the window delay
        # "The tracker_range_20hz is the range measured by the onboard tracker
        #  as the window delay, corrected for instrumental effects and
        #  CoG offset"
        tracker_range = self.sral.nc.range_ku_l1b_echo_sar_ku

        # Compute the range for each range bin
        wfm_range = np.ndarray(shape=shape, dtype=np.float32)
        rbw = self.settings.range_bin_width
        ntb = self.settings.nominal_tracking_bin
        rbi = np.arange(n_range_bins)

        # Loop over each waveform
        for i in np.arange(n_records):
            wfm_range[i, :] = tracker_range[i] + (rbi*rbw) - (ntb*rbw)

        # Transfer to l1b object
        radar_mode = self.settings.radar_mode
        self.l1b.waveform.set_waveform_data(wfm_power, wfm_range, radar_mode)

    def _transfer_range_corrections(self):
        """ Retrieve the geophysical range corrections """

        # The definition of which parameter to choose is set in
        # config/mission_def.yaml
        # (see sentinel3x.options.input_dataset.range_correction_target_dict)
        grc_target_dict = self.settings.range_correction_targets
        dummy_val = np.zeros(shape=(self.l1b.n_records), dtype=np.float32)
        for name in grc_target_dict.keys():
            target_parameter = grc_target_dict[name]
            if target_parameter is None:
                correction = dummy_val
            else:
                correction = getattr(self.sral.nc, target_parameter)
            self.l1b.correction.set_parameter(name, correction)

    def _transfer_classifiers(self):
        """
        Transfer the classifiers from L2, L1 (if available) and from
        additional waveform shape analysis
        XXX: This is a makeshift implementation for the expert user
             assessment of early access data
        """

        # Try to get l1 stack parameters (might not be present)
        l1_classifier_target_dict = self.settings.classifier_targets
        for name in l1_classifier_target_dict.keys():
            target_parameter = l1_classifier_target_dict[name]
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
        surface_type = self.sral.nc.surf_type_l1b_echo_sar_ku
        for key in ESA_SURFACE_TYPE_DICT.keys():
            flag = surface_type == ESA_SURFACE_TYPE_DICT[key]
            self.l1b.surface_type.add_flag(flag, key)


class L1bAdapterSentinel3A(L1bAdapterSentinel3):
    """ Class for Sentinel-3A """

    def __init__(self, config):
        super(L1bAdapterSentinel3A, self).__init__(config, "sentinel3a")
