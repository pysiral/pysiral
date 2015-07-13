# -*- coding=utf-8 -*-
"""
Created on Fri Jul 10 21:00:37 2015

@author=Stefan
"""

from construct import (Struct, Array, Padding, Adapter, Bit, BitStruct,
                       SBInt16, UBInt16, SBInt32, UBInt32, SBInt64, UBInt64)


class OneHundredth(Adapter):
    def _decode(self, obj, context):
        return float(obj)/100


class OneHundredthDecibel(Adapter):
    def _decode(self, obj, context):
        return float(obj)/100


class OneTenthMicroDeg(Adapter):
    def _decode(self, obj, context):
        return float(obj)*1e-7


class MilliMeter(Adapter):
    def _decode(self, obj, context):
        return float(obj)*1e-3


class Micrometer(Adapter):
    def _decode(self, obj, context):
        return float(obj)*1e-6


class PicoSecond(Adapter):
    def _decode(self, obj, context):
        return float(obj)*1e-12


class MicroWatts(Adapter):
    def _decode(self, obj, context):
        return float(obj)*1e-6


class MicroRadians(Adapter):
    def _decode(self, obj, context):
        return float(obj)*1e-6


class OneTenthMicroRadians(Adapter):
    def _decode(self, obj, context):
        return float(obj)*1e-6


class Cryosat2L1bMDSDefinition(object):
    """ Holds the Definition of L1b MDS records """

    def __init__(self):

        self.baseline = None
        self.radar_mode = None
        self.n_records = 0

    def get_mds_parser(self):

        self.time_orbit_group = Struct(
            "time_orbit",
            SBInt32("day"),
            UBInt32("sec"),
            UBInt32("msec"),
            UBInt16("uso_corr"),
            UBInt32("mode_id"),
            UBInt16("source_sequence_counter"),
            UBInt32("instrument_configuration"),
            UBInt32("burst_counter"),
            OneTenthMicroDeg(SBInt32("atitude")),
            OneTenthMicroDeg(SBInt32("longitude")),
            MilliMeter(SBInt32("altitude_cog")),
            MilliMeter(SBInt32("altitude_rate")),
            Array(3, MilliMeter(SBInt32("satellite_velocity"))),
            Array(3, Micrometer(SBInt32("real_beam"))),
            Array(3, Micrometer(SBInt32("interferometry_baseline"))),
            UBInt16("star_tracker_usage"),
            OneTenthMicroDeg(SBInt32("antenna_roll")),
            OneTenthMicroDeg(SBInt32("antenna_pitch")),
            OneTenthMicroDeg(SBInt32("antenna_yaw")),
            UBInt32("fbr_confidence_flag"),
            Padding(4))

        self.measurement_group = Struct(
            "measurement",
            PicoSecond(SBInt64("window_delay")),
            SBInt32("h0"),
            SBInt32("cor2"),
            SBInt32("lai"),
            SBInt32("fai"),
            OneHundredthDecibel(SBInt32("agc_1")),
            OneHundredthDecibel(SBInt32("agc_2")),
            OneHundredthDecibel(SBInt32("fixed_gain_1")),
            OneHundredthDecibel(SBInt32("fixed_gain_2")),
            MicroWatts(SBInt32("tx_power")),
            MilliMeter(SBInt32("doppler_range_correction")),
            MilliMeter(SBInt32("instrument_range_correction_txrx")),
            MilliMeter(SBInt32("instrument_range_correction_rx")),
            OneHundredthDecibel(SBInt32("instrument_gain_correction_txrx")),
            OneHundredthDecibel(SBInt32("instrument_gain_correction_rx")),
            MicroRadians(SBInt32("phase_correction_internal")),
            MicroRadians(SBInt32("phase_correction_external")),
            OneHundredthDecibel(SBInt32("noise_power")),
            MicroRadians(SBInt32("phase_slope_correction")),
            Padding(4))

        self.corrections_group = Struct(
            "corrections",
            MilliMeter(SBInt32("dry_troposphere")),
            MilliMeter(SBInt32("wet_troposphere")),
            MilliMeter(SBInt32("inverse_barometric")),
            MilliMeter(SBInt32("dynamic_atmosphere")),
            MilliMeter(SBInt32("ionosphere_gim")),
            MilliMeter(SBInt32("ionosphere_mod")),
            MilliMeter(SBInt32("ocean_tide_elastic")),
            MilliMeter(SBInt32("ocean_tide_long_period")),
            MilliMeter(SBInt32("ocean_loading_tide")),
            MilliMeter(SBInt32("solid_earth_tide")),
            MilliMeter(SBInt32("geocentric_polar_tide")),
            UBInt32("surface_type"),
            Padding(4),
            UBInt32("correction_status"),
            UBInt32("correction_error"),
            Padding(4))

        self.onehz_waveform_group = Struct(
            "avg_waveform",
            SBInt32("day"),
            UBInt32("sec"),
            UBInt32("msec"),
            OneTenthMicroDeg(SBInt32("latitude")),
            OneTenthMicroDeg(SBInt32("longitude")),
            MilliMeter(SBInt32("altitude_cog")),
            PicoSecond(UBInt64("window_delay")),
            Array(128, UBInt16("wfm")),
            SBInt32("linear_scale"),
            SBInt32("power_scale"),
            UBInt16("num_avg_echoes"),
            UBInt16("flags"))

        self.waveform_flag_group = BitStruct(
            "flag",
            Bit("approximate_beam_steering"),
            Bit("exact_beam_steering"),
            Bit("doppler_weighting_computed"),
            Bit("dopple_weighting_applied_before_stack"),
            Bit("multi_look_incomplete"),
            Bit("beam_angle_steering_error"),
            Bit("anti_aliased_power_echoes"),
            Bit("auto_beam_steering"),
            Padding(8))

        self.waveform_beam_group = Struct(
            "beam",
            OneHundredth(UBInt16("stack_standard_deviation")),
            OneHundredth(UBInt16("stack_centre")),
            OneHundredthDecibel(SBInt16("stack_scaled_amplitude")),
            OneHundredth(SBInt16("stack_skewness")),
            OneHundredth(SBInt16("stack_kurtosis")),
            MicroRadians(SBInt16("stack_standard_deviation_boresight")),
            MicroRadians(SBInt16("stack_centre_angle")),
            OneTenthMicroRadians(SBInt32("doppler_angle_start")),
            OneTenthMicroRadians(SBInt32("doppler_angle_stop")),
            OneTenthMicroRadians(SBInt32("look_angle_start")),
            OneTenthMicroRadians(SBInt32("look_angle_stop")),
            UBInt16("n_beams_after_weighing"),
            UBInt16("n_beams_before_weighing"),
            Padding(66))

        self.waveform_group = Struct(
            "waveform",
            Array(256, UBInt16("wfm")),
            SBInt32("linear_scale"),
            SBInt32("power_scale"),
            UBInt16("num_avg_echoes"),
            self.waveform_flag_group,
            self.waveform_beam_group)

        self.mds_record = Struct(
            "mds_record",
            Array(20, self.time_orbit_group),
            Array(20, self.measurement_group),
            self.corrections_group,
            self.onehz_waveform_group,
            Array(20, self.waveform_group))

        self.mds = Array(self.n_records, self.mds_record)
        return self.mds


def cryosat2_get_mds_def(radar_mode, baseline, n_records):
    # XXX: Testing purposes only, needs functionality

    definition = Cryosat2L1bMDSDefinition()
    definition.baseline = baseline
    definition.radar_mode = radar_mode
    definition.n_records = n_records

    return definition.get_mds_parser()
