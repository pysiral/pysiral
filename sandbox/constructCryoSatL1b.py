# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 09:00:30 2015

@author: Stefan
"""

from pysiral.config import ConfigInfo
from construct import (Struct, Array, Padding, Adapter, Bit, BitStruct,
                       SBInt16, UBInt16, SBInt32, UBInt32, SBInt64, UBInt64)

import os
import glob


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


def contruct_cryosat_l1b():

    # Get a CryoSat_l1b File
    info = ConfigInfo()
    l1b_directory = info.local_machine.local_l1b_repository.cryosat2.sar
    l1b_directory = os.path.join(l1b_directory, "2015", "04")
    l1b_files = glob.glob(os.path.join(l1b_directory, "*.DBL"))

    # Just select one
    l1b_file = l1b_files[0]

    time_orbit_group = Struct(
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

    measurement_group = Struct(
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

    corrections_group = Struct(
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

    onehz_waveform_group = Struct(
        "avg_waveform",
        SBInt32("day"),
        UBInt32("sec"),
        UBInt32("msec"),
        OneTenthMicroDeg(SBInt32("latitude")),
        OneTenthMicroDeg(SBInt32("longitude")),
        MilliMeter(SBInt32("altitude_cog")),
        PicoSecond(UBInt64("window_delay")),
        Array(128, UBInt16("wfm")),
        SBInt32("echo_scale_watts"),
        SBInt32("echo_scale_power"),
        UBInt16("num_avg_echoes"),
        UBInt16("flags"))

    waveform_flag_group = BitStruct(
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

    waveform_beam_group = Struct(
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

    waveform_group = Struct(
        "waveform",
        Array(256, UBInt16("wfm")),
        SBInt32("echo_scale_watts"),
        SBInt32("echo_scale_power"),
        UBInt16("num_avg_echoes"),
        waveform_flag_group,
        waveform_beam_group)

    mds_record = Struct(
        "mds_record",
        Array(20, time_orbit_group),
        Array(20, measurement_group),
        corrections_group,
        onehz_waveform_group,
        Array(20, waveform_group))

    mds = Array(128, mds_record)

    print mds.sizeof()

    with open(l1b_file, "rb") as fh:
        fh.seek(14679)
        mds = mds.parse(fh.read(mds.sizeof()))

    print [record.corrections.surface_type for record in mds]


if __name__ == "__main__":
    contruct_cryosat_l1b()
