# -*- coding: utf-8 -*-

from pysiral.units import (OneHundredth, OneHundredthDecibel, MilliMeter,
                           MicroDeg, TenMicroDeg, Centimeter, TenThousands,
                           TenPascal, OneTenths, OneThousands, Per256, Per2048,
                           Per8096)

from construct import (Struct, Array, Padding, Bit, BitStruct, SBInt8, UBInt8,
                       SBInt16, UBInt16, SBInt32, UBInt32)


class EnvisatSGDRMDS(object):

    def __init__(self):
        self.n_records = 0
        self.n_blocks = 20


class EnvisatSGDRMDSRA2(EnvisatSGDRMDS):

    def __init__(self):
        super(EnvisatSGDRMDSRA2, self).__init__()

    def get_mds_parser(self):

        self.timestamp = Struct(
            "utc_timestamp",
            SBInt32("day"),
            UBInt32("sec"),
            UBInt32("msec"))

        self.time_orbit = Struct(
            "time_orbit",
            self.timestamp,
            SBInt8("quality_indicator"),
            Padding(3),
            MicroDeg(SBInt32("latitude")),
            MicroDeg(SBInt32("longitude")),
            UBInt32("source_packet_counter"),
            UBInt32("instrument_mode_id"),
            UBInt32("measurement_confidence_data"),
            MilliMeter(UBInt32("altitude")),
            Array(20, MilliMeter(SBInt16("18hz_altitude_differences"))),
            MilliMeter(SBInt16("instantaneous_altitude_rate")),
            Padding(50))

        self.range_information = Struct(
            "range_information",
            Array(20, MilliMeter(UBInt32("18hz_tracker_range_no_doppler_ku"))),
            Array(20, MilliMeter(UBInt32("18hz_tracker_range_no_doppler_s"))),
            UBInt32("map_valid_points_tracker_ku"),
            Padding(4),
            MilliMeter(UBInt32("ocean_range_ku")),
            MilliMeter(UBInt32("band_ocean_range_s")),
            Array(20, MilliMeter(UBInt32("18hz_ocean_range_ku"))),
            Array(20, MilliMeter(UBInt32("18hz_ocean_range_s"))),
            MilliMeter(UBInt16("sd_18hz_ocean_range_ku")),
            MilliMeter(UBInt16("sd_18hz_ocean_range_s")),
            UBInt16("18hz_n_valid_ocean_ku"),
            UBInt16("18hz_n_valid_ocean_s"),
            BitStruct("18hz_valid_points_ocean_ku", Array(32, Bit("flag"))),
            BitStruct("18hz_valid_points_ocean_s", Array(32, Bit("flag"))),
            Array(20, MilliMeter(UBInt32("18Hz_ice1_range_ku"))),
            Array(20, MilliMeter(UBInt32("18Hz_ice1_range_s"))),
            Array(20, MilliMeter(UBInt32("18Hz_ice2_range_ku"))),
            Array(20, MilliMeter(UBInt32("18Hz_ice2_range_s"))),
            Array(20, MilliMeter(UBInt32("18Hz_sea_ice_range_ku"))),
            Array(20, TenMicroDeg(SBInt16("18hz_latitude_differences"))),
            Array(20, TenMicroDeg(SBInt16("18hz_longitude_differences"))))

        self.range_correction = Struct(
            "range_correction",
            Array(20, MilliMeter(SBInt16("18Hz_ku_range_instrumental"))),
            Array(20, MilliMeter(SBInt16("18Hz_s_range_instrumental"))),
            Array(20, MilliMeter(SBInt16("18Hz_ku_range_doppler"))),
            Array(20, MilliMeter(SBInt16("18Hz_s_range_doppler"))),
            Array(20, MilliMeter(SBInt16("18Hz_ku_range_doppler_slope"))),
            Array(20, MilliMeter(SBInt16("18Hz_s_range_doppler_slope"))),
            MilliMeter(SBInt16("dry_troposphere")),
            MilliMeter(SBInt16("inverse_barometric")),
            MilliMeter(SBInt16("wet_tropospheric_model")),
            MilliMeter(SBInt16("wet_tropospheric_mwr")),
            MilliMeter(SBInt16("ra2_ionospheric_ku")),
            MilliMeter(SBInt16("ra2_ionospheric_s")),
            MilliMeter(SBInt16("doris_ku_ionospheric")),
            MilliMeter(SBInt16("doris_s_ionospheric")),
            MilliMeter(SBInt16("model_ku_ionospheric")),
            MilliMeter(SBInt16("model_s_ionospheric")),
            MilliMeter(SBInt16("sea_state_bias_ku")),
            MilliMeter(SBInt16("sea_state_bias_s")),
            MilliMeter(SBInt16("dib_hf")),
            Padding(10))

        self.significant_wave_height = Struct(
            "significant_wave_height",
            MilliMeter(SBInt32("swh_square_ku")),
            MilliMeter(SBInt32("swh_square_s")),
            MilliMeter(SBInt16("swh_ku")),
            MilliMeter(SBInt16("swh_s")),
            MilliMeter(SBInt16("swh_sd_18hz_ku")),
            MilliMeter(SBInt16("swh_sd_18hz_s")),
            UBInt16("shw_18hz_n_valid_ku"),
            UBInt16("shw_18hz_n_valid_s"),
            BitStruct("slope_model_present", Array(32, Bit("flag"))),
            Centimeter(SBInt32("1hz_elevation_echo_point")),
            Array(20, Centimeter(SBInt16("18hz_elevation_difference"))),
            Array(20, TenMicroDeg(SBInt16("18hz_slope_latitude_diff"))),
            Array(20, TenMicroDeg(SBInt16("18hz_slope_longitude_diff"))),
            Array(20, MilliMeter(SBInt16("18hz_ice2_leading_edge_width_ku"))),
            Array(20, MilliMeter(SBInt16("18hz_ice2_leading_edge_width_s"))),
            Padding(40))

        self.backscatter = Struct(
            "backscatter",
            Array(20, OneHundredthDecibel(SBInt16("18hz_ku_k_cal_ku"))),
            Array(20, OneHundredthDecibel(SBInt16("18hz_s_k_cal_s"))),
            BitStruct("map_valid_k_cal_ku", Array(32, Bit("flag"))),
            Padding(4),
            OneHundredthDecibel(SBInt16("ocean_sigma_corr_ku")),
            OneHundredthDecibel(SBInt16("ocean_sigma_corr_s")),
            OneHundredthDecibel(SBInt16("sd_18hz_cean_sigma_ku")),
            OneHundredthDecibel(SBInt16("sd_18hz_ocean_sigma_s")),
            UBInt16("n_valid_18hz_ocean_sigma_ku"),
            UBInt16("n_valid_18hz_ocean_sigma_s"),
            Array(20, OneHundredthDecibel(SBInt16("18hz_ice1_sigma_ku"))),
            Array(20, OneHundredthDecibel(SBInt16("18hz_ice1_sigma_s"))),
            Array(20, OneHundredthDecibel(SBInt16("18hz_ice2_le_sigma_ku"))),
            Array(20, OneHundredthDecibel(SBInt16("18hz_ice2_le_sigma_s"))),
            Array(20, OneHundredthDecibel(SBInt16("18hz_ice2_sigma_ku"))),
            Array(20, OneHundredthDecibel(SBInt16("18hz_ice2_sigma_s"))),
            Array(20, OneHundredthDecibel(SBInt16("18hz_sea_ice_sigma_ku"))),
            Padding(40))

        self.backscatter_correction = Struct(
            "backscatter_correction",
            OneHundredthDecibel(SBInt16("agc_net_instrumental_ku")),
            OneHundredthDecibel(SBInt16("agc_net_instrumental_s")),
            OneHundredthDecibel(SBInt16("atmospheric_attenuation_ku")),
            OneHundredthDecibel(SBInt16("atmospheric_attenuation_s")),
            OneHundredthDecibel(SBInt32("rain_attenuation_ku")))

        self.off_nadir_information = Struct(
            "off_nadir_angle",
            TenThousands(SBInt16("square_angle_from_platform_data")),
            TenThousands(SBInt16("square_angle_from_waveform_data")),
            Array(20, SBInt32("slope_ice2_first_trailing_edge_ku")),
            Array(20, SBInt32("slope_ice2_first_trailing_edge_s")),
            Array(20, SBInt32("slope_ice2_second_trailing_edge_ku")),
            Array(20, SBInt32("slope_ice2_second_trailing_edge_s")),
            Padding(40))

        self.geophysical_information = Struct(
            "geophysical_information",
            MilliMeter(SBInt32("mean_sea_surface_height")),
            MilliMeter(SBInt32("geoid_height")),
            MilliMeter(SBInt32("ocean_depth_land_elevation")),
            MilliMeter(SBInt16("total_geopcentric_ocean_tide_1")),
            MilliMeter(SBInt16("total_geopcentric_ocean_tide_2")),
            MilliMeter(SBInt16("long_period_tide")),
            MilliMeter(SBInt16("tide_loading_2")),
            MilliMeter(SBInt16("solid_earth_tide")),
            MilliMeter(SBInt16("geocentric_pole_tide")),
            TenPascal(SBInt16("model_surface_atmospheric_pressure")),
            OneHundredth(SBInt16("mwr_water_vapour_content")),
            OneHundredth(SBInt16("mwr_liquid_water_content")),
            OneTenths(SBInt16("ra2_total_electron_content")),
            MilliMeter(SBInt16("ra2_wind_speed")),
            MilliMeter(SBInt16("model_wind_vector_u")),
            MilliMeter(SBInt16("model_wind_vector_v")),
            MilliMeter(SBInt16("tide_loading_1")),
            Padding(8))

        self.mwr_information = Struct(
            "mwr_information",
            OneHundredth(SBInt16("mwr_brightness_temp_23_8ghz")),
            OneHundredth(SBInt16("mwr_brightness_temp_36_5ghz")),
            OneHundredth(SBInt16("mwr_sd_brightness_temp_23_8ghz")),
            OneHundredth(SBInt16("mwr_sd_brightness_temp_36_5ghz")),
            Padding(2))

        self.flags = Struct(
            "flag",
            UBInt16("average_ku_chirp_band"),
            BitStruct("ku_chirp_band_id", Array(64, Bit("flag"))),
            BitStruct("error_chirp_band_id", Array(32, Bit("flag"))),
            BitStruct("instrument", Array(32, Bit("flag"))),
            BitStruct("fault_identifier", Array(64, Bit("flag"))),
            Padding(8),
            BitStruct("waveform_fault_identifier", Array(64, Bit("flag"))),
            BitStruct("instrument_mode_id", Array(96, Bit("flag"))),
            UBInt16("n_measures_flight_calibration_s"),
            UBInt16("n_measures_flight_calibration_ku"),
            BitStruct("mwr_instrument_quality", Array(16, Bit("flag"))),
            Padding(6),
            Padding(8),
            Padding(8),
            BitStruct("ocean_retracking_quality_ku", Array(32, Bit("flag"))),
            BitStruct("ocean_retracking_quality_s", Array(32, Bit("flag"))),
            BitStruct("ice1_retracking_quality_ku", Array(32, Bit("flag"))),
            BitStruct("ice1_retracking_quality_s", Array(32, Bit("flag"))),
            BitStruct("ice2_retracking_quality_ku", Array(32, Bit("flag"))),
            BitStruct("ice2_retracking_quality_s", Array(32, Bit("flag"))),
            BitStruct("sea_ice_retracking_quality_ku", Array(32, Bit("flag"))),
            OneThousands(UBInt16("1hz_pulse_peakiness_ku")),
            OneThousands(UBInt16("1hz_pulse_peakiness_s")),
            UBInt16("altimeter_surface_type"),
            UBInt16("radiometer_land_ocean"),
            UBInt16("mwr_quality_interpolation"),
            UBInt16("altimeter_rain"),
            UBInt16("interpolation"),
            UBInt8("sea_ice"),
            BitStruct("membership_01", Array(8, Bit("flag"))),
            BitStruct("membership_02", Array(8, Bit("flag"))),
            BitStruct("membership_03", Array(8, Bit("flag"))),
            BitStruct("membership_04", Array(8, Bit("flag"))),
            Padding(1))

        self.mds_record = Struct(
            "mds_record",
            self.time_orbit,
            self.range_information,
            self.range_correction,
            self.significant_wave_height,
            self.backscatter,
            self.backscatter_correction,
            self.off_nadir_information,
            self.geophysical_information,
            self.mwr_information,
            self.flags)

        self.mds = Array(self.n_records, self.mds_record)

        return self.mds


class EnvisatSGDRMDSWFM18HZ(EnvisatSGDRMDS):

    def __init__(self):
        super(EnvisatSGDRMDSWFM18HZ, self).__init__()

    def get_mds_parser(self):

        self.timestamp = Struct(
            "utc_timestamp",
            SBInt32("day"),
            UBInt32("sec"),
            UBInt32("msec"))

        self.waveform_data = Struct(
            "wfm",
            Array(128, Per2048(UBInt16("average_wfm_if_corr_ku"))),
            Array(2, Per2048(UBInt16("central_filters_if_corr_ku"))),
            Array(64, Per8096(UBInt16("average_wfm_if_corr_s"))),
            Array(2, SBInt16("indexes_of_2_dft_samples")),
            Per256(SBInt16("delta_offset_fft_filter_units")),
            Padding(18),
            Per2048(SBInt16("noise_power_measurement")),
            OneHundredth(SBInt16("agc_of_noise_power_measurement")),
            OneHundredthDecibel(UBInt16("reference_power_value")),
            Padding(10))

        self.mds_record = Struct(
            "mds_record",
            self.timestamp,
            SBInt8("quality_indicator"),
            Padding(3),
            UBInt32("source_packet_counter"),
            Padding(8),
            Array(20, self.waveform_data))

        self.mds = Array(self.n_records, self.mds_record)
        return self.mds


def envisat_get_mds_def(n_records, mds_target):
    """
    Returns the a construct parser for MDS records in Envisat SGDR files
    """
    valid_mds_targets = ["ra2", "wfm18hz"]
    if mds_target not in valid_mds_targets:
        raise ValueError("Unrecognized MDS target: "+mds_target)
    # Retrieve the corresponding parser class for the mds target
    mds_parser_class_name = "EnvisatSGDRMDS"+mds_target.upper()
    mds_parser_class = globals()[mds_parser_class_name]()
    mds_parser_class.n_records = n_records
    return mds_parser_class
