# -*- coding=utf-8 -*-
"""
Created on Fri Jul 10 21:00:37 2015

@author=Stefan
"""

from collection import namedtuple


class Cryosat2MSDDef_C(object):
    """ Holds the Definition of L1b MDS records """

    def __init__(self):

        # Default definition of fields
        mds_def = namedtuple('mds_def', ['name', 'fmt', 'size'])

        self.time_orbit_group = [
            mds_def(name="day", fmt="l", size=4),
            mds_def(name="sec", fmt="l", size=4),
            mds_def(name="msec", fmt="l", size=4),
            mds_def(name="uso_corr", fmt="l", size=4),
            mds_def(name="mode_id", fmt="H", size=2),
            mds_def(name="source_sequence_counter", fmt="H", size=2),
            mds_def(name="instrument_configuration", fmt="H", size=2),
            mds_def(name="burst_counter", fmt="l", size=4),
            mds_def(name="measurement_latitude", fmt="l", size=4),
            mds_def(name="measurement_longtitude", fmt="l", size=4),
            mds_def(name="altitude_cog", fmt="l", size=4),
            mds_def(name="altitude_rate", fmt="l", size=4),
            mds_def(name="satellite_velocity", fmt="3l", size=12),
            mds_def(name="real_beam", fmt="3l", size=12),
            mds_def(name="interferometry_baseline", fmt="3l", size=12),
            mds_def(name="star_tracker_usage", fmt="H", size=2),
            mds_def(name="antenna_roll", fmt="l", size=4),
            mds_def(name="antenna_pitch", fmt="l", size=4),
            mds_def(name="antenna_yaw", fmt="l", size=4),
            mds_def(name="fbr_confidence_flag", fmt="l", size=4)]

        # measurement group
        self.measurement_group = [
            - {name: td, struct_fmt: "q"}
            - {name: h0, struct_fmt: "l"}
            - {name: cor2, struct_fmt: "l"}
            - {name: lai, struct_fmt: "l"}
            - {name: fai, struct_fmt: "l"}
            - {name: agc_ch1, struct_fmt: "l"}
            - {name: agc_ch2, struct_fmt: "l"}
            - {name: tr_gain_ch1, struct_fmt: "l"}
            - {name: tr_gain_ch2, struct_fmt: "l"}
            - {name: tx_power, struct_fmt: "l"}
            - {name: dopp_rc, struct_fmt: "l"}
            - {name: tr_inst_rc, struct_fmt: "l"}
            - {name: r_inst_rc, struct_fmt: "l"}
            - {name: tr_inst_gain_c, struct_fmt: "l"}
            - {name: r_inst_gain_c, struct_fmt: "l"}
            - {name: int_phase_c, struct_fmt: "l"}
            - {name: ext_phase_c, struct_fmt: "l"}
            - {name: noise_pwr, struct_fmt: "l"}
            - {name: phase_slope_c, struct_fmt: "l"}
            - {name: spares, struct_fmt: "4b"}

        # Geocorrections group
        geocorrections_group:
            - {name: dry_c, struct_fmt: "l"}
            - {name: wet_c, struct_fmt: "l"}
            - {name: ib_c, struct_fmt: "l"}
            - {name: dac_c, struct_fmt: "l"}
            - {name: iono_gim, struct_fmt: "l"}
            - {name: iono_mod, struct_fmt: "l"}
            - {name: h_ot, struct_fmt: "l"}
            - {name: h_lpeot, struct_fmt: "l"}
            - {name: h_olt, struct_fmt: "l"}
            - {name: h_set, struct_fmt: "l"}
            - {name: h_gpt, struct_fmt: "l"}
            - {name: surf_type, struct_fmt: "L"}
            - {name: spare1, struct_fmt: "4b"}
            - {name: corr_status, struct_fmt: "L"}
            - {name: corr_error, struct_fmt: "L"}
            - {name: spare2, struct_fmt: "4b"}

        # 1Hz Average SAR Waveform group
        onehz_wavefrom_group:
            - {name: day_1hz, struct_fmt: "l"}
            - {name: sec_1hz, struct_fmt: "L"}
            - {name: micsec_1hz, struct_fmt: "L"}
            - {name: lat_1hz, struct_fmt: "l"}
            - {name: lon_1hz, struct_fmt: "l"}
            - {name: alt_1hz, struct_fmt: "l"}
            - {name: td_1hz, struct_fmt: "q"}
            - {name: avg_wfm, struct_fmt: "128I"}
            - {name: linear_wfm_mult, struct_fmt: "l"}
            - {name: power2_wfm_mult, struct_fmt: "l"}
            - {name: num_avg_echoes, struct_fmt: "H"}
            - {name: flags, struct_fmt: "H"}

        # 20 Hz SAR Waveform group
        sar_waveform_group:
            - {name: wfm, struct_fmt: "256I"}
            - {name: linear_wfm_multiplier, struct_fmt: "l"}
            - {name: power2_wfm_multiplier, struct_fmt: "l"}
            - {name: num_avg_echoes, struct_fmt: "H"}
            - {name: flags, struct_fmt: "H"}
            - {name: beam_sd, struct_fmt: "H"}
            - {name: beam_centre, struct_fmt: "H"}
            - {name: beam_amplitude, struct_fmt: "H"}
            - {name: beam_skew, struct_fmt: "h"}
            - {name: beam_kurt, struct_fmt: "h"}
            - {name: beam_spare, struct_fmt: "50b"}