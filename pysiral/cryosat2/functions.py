# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 14:48:01 2015

@author: Stefan
"""

import parse
import xmltodict
import numpy as np
import numpy.typing as npt
from loguru import logger
from attrdict import AttrDict
from pathlib import Path
from dateutil import parser as dtparser
from datetime import datetime, timedelta

from pysiral.l1bdata import Level1bData
from pysiral.l1preproc.procitems import L1PProcItem
from pysiral.waveform import get_waveforms_peak_power, get_footprint_sar, get_sigma0_sar


class L1PWaveformResampleSIN(L1PProcItem):
    """
    A post-processing class for CryoSat-2 l1p data that resamples the SIN waveform group
    to the same size as SAR waveform group
    """

    def __init__(self, **cfg) -> None:
        super(L1PWaveformResampleSIN, self).__init__(**cfg)

    def apply(self, l1: Level1bData) -> None:
        """
        API class for the Level-1 pre-processor. Functionality is reduce the size of the waveform power and
        range arrays for SIN data to the one for SAR data.
        :param l1: A Level-1 data instance
        :return: None, Level-1 object is change in place
        """

        # The functionality to reduce the waveform bins is built-in in the Level-1 data object
        if l1.radar_modes == "sin" and self.sin_target_bins is not None:
            l1.reduce_waveform_bin_count(self.sin_target_bins)


class L1PWaveformPadLRM(L1PProcItem):
    """
    A post-processing class for CryoSat-2 l1p data that resamples the SIN waveform group
    to the same size as SAR waveform group
    """

    def __init__(self, **cfg) -> None:
        super(L1PWaveformPadLRM, self).__init__(**cfg)

    def apply(self, l1: Level1bData) -> None:
        """
        API class for the Level-1 pre-processor. Functionality is to reduce the size
        of the waveform power and range arrays for SIN data to the one for SAR data.

        :param l1: A Level-1 data instance

        :return: None, Level-1 object is change in place
        """

        # The functionality to reduce the waveform bins is built-in in the Level-1 data object
        if l1.radar_modes == "lrm" and self.lrm_target_bins is not None:
            l1.increase_waveform_bin_count(self.lrm_target_bins)


class L1PCryoSat2Sigma0(L1PProcItem):
    """
    Class to compute radar backscatter coefficient (sigma0) directly
    from waveform data.
    """

    def __init__(self, **cfg):
        super(L1PCryoSat2Sigma0, self).__init__(**cfg)

    def apply(self, l1: Level1bData):
        """
        Compute the radar backscatter coefficient across all CryoSat-2 radar modes.
        This method first gets the variables needed for sigma0 computation of all
        radar modes and the leave the specifics to the respective sub_methods.

        Finally, two variables (sigma0 and peak_power) are added to the classifier
        container of the Level-1 data container.

        :param l1: The Level-1 data container

        :return: None, Level-1 data container is changed in place
        """

        # The output variable
        sigma0 = np.full(l1.n_records, np.nan)

        # Get the radar mode
        radar_mode = l1.get_parameter_by_name("waveform", "radar_mode")

        # Get the peak power as measure for the Rx Power
        rx_power = get_waveforms_peak_power(l1.waveform.power)
        tx_power = l1.get_parameter_by_name("classifier", "transmit_power")
        if tx_power is None:
            msg = "classifier `transmit_power` must exist for this pre-processor item -> aborting"
            logger.warning(msg)
            return

        # The computation of sigma0 requires the range to the surface as input
        # and is a step after the retracker. Here we use the satellite altitude
        # as an approximation for sea ice and ocean surfaces.
        altitude = l1.time_orbit.altitude

        # Sigma0 for LRM data
        is_lrm = radar_mode == 0
        if np.any(is_lrm):
            lrm_idxs = np.where(is_lrm)
            sigma0_lrm = self._get_sigma0_lrm(l1, lrm_idxs, rx_power, tx_power, altitude)
            sigma0[lrm_idxs] = sigma0_lrm[lrm_idxs]

        # Sigma0 for SAR/SARin data
        is_sar_sin = np.logical_not(is_lrm)
        if np.any(is_sar_sin):
            sar_sin_idxs = np.where(is_sar_sin)
            sigma0_sar_sin = self._get_sigma0_sar_sin(l1, sar_sin_idxs, rx_power, tx_power, altitude)
            sigma0[sar_sin_idxs] = sigma0_sar_sin[sar_sin_idxs]

        import matplotlib.pyplot as plt
        plt.figure(dpi=150)
        plt.plot(sigma0)
        plt.figure(dpi=150)
        plt.plot(radar_mode)
        plt.show()


        # Add the classifier
        l1.classifier.add(rx_power, "peak_power")
        l1.classifier.add(sigma0, "sigma0")

    def _get_sigma0_sar_sin(self,
                            l1: Level1bData,
                            sar_sin_idxs: npt.NDArray,
                            rx_power: npt.NDArray,
                            tx_power: npt.NDArray,
                            rng: npt.NDArray,
                            ) -> npt.NDArray:
        """
        Compute sigma0 for SAR/SARin

        :param l1: The Level-1 data container
        :param rx_power:
        :param tx_power:

        :return:
        """

        footprint_sar_kwargs = self.cfg.get("footprint_sar_kwargs", {})
        sigma0_sar_kwargs = self.cfg.get("sigma0_sar_kwargs", {})
        velocity = self._get_orbit_velocity_from_l1(l1)
        sigma0 = np.full(rx_power.shape, np.nan)

        # TODO: Can this be vectorized?
        for i in sar_sin_idxs:
            footprint_area = get_footprint_sar(rng[i], velocity[i], **footprint_sar_kwargs)
            sigma0[i] = get_sigma0_sar(rx_power[i], tx_power[i], rng[i], footprint_area,
                                       **sigma0_sar_kwargs)

        return sigma0

    @staticmethod
    def _get_sigma0_lrm(l1: Level1bData,
                        lrm_idxs: npt.NDArray,
                        rx_power: npt.NDArray,
                        tx_power: npt.NDArray,
                        rng: npt.NDArray,
                        ) -> npt.NDArray:
        """
        Compute sigma0 for LRM

        :param l1: The Level-1 data container
        :param rx_power:
        :param tx_power:

        :return:
        """
        # sigma0_lrm_kwargs = self.cfg.get("sigma0_lrm_kwargs", {})
        sigma0 = np.full(rx_power.shape, np.nan)
        roll_deg = l1.get_parameter_by_name("time_orbit", "antenna_roll")
        pitch_deg = l1.get_parameter_by_name("time_orbit", "antenna_pitch")

        # TODO: Can this be vectorized?
        for i in lrm_idxs:
            sigma0[i] = cryosat2_sigma0_lrm(rng[i], rx_power[i], roll_deg[i], pitch_deg[i], tx_power[i])

        return sigma0

    @staticmethod
    def _get_orbit_velocity_from_l1(l1: Level1bData) -> npt.NDArray:
        """
        Computes total orbital velocity

        :param l1:

        :return: velocity
        """
        sat_vel_x = l1.get_parameter_by_name("classifier", "satellite_velocity_x")
        sat_vel_y = l1.get_parameter_by_name("classifier", "satellite_velocity_y")
        sat_vel_z = l1.get_parameter_by_name("classifier", "satellite_velocity_z")
        if sat_vel_x is None or sat_vel_y is None or sat_vel_z is None:
            msg = "classifier `satellite_velocity_[x|y|z]` must exist for this pre-processor item -> aborting"
            raise ValueError(msg)
        return np.sqrt(sat_vel_x**2. + sat_vel_y**2. + sat_vel_z**2.)


# def get_sar_sigma0(wf_peak_power_watt, tx_pwr, r, v_s, **sigma0_par_dict):
#     """ Wrapper function to compute sigma nought for all waveforms """
#     n_records = wf_peak_power_watt.shape[0]
#     sigma0 = np.ndarray(shape=(n_records,))
#     for i in np.arange(n_records):
#         sigma0[i] = sar_sigma0(wf_peak_power_watt[i], tx_pwr[i], r[i], v_s[i], **sigma0_par_dict)
#     return sigma0
#
#
# def sar_sigma0(wf_peak_power_watt, tx_pwr, r, v_s,
#                wf_thermal_noise_watt=0.0, ptr_width=2.819e-09, tau_b=0.00352,
#                lambda_0=0.022084, wf=1., g_0=19054.607179632483,
#                bias_sigma0=0.0, l_atm=1.0, l_rx=1.0,
#                c_0=299792458.0, r_mean=6371000.0):
#     """
#     returns sigma nought (sigma0) for sar waveforms.
#     Applicable Documents
#     --------------------
#     Guidelines for reverting Waveform Power to Sigma Nought for
#     CryoSat-2 in SAR mode (v2.2), Salvatore Dinardo, 23/06/2016
#     XCRY-GSEG-EOPS-TN-14-0012
#     Arguments
#     ---------
#         wf_peak_power_watt (float)
#             waveform peak power in watt
#         tx_pwr (float)
#             transmitted peak power in watt
#         r (float)
#             range from satellite center of mass to surface reflection point
#             (to be appoximated by satellite altitude if no retracker range
#             available)
#         v_s (float)
#             satellite along track velocity in meter/sec
#     Keywords
#     --------
#         wf_thermal_noise_watt (float)
#             estimate of thermal noise power in watt (default: 0.0)
#             will be used to estimate waveform amplitude (Pu)
#         ptr_width (float)
#             3dB range point target response temporal width in seconds
#             (default: 2.819e-09 sec for CryoSat-2 SAR)
#         tau_b (float)
#             burst length in seconds
#             (default: 0.00352 sec for CryoSat-2 SAR)
#         lambda_0 (float)
#             radar wavelength in meter
#             (default: 0.022084 m for CryoSat-2 Ku Band altimeter)
#         wf (float)
#             footprint widening factor (1.486 * rv in case of Hamming window
#             application on burst data; rv: unspecified empirical factor)
#             (default: 1 no weighting window application)
#         g_0 (float)
#             antenna gain at boresight
#             (default: 10^(4.28) from document)
#         bias_sigma0 (float)
#             sigma nought bias
#             (default: 0.0)
#         l_atm (float)
#             two ways atmosphere losses (to be modelled)
#             (default: 1.0 (no loss))
#         l_rx (float)
#             receiving chain (RX) waveguide losses (to be characterized)
#             (default: 1.0 (no loss))
#         c_0 (float)
#             vacuum light speed in meter/sec
#         r_mean (float)
#             mean earth radius in meter
#     Returns
#     -------
#         sigma_0 (float)
#     """
#
#     # XXX: The definition of the variable Pu is not quite clear to me (Stefan)
#     # In the document it is referred to as "waveform power value in output
#     # of the re-tracking stage", however generally it is referred to as
#     # "waveform amplitude" that is obtained by a waveform function fit
#     # It is the scope of this function to provide a sigma0 estimate without
#     # proper retracking, therefore Pu is simply defined by the peak power
#     # and the thermal noise_power in watt
#     pu = wf_peak_power_watt + wf_thermal_noise_watt
#
#     # Intermediate steps & variables
#     pi = np.pi
#     alpha_earth = 1. + (r/r_mean)
#     lx = (lambda_0 * r)/(2. * v_s * tau_b)
#     ly = np.sqrt((c_0 * r * ptr_width)/alpha_earth)
#     a_sar = (2. * ly) * (wf * lx)
#     k = ((4.*pi)**3. * r**4. * l_atm * l_rx)/(lambda_0**2. * g_0**2. * a_sar)
#
#     # Final computation of sigma nought
#     sigma0 = 10. * np.log10(pu/tx_pwr) + 10. * np.log10(k) + bias_sigma0
#
#     return sigma0

def get_cryosat2_wfm_power(counts, linear_scale, power_scale):
    """
    Converts the Cryosat-2 waveform counts into a physical unit (Watts)

    Applicable Documents
    --------------------

    - CryoSat-2 User Handbook, L1B Parameters

    Arguments
    ---------
        counts (int list)
            unscaled echo waveform
            (typically: ``l1b.waveform.wfm``)
        linear_scale (int)
            linear scale factor
            (typically: ``l1b.waveform.linear_scale``)
        power_scale (int)
            power scale factor
            (typically: ``l1b.waveform.power_scale``)

    Returns
    -------
        float list with echo power in Watts
    """
    return counts*(linear_scale*1e-9)*2.0**power_scale


def get_cryosat2_wfm_range(window_delay, n_range_bins):
    """
    Calculates the range value of each range bin in the waveform

    Applicable Documents
    --------------------
    - CryoSat-2 User Handbook, Range Window and Window Delay

    Arguments
    ---------
    window_delay (float)
        two-way delay time in seconds from the CryoSat-2 data L1b file

    n_range_bins (int)
        Number of range bins (the length of the waveform list)

    Returns
    -------
        wfm_range
            range value in meterfrom the satellite to each range bin in the
            waveform

    Notes
    -----

    - It is assumed that the two way delay time (window delay) already
      contains instrument related range corrections
    - No other range corrections are applied here
    - The bandwidth of CryoSat-2 (320Mhz) is hard-coded here

    """
    lightspeed = 299792458.0
    bandwidth = 320000000.0
    # The two way delay time give the distance to the central bin
    central_window_range = window_delay*lightspeed/2.0
    # Calculate the offset from the center to the first range bin
    window_size = (n_range_bins*lightspeed)/(4.0*bandwidth)
    first_bin_offset = window_size/2.0
    # Calculate the range increment for each bin
    range_increment = np.arange(n_range_bins)*lightspeed/(4.0*bandwidth)
    wfm_range = central_window_range - first_bin_offset + range_increment
    return wfm_range


def get_cryosat2_wfm_range_userhandbook(window_delay, n_range_bins):
    lightspeed = 299792458.0
    bandwidth = 320000000.0
    range_bin_index = np.arange(n_range_bins)
    wfm_range = lightspeed/(4. * bandwidth) * \
        (2*window_delay*bandwidth - n_range_bins/2.0 + range_bin_index)
    return wfm_range


def get_tai_datetime_from_timestamp(mdsr_timestamp):
    """
    Converts the TAI MDSR timestamp into a datetime object

    Background
    ----------
    The timestamp in the CryoSat-2 data files is given as
    days since January 1, 2000, seconds of the day and microseconds of the
    day.

    Arguments
    ---------

    mdsr_timestamp (object or object list)
        any class object with attributes day, sec, msec
        (attributes must be of type int)

    Return
    ------

    datetime object with date and time in TAI

    """
    timestamps = np.asarray(mdsr_timestamp)
    output = np.ndarray(shape=timestamps.shape, dtype=object)
    for i, timestamp in enumerate(timestamps):
        output[i] = datetime(2000, 1, 1) + timedelta(
            timestamp.day, timestamp.sec, timestamp.msec)
    if len(output) > 0:
        return output
    else:
        return output[0]


def parse_cryosat_l1b_filename(filename):
    """
    Returns the information in the CryoSat-2 l1b filename
    """
    # Strip path and file extension
    filename = Path(filename).stem
    # Construct the parser
    parser_str = "CS_{proc_stage}_"
    parser_str += "{instrument}_"
    parser_str += "{radar_mode}_"
    parser_str += "{data_level}_"
    parser_str += "{start_dt}_"
    parser_str += "{stop_dt}_"
    parser_str += "{baseline}"
    parser = parse.compile(parser_str)
    # Parse the filename
    result = parser.parse(filename)
    # Do some post-editing
    # - parse is not that smart when it comes to work with date strings
    # - naming conventions as the rest of pysiral
    info = {"mission": "cryosat2", "instrument": result["instrument"].lower(),
            "radar_mode": result["radar_mode"].lower(), "data_level": "L" + result["data_level"],
            "start_dt": dtparser.parse(result["start_dt"]), "stop_dt": dtparser.parse(result["stop_dt"]),
            "baseline": result["baseline"]}
    return AttrDict(info)


def parse_cryosat_l1b_xml_header(filename):
    """
    Reads the XML header file of a CryoSat-2 L1b Data set
    and returns the contents as an OrderedDict
    """
    with open(str(filename)) as fd:
        content_odereddict = xmltodict.parse(fd.read())
    return content_odereddict[u'Earth_Explorer_Header']


# def get_sar_sigma0(wf_peak_power_watt, tx_pwr, r, v_s, **sigma0_par_dict):
#     """ Wrapper function to compute sigma nought for all waveforms """
#     n_records = wf_peak_power_watt.shape[0]
#     sigma0 = np.ndarray(shape=(n_records))
#     for i in np.arange(n_records):
#         sigma0[i] = sar_sigma0(wf_peak_power_watt[i], tx_pwr[i], r[i], v_s[i], **sigma0_par_dict)
#     return sigma0


def get_footprint_lrm(r: float, band_width: float = 320000000.0) -> float:
    """
    Compute the CryoSat-2 LRM footprint for variable range to the surface.

    Applicable Documents:

        Michele Scagliola, CryoSat Footprints ESA/Aresys, v1.2
        ESA document ref: XCRY-GSEG-EOPG-TN-13-0013

    :param r: range from satellite center of mass to surface reflection point
    (to be appoximated by satellite altitude if no retracker range available)
    :param band_width: CryoSat-2 pulse bandwidth in Hz

    :return sigma_0: Radar backscatter coefficient
    """

    c_0 = 299792458.0
    footprint_radius = np.sqrt(r * c_0 / band_width) / 1000.
    area_lrm = np.pi * footprint_radius ** 2.

    return area_lrm

def get_footprint_sar(r: float,
                      v_s: float,
                      ptr_width: float = 2.819e-09,
                      tau_b: float = 0.00352,
                      lambda_0: float = 0.022084,
                      wf: float = 1.0,
                      r_mean: float = 6371000.0
                      ) -> float:
    """
    Compute the CryoSat-2 SAR footprint for variable range to the surface.

    Applicable Documents:

        Michele Scagliola, CryoSat Footprints ESA/Aresys, v1.2
        ESA document ref: XCRY-GSEG-EOPG-TN-13-0013

    :param r: range from satellite center of mass to surface reflection point
    (to be appoximated by satellite altitude if no retracker range available)
    :param r: range from satellite center of mass to surface reflection point
        (to be appoximated by satellite altitude if no retracker range available)
    :param v_s: satellite along track velocity in meter/sec
    :param ptr_width: 3dB range point target response temporal width in seconds
        (default: 2.819e-09 sec for CryoSat-2 SAR)
    :param tau_b: burst length in seconds (default: 0.00352 sec for CryoSat-2 SAR)
    :param lambda_0: radar wavelength in meter (default: 0.022084 m for CryoSat-2 Ku Band altimeter)
    :param wf: footprint widening factor (1.486 * rv in case of Hamming window application on burst data;
        rv: unspecified empirical factor) (default: 1 no weighting window application)
    :param r_mean: mean earth radius in meter

    :return area_sar: The SAR footprint in square meters
    """
    c_0 = 299792458.0
    alpha_earth = 1. + (r / r_mean)
    lx = (lambda_0 * r) / (2. * v_s * tau_b)
    ly = np.sqrt((c_0 * r * ptr_width) / alpha_earth)
    return (2. * ly) * (wf * lx)


def cryosat2_sigma0_lrm(
        rng: float,
        p_rx: float,
        roll_deg: float,
        pitch_deg: float,
        p_tx: float = 22.4,
        sigma_bias: float = 0.0
        ) -> float:
    """
    Compute LRM backscatter with optional mispointing and bias correction.
    The code was supplied by David Brockley <d.brockley@ucl.ac.uk>.

    :param rng: range to the surface
    :param p_rx: rx power
    :param roll_deg:
    :param pitch_deg:
    :param p_tx: tx power
    :param sigma_bias: sigma0 bias (default: 0.0)

    :raises None:

    :return:
    """

    ellipse_semi_major_m = 6378137.0
    speed_of_light_ms = 2.99792458E+08
    effective_pulse_len_s = 4.183e-9
    ant_gain_linear = 18197.008586099826
    wavelength_m = 2.2084e-2
    r_nom = 720000.0
    misp_roll_bias_deg = 0.
    misp_pitch_bias_deg = 0.
    beam_angle_az_deg = 1.06
    beam_angle_el_deg = 1.1992

    gamma = (2.0 / np.log(2.0)) * np.power(np.sin(1.0 /
                                                  (1.0 / np.radians(beam_angle_az_deg) +
                                                   1.0 / np.radians(beam_angle_el_deg))), 2.0)
    power_ratio_db = 10.0 * np.log10(p_rx / p_tx)
    top = (speed_of_light_ms * np.pi * np.power(wavelength_m, 2.0) *
           np.power(ant_gain_linear, 2.0) * effective_pulse_len_s)
    earth_cor = (r_nom+ellipse_semi_major_m)/ellipse_semi_major_m
    c_const = 10.0 * np.log10(top/(np.power(4.0*np.pi*r_nom,3.0)*earth_cor))
    mispoint_rad = np.sqrt(
        np.power(np.radians(roll_deg) + np.radians(misp_roll_bias_deg), 2.0) +
        np.power(np.radians(pitch_deg) + np.radians(misp_pitch_bias_deg), 2.0)
    )
    mispoint_cor_lin = -4.0 * np.power(np.sin(mispoint_rad), 2.0) / gamma
    misp_db = 10.0 * mispoint_cor_lin / np.log(10.0)
    return power_ratio_db - c_const + 30.0 * np.log10(rng / r_nom) - misp_db + sigma_bias
