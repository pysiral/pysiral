# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 14:48:01 2015

@author: Stefan
"""

import parse
import xmltodict
import numpy as np
from loguru import logger
from attrdict import AttrDict
from pathlib import Path
from dateutil import parser as dtparser
from datetime import datetime, timedelta

from pysiral.core import DefaultLoggingClass


class L1PWaveformResampleSIN(DefaultLoggingClass):
    """
    A post-processing class for CryoSat-2 l1p data that resamples the SIN waveform group
    to the same size as SAR waveform group
    """

    def __init__(self, **cfg):
        cls_name = self.__class__.__name__
        super(L1PWaveformResampleSIN, self).__init__(cls_name)
        self.sin_target_bins = cfg.get("sin_target_bins", None)
        if self.sin_target_bins is None:
            msg = "Missing option `sin_target_bins` -> SIN waveform will not be resampled!"
            logger.warning(msg)

    def apply(self, l1):
        """
        API class for the Level-1 pre-processor. Functionality is reduce the size of the waveform power and
        range arrays for SIN data to the one for SAR data.
        :param l1: A Level-1 data instance
        :return: None, Level-1 object is change in place
        """

        # The functionality to reduce the waveform bins is built-in in the Level-1 data object
        if l1.radar_modes == "sin" and self.sin_target_bins is not None:
            l1.reduce_waveform_bin_count(self.sin_target_bins)


class L1PWaveformPadLRM(DefaultLoggingClass):
    """
    A post-processing class for CryoSat-2 l1p data that resamples the SIN waveform group
    to the same size as SAR waveform group
    """

    def __init__(self, **cfg):
        cls_name = self.__class__.__name__
        super(L1PWaveformPadLRM, self).__init__(cls_name)
        self.lrm_target_bins = cfg.get("lrm_target_bins", None)
        if self.lrm_target_bins is None:
            msg = "Missing option `lrm_target_bins` -> LRM waveform will not be padded!"
            logger.warning(msg)

    def apply(self, l1):
        """
        API class for the Level-1 pre-processor. Functionality is reduce the size of the waveform power and
        range arrays for SIN data to the one for SAR data.
        :param l1: A Level-1 data instance
        :return: None, Level-1 object is change in place
        """

        # The functionality to reduce the waveform bins is built-in in the Level-1 data object
        if l1.radar_modes == "lrm" and self.lrm_target_bins is not None:
            l1.increase_waveform_bin_count(self.lrm_target_bins)


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
    area_sar = (2. * ly) * (wf * lx)

    return area_sar
