# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 14:48:01 2015

@author: Stefan
"""

import parse
import xmltodict
import numpy as np
from attrdict import AttrDict
from pathlib import Path
from dateutil import parser as dtparser
from datetime import datetime, timedelta

from pysiral.logging import DefaultLoggingClass

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
            self.log.warning(msg)

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
    return counts*(linear_scale*1e-9)*2.0**(power_scale)


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
    output = np.ndarray(shape=(len(timestamps)), dtype=object)
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
    info = {}
    info["mission"] = "cryosat2"
    info["instrument"] = result["instrument"].lower()
    info["radar_mode"] = result["radar_mode"].lower()
    info["data_level"] = "L"+result["data_level"]
    info["start_dt"] = dtparser.parse(result["start_dt"])
    info["stop_dt"] = dtparser.parse(result["stop_dt"])
    info["baseline"] = result["baseline"]
    return AttrDict(info)


def parse_cryosat_l1b_xml_header(filename):
    """
    Reads the XML header file of a CryoSat-2 L1b Data set
    and returns the contents as an OrderedDict
    """
    with open(filename) as fd:
        content_odereddict = xmltodict.parse(fd.read())
    return content_odereddict[u'Earth_Explorer_Header']

