# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 14:48:01 2015

@author: Stefan
"""

import os
import parse
import xmltodict
import numpy as np
from treedict import TreeDict
from dateutil import parser as dtparser
from datetime import datetime, timedelta


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
    filename = os.path.basename(filename)
    filename = os.path.splitext(filename)[0]
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
    info = TreeDict()
    info.mission = "cryosat2"
    info.instrument = result["instrument"].lower()
    info.radar_mode = result["radar_mode"].lower()
    info.data_level = "L"+result["data_level"]
    info.start_dt = dtparser.parse(result["start_dt"])
    info.stop_dt = dtparser.parse(result["stop_dt"])
    info.baseline = result["baseline"]
    return info


def parse_cryosat_l1b_xml_header(filename):
    """
    Reads the XML header file of a CryoSat-2 L1b Data set
    and returns the contents as an OrderedDict
    """
    with open(filename) as fd:
        content_odereddict = xmltodict.parse(fd.read())
    return content_odereddict[u'Earth_Explorer_Header']


def tai2utc(tai_datetime):
    """
    Converts a datetime object in TAI into UTC
    TODO: No real functionality here yet, maybe parse leap second list from
          IETF?: http://www.ietf.org/timezones/data/leap-seconds.list
    """
    utc_offset_seconds = 33
    tai_times = np.asarray(tai_datetime)
    output = np.ndarray(shape=(len(tai_times)), dtype=object)
    for i, tai_time in enumerate(tai_times):
        output[i] = tai_time+timedelta(seconds=utc_offset_seconds)
    if len(output) > 0:
        return output
    else:
        return output[0]
