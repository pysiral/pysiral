# -*- coding: utf-8 -*-
"""
Created on Tue Jul 07 13:53:17 2015

@author: Stefan
"""

from pysiral.config import ConfigInfo
from pysiral.cryosat2.cryosat2_l1b import CryoSatL1B

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import glob
import os


def parse_cryosat_l1b():
    """ Sandbox script to parse a cryosat-l1b file """

    # Get Path information
    info = ConfigInfo()

    # Get an L1B SAR file
    l1b_directory = info.local_machine.local_l1b_repository.cryosat2.sar
    l1b_directory = os.path.join(l1b_directory, "2015", "04")
    l1b_files = glob.glob(os.path.join(l1b_directory, "*.DBL"))

    # Just select one
    l1b_file = l1b_files[0]

    # Parse the file into an L1B object
    l1b = CryoSatL1B()
    l1b.filename = l1b_file
    l1b.parse()

    # Get timing informatio
    timeorbit = l1b.mds[0].time_orbit[0]
    time_tai = get_tai_time(timeorbit.day, timeorbit.sec, timeorbit.msec)
    time_utc = tai2utc(time_tai)
    print time_tai

    # Get Waveform
    record = l1b.mds[0].waveform[0]
    measurement = l1b.mds[0].measurement[0]
    wfm = np.array(record.wfm).astype(np.float32)
    wfm_power = get_wfm_power(wfm, record.linear_scale, record.power_scale)

    # Calculate range for wavforme
    wfm_range = get_wfm_range(measurement.window_delay, len(wfm))

    # Plot the figure
    plt.figure()
    plt.plot(wfm_range, wfm_power)
    plt.show()


def get_wfm_power(counts, linear_scale, power_scale):
    """ Test function to return the waveform power in Watt """
    return counts*(linear_scale*1e-9)*2.0**(power_scale)


def get_wfm_range(window_delay, n_counts):
    """ Test function to return range for the waveform bins """
    lightspeed = 299792458.0
    bandwidth = 320000000.0
    # The two way delay time give the distance to the central bin
    central_bin_range = window_delay*lightspeed/2.0
    # Calculate the offset from the center to the first range bin
    window_size = (n_counts*lightspeed)/(4.0*bandwidth)
    first_bin_offset = window_size/2.0
    # Calculate the range increment for each bin
    range_increment = np.arange(n_counts)*lightspeed/(4.0*bandwidth)
    wfm_range = central_bin_range - first_bin_offset + range_increment
    return wfm_range


def get_tai_time(day, sec, msec):
    return datetime(2000, 1, 1)+timedelta(day, sec, msec)


def tai2utc(time_tai):
    fraction_of_day = float(time_tai.seconds)/86400.
    offset = pysofa.dat(
        time_tai.year,
        time_tai.month,
        time_tai.day,
        fraction_of_day)
    print offset

if __name__ == "__main__":
    parse_cryosat_l1b()
