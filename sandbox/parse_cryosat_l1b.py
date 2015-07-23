# -*- coding: utf-8 -*-
"""
Created on Tue Jul 07 13:53:17 2015

@author: Stefan
"""

from pysiral.config import ConfigInfo
from pysiral.cryosat2.l1bfile import CryoSatL1B
from pysiral.cryosat2.functions import (
    get_tai_datetime_from_timestamp, get_cryosat2_wfm_power,
    get_cryosat2_wfm_range, tai2utc)
import numpy as np
import matplotlib.pyplot as plt

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
    l1b.post_processing()

    plt.figure()
    plt.plot(get_struct_field(l1b.corrections, "dry_troposphere"))
    plt.show()

    # Get timing information
    timeorbit = l1b.mds[0].time_orbit[0]
    time_tai = get_tai_datetime_from_timestamp(timeorbit.tai_timestamp)
    time_utc = tai2utc(time_tai)

    # Get Waveform
    record = l1b.mds[0].waveform[0]
    measurement = l1b.mds[0].measurement[0]
    wfm = np.array(record.wfm).astype(np.float32)
    wfm_power = get_cryosat2_wfm_power(
        wfm, record.linear_scale, record.power_scale)

    # Calculate range for wavforme
    wfm_range = get_cryosat2_wfm_range(measurement.window_delay, len(wfm))

    # Plot the figure
    plt.figure()
    plt.plot(wfm_range, wfm_power)
    plt.show()


def get_struct_field(struct, field):
    return np.array([record[field] for record in struct])

if __name__ == "__main__":
    parse_cryosat_l1b()
