# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 11:32:49 2015

@author: shendric
"""

import os
import glob

import numpy as np

from pysiral.config import ConfigInfo
from pysiral.l1bdata import L1bConstructor

from tfmra_wrapper import pytfmra

import matplotlib.pyplot as plt


def test_tfmra_wrapper():
    # Get Configuration Information
    config = ConfigInfo()

    # Get an L1B SAR file
    l1b_directory = config.local_machine.l1b_repository.cryosat2.sar
    l1b_directory = os.path.join(l1b_directory, "2015", "04")
    l1b_files = glob.glob(os.path.join(l1b_directory, "*.DBL"))

    # Read the file
    l1b = L1bConstructor(config)
    l1b.mission = "cryosat2"
    l1b.filename = l1b_files[0]
    l1b.construct()

    retracker_range = np.ndarray(shape=(l1b.n_records), dtype=np.float32)
    # This is the relevant call for the retracker wrapper funtion
    for index in np.arange(l1b.n_records):
        retracker_range[index] = pytfmra(l1b.waveform.range[index, :],
                                         l1b.waveform.power[index, :])

    elevation = l1b.time_orbit.altitude - retracker_range
    plt.figure()
    plt.plot(elevation)
    plt.show()

if __name__ == "__main__":
    test_tfmra_wrapper()
