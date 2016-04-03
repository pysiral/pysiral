# -*- coding: utf-8 -*-
"""
Created on Sun Apr 03 20:47:45 2016

@author: Stefan
"""

from pysiral.config import ConfigInfo
from pysiral.l1bdata import L1bConstructor
from pysiral.logging import stdout_logger
from pysiral.cryosat2.preproc import CryoSat2PreProcJob

import numpy as np

import glob
import os


def test_cryosat2_sin_bincount_reduction():

    # Init the logger
    log = stdout_logger("cs2-sin-bin-count-reduction")

    # Get Configuration Information
    config = ConfigInfo()

    # Get an L1B SAR file
    l1b_directory = config.local_machine.l1b_repository.cryosat2.sin
    l1b_directory = os.path.join(l1b_directory, "2015", "03")
    l1b_files = glob.glob(os.path.join(l1b_directory, "*.DBL"))

    # Read the file
    l1b = L1bConstructor(config)
    l1b.mission = "cryosat2"
    l1b.filename = l1b_files[1]
    l1b.construct()

    # Cut to the ocean section (nothing else matters)
    job = CryoSat2PreProcJob()
    job.log = log
    job.extract_polar_ocean_segments(l1b)

    # Show the initial waveform
    waveform_test_plot(l1b, pre=True)

    # Now reduce the bin count to SAR
    l1b.reduce_waveform_bin_count(256)

    # Show the reduced waveform
    waveform_test_plot(l1b)


def waveform_test_plot(l1b, pre=False):

    import matplotlib.pyplot as plt

    plt.figure("Waveform Test", figsize=(16, 8), facecolor="white")
    plt.imshow(np.log(l1b.waveform.power.T), cmap=plt.get_cmap("magma"),
               interpolation="none")
    if pre:
        lead_bins = int(0.4*256)
        trail_bins = 256-lead_bins
        max_index = np.argmax(l1b.waveform.power, axis=1)
        top = max_index - lead_bins
        bottom = max_index + trail_bins
        x = np.arange(len(max_index))
        plt.scatter(x, max_index, marker="o", facecolor="none", color="white")
        plt.scatter(x, top, marker="_", color="white")
        plt.scatter(x, bottom, marker="_", color="white")
    plt.show(block=True)

if __name__ == "__main__":
    test_cryosat2_sin_bincount_reduction()
