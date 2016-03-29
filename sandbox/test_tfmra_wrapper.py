# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 11:32:49 2015

@author: shendric
"""

import os
import glob
import cPickle as pickle

import numpy as np

from pysiral.config import ConfigInfo
from pysiral.l1bdata import L1bConstructor

from tfmra_wrapper import pytfmra

import matplotlib.pyplot as plt


def test_tfmra_wrapper():
    # # Get Configuration Information
    config = ConfigInfo()

    # # Get an L1B SAR file
    l1b_directory = config.local_machine.l1b_repository.cryosat2.sar
    l1b_directory = os.path.join(l1b_directory, "2015", "04")
    l1b_files = glob.glob(os.path.join(l1b_directory, "*.DBL"))

    # Read the file
    l1b = L1bConstructor(config)
    l1b.mission = "cryosat2"
    l1b.filename = l1b_files[0]
    l1b.construct()

    pickle_file = r"test_tfmra_l1b_object.pickle"
    with open(pickle_file, "wb") as fh:
        pickle.dump(
            [l1b.n_records, l1b.waveform.range, l1b.waveform.power], fh)

    with open(pickle_file, "r") as fh:
        l1b_n_records, l1b_waveform_range, l1b_waveform_power = \
            pickle.load(fh)

    retracker_range = np.ndarray(shape=(l1b_n_records), dtype=np.float32)
    # This is the relevant call for the retracker wrapper funtion
    for index in np.arange(l1b_n_records):
        retracker_range[index] = pytfmra(l1b_waveform_range[index, :],
                                         l1b_waveform_power[index, :])

    elevation = retracker_range
    plt.figure()
    plt.plot(elevation)
    plt.show()

if __name__ == "__main__":
    test_tfmra_wrapper()
