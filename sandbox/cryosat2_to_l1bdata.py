# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 15:18:59 2015

@author: Stefan
"""

import os
import glob

from pysiral.config import ConfigInfo
from pysiral.l1bdata import L1bConstructor

import matplotlib.pyplot as plt


def cryosat2_to_l1bdata():
    """
    Sandbox script to create an L1BData object from a CryoSat-2 L1b data file
    """

    # Get Path information
    info = ConfigInfo()

    # Get an L1B SAR file
    l1b_directory = info.local_machine.local_l1b_repository.cryosat2.sar
    l1b_directory = os.path.join(l1b_directory, "2015", "04")
    l1b_files = glob.glob(os.path.join(l1b_directory, "*.DBL"))

    l1b = L1bConstructor()
    l1b.mission = "cryosat2"
    l1b.filename = l1b_files[0]
    l1b.construct()

    plt.figure()
    plt.plot(l1b.time_orbit.longitude, l1b.time_orbit.latitude)
    plt.show()

if __name__ == "__main__":
    cryosat2_to_l1bdata()
