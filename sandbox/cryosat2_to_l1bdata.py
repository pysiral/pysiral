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
from mpl_toolkits.basemap import Basemap
import numpy as np


def cryosat2_to_l1bdata():
    """
    Sandbox script to create an L1BData object from a CryoSat-2 L1b data file
    """

    # Get Configuration Information
    config = ConfigInfo()

    # Get an L1B SAR file
    l1b_directory = config.local_machine.local_l1b_repository.cryosat2.sar
    l1b_directory = os.path.join(l1b_directory, "2015", "04")
    l1b_files = glob.glob(os.path.join(l1b_directory, "*.DBL"))

    # Read the file
    l1b = L1bConstructor(config)
    l1b.mission = "cryosat2"
    l1b.filename = l1b_files[0]
    l1b.construct()

    # Quick plots
    cryosat2_l1b_orbit_plot(l1b)
    cryosat2_l1b_corrections_plot(l1b)


def cryosat2_l1b_orbit_plot(l1b):

    # AWI eisblau #00ace5
    # AWI tiefblau #003e6e
    # AWI grau 1 #4b4b4d
    # AWI grau 2 #bcbdbf

    grid_keyw = {"dashes": (None, None),
                 "color": "#bcbdbf",
                 "linewidth": 0.5,
                 "latmax": 88}

    lat_0 = 80.0
    if np.nanmean(l1b.time_orbit.latitude) < 0:
        lat_0 *= -1.0

    plt.figure(facecolor="white")
    m = Basemap(projection='ortho', lon_0=0, lat_0=lat_0, resolution='i')
    m.fillcontinents(color='#00ace5', lake_color='#00ace5')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-90., 120., 10.), **grid_keyw)
    m.drawmeridians(np.arange(0., 420., 15.), **grid_keyw)

    x, y = m(l1b.time_orbit.longitude, l1b.time_orbit.latitude)
    m.plot(x, y, color="#003e6e", linewidth=2.0)

    plt.title("")
    plt.show(block=False)


def cryosat2_l1b_corrections_plot(l1b):

    from matplotlib.ticker import MultipleLocator

    n = len(l1b.correction.list)
    f, ax = plt.subplots(n, sharex=True, facecolor="white", figsize=(10, 16))
    for i in np.arange(n):
        correction, name = l1b.correction.get_parameter_by_index(i)
        ax[i].plot(correction, lw=2, color="#00ace5")
        ax[i].set_title(name)
        ax[i].yaxis.set_minor_locator(MultipleLocator(0.02))
        ax[i].yaxis.grid(True, which='minor')
        ax[i].yaxis.set_tick_params(direction='out')
        ax[i].yaxis.set_ticks_position('left')
        ax[i].xaxis.set_ticks([])

        spines_to_remove = ["top", "right", "bottom"]
        for spine in spines_to_remove:
            ax[i].spines[spine].set_visible(False)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    cryosat2_to_l1bdata()
