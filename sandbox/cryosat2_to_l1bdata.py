# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 15:18:59 2015

@author: Stefan
"""

import os
import glob
import argparse

from pysiral.config import ConfigInfo
from pysiral.l1bdata import L1bConstructor

import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
import numpy as np
import time


def cryosat2_to_l1bdata():
    """
    Sandbox script to create an L1BData object from a CryoSat-2 L1b data file
    """

    """ parse command line arguments """
    parser = get_cs2l1bdata_argparser()
    args = parser.parse_args()

    # Get Configuration Information
    config = ConfigInfo()

    # Get an L1B SAR file
    l1b_directory = config.local_machine.l1b_repository.cryosat2[args.mode]
    l1b_directory = os.path.join(l1b_directory, args.month[0], args.month[1])
    l1b_files = glob.glob(os.path.join(l1b_directory, "*.DBL"))

    # Read the file
    t0 = time.time()
    l1b = L1bConstructor(config)
    l1b.mission = "cryosat2"
    l1b.filename = l1b_files[args.file_number]
    l1b.get_header_info()
    print l1b.filename
    print l1b.info.open_ocean_percent
    print l1b.info.lat_min
    print l1b.info.lat_max
    print l1b.info.lon_min
    print l1b.info.lon_max
    l1b.construct()

    t1 = time.time()
    print "Constructing CryoSat-2 l1bdata object in %.3g seconds" % (t1 - t0)

    # Quick plots
    cryosat2_l1b_orbit_plot(l1b)
    cryosat2_l1b_corrections_plot(l1b)
    cryosat2_l1b_waveform_plot(l1b)


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
    m.scatter(x[0], y[0], color="#003e6e")

    plt.title("")
    plt.show(block=False)


def cryosat2_l1b_corrections_plot(l1b):

    from matplotlib.ticker import MultipleLocator

    n = len(l1b.correction.parameter_list)
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
    plt.show(block=False)


def cryosat2_l1b_waveform_plot(l1b):

    plt.figure(facecolor="white", figsize=(12, 6))
    ax = plt.gca()

    image, elevation_range = align_echo_power(
        l1b.waveform.power, l1b.waveform.range, l1b.time_orbit.altitude)
    image_extent = (0, len(l1b.time_orbit.altitude),
                    np.amin(elevation_range), np.amax(elevation_range))

    ref_aspect = 0.4
    num_recs, num_range_bins = image.shape
    waveform_aspect = float(num_recs)/float(num_range_bins)
    aspect = waveform_aspect / ref_aspect

    im = ax.imshow(
        np.log10(image).transpose(), cmap=plt.get_cmap("magma"),
        interpolation='none', origin='lower', extent=image_extent,
        aspect=aspect)

    ax.yaxis.set_tick_params(direction='out')
    ax.yaxis.set_ticks_position('left')
    ax.set_ylabel("Elevation (m)")
    ax.xaxis.set_ticks([])
    ax.spines["left"].set_position(("data", -20))

    spines_to_remove = ["top", "right", "bottom"]
    for spine in spines_to_remove:
        ax.spines[spine].set_visible(False)

    cb = plt.colorbar(im, orientation="vertical", shrink=0.5)
    cb.set_label("Echo Power log(dB)")
    plt.tight_layout()
    plt.show()


def align_echo_power(power, range, altitude):

    from scipy.interpolate import interp1d

    n_range_bins = len(power[0, :])
    n_records = len(power[:, 0])

    elevation = np.repeat(altitude, n_range_bins)
    elevation = elevation.reshape(n_records, n_range_bins)
    elevation -= range

    range_step = range[0, 1] - range[0, 0]

    elevation_range = np.arange(np.amin(elevation)-range_step*0.5,
                                np.amax(elevation)+range_step*0.5,
                                range_step)

    aligned_power = np.ndarray(shape=(n_records, len(elevation_range)))
    aligned_power *= np.nan

    for i in np.arange(n_records):
        f = interp1d(elevation[i, :].flatten(), power[i, :].flatten(),
                     bounds_error=False, fill_value=np.nan)
        aligned_power[i, :] = f(elevation_range)

    return aligned_power, elevation_range


def get_cs2l1bdata_argparser():
    """ Handle command line arguments """
    parser = argparse.ArgumentParser()
    # Mission id string: cryosat2, envisat, ...
    # Start month as list: [yyyy, mm]
    parser.add_argument(
        '-month', action='store', dest='month',
        nargs='+', type=str, required=True,
        help='start date as year and month (-t0 yyyy mm)')

    parser.add_argument(
        '-radarmode', action='store', dest='mode',
        type=str, required=True,
        help='radar mode [sar|sin]')

    # Add an Option to skip a number of files (e.g. for a restart)
    parser.add_argument(
        "-f", "-file", action='store', type=int, const=0, nargs='?',
        dest='file_number', help='file number in month')

    return parser

if __name__ == "__main__":
    mpl.rcParams['font.sans-serif'] = "arial"

    for target in ["xtick.color", "ytick.color", "axes.edgecolor",
                   "axes.labelcolor"]:
        mpl.rcParams[target] = "#4b4b4d"

    cryosat2_to_l1bdata()
