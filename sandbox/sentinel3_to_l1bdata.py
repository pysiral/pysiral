# -*- coding: utf-8 -*-
"""
Created on Tue Mar 08 18:42:01 2016

@author: Stefan
"""

from pysiral.sentinel3.iotools import get_sentinel3_l1b_filelist
from pysiral.config import ConfigInfo
from pysiral.l1bdata import L1bConstructor

import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap

import os
import argparse
import numpy as np
import time


def sentinel3_to_l1bdata():
    """
    Sandbox script to create an L1BData object from a ERS-2 L1b data file
    """

    """ parse command line arguments """
    parser = get_s3l1bdata_argparser()
    args = parser.parse_args()

    # Get Configuration Information
    config = ConfigInfo()

    # Get an L1B SAR file
    l1b_directory = config.local_machine.l1b_repository.sentinel3a.sral
    l1b_directory = os.path.join(l1b_directory, args.month[0], args.month[1])
    l1b_files = get_sentinel3_l1b_filelist(l1b_directory)

    # Read the file
    t0 = time.time()
    l1b = L1bConstructor(config)
    l1b.mission = "sentinel3a"
    l1b.filename = l1b_files[args.file_number]
    l1b.construct()
    t1 = time.time()
    print "Constructing Sentinel-3A l1bdata object in %.3g seconds" % (t1 - t0)

    # Quick plots
    # s3_l1b_fullorbit_plot(l1b)

    # Limit amounts of waveforms to plot => memory error
    # l1b.trim_to_subset(np.arange(1000)+9000)

    # Quick plots
    s3_l1b_orbit_plot(l1b)
    s3_l1b_corrections_plot(l1b)
    s3_l1b_classifier_plot(l1b)
    s3_l1b_waveform_plot(l1b)


def s3_l1b_orbit_plot(l1b):

    # AWI eisblau #00ace5
    # AWI tiefblau #003e6e
    # AWI grau 1 #4b4b4d
    # AWI grau 2 #bcbdbf

    projection = {
        "north": {
            "projection": "ortho",
            "lon_0": -45,
            "lat_0": 80,
            "resolution": "i"},
        "south": {
            "projection": "ortho",
            "lon_0": 0,
            "lat_0": -70,
            "resolution": "i"}
            }

    hemisphere = "north"
    if l1b.time_orbit.latitude[0] < 0:
        hemisphere = "south"

    xsize, ysize = 4000, 2000

    lats = np.linspace(-0.5*np.pi, 0.5*np.pi, num=ysize)

    if hemisphere == "south":
        shading = np.zeros(shape=(ysize, xsize))
        indices = np.where(lats > np.deg2rad(-90))[0]
        colormap_name = "Greys"
    elif hemisphere == "north":
        shading = np.ones(shape=(ysize, xsize)) * 2.0
        indices = np.where(lats < np.deg2rad(90))[0]
        colormap_name = "Greys_r"
    for i in indices:
        shading[i, :] = (np.sin(lats[i])+1.0)**3.

#    plt.figure()
#    plt.imshow(np.flipud(shading), cmap=plt.get_cmap(colormap_name))
#    plt.colorbar()
#    plt.show()
#    stop

    from pysiral.iotools import get_temp_png_filename
    temp_filename = get_temp_png_filename()
    fig = plt.figure(figsize=(10, 5))
    ax = plt.axes([0, 0, 1, 1])
    plt.axis('off')
    ax.imshow(np.flipud(shading), cmap=plt.get_cmap(colormap_name),
              vmin=0, vmax=1.1*np.amax(shading))
    plt.savefig(temp_filename, dpi=600)
    plt.close(fig)

    grid_keyw = {"dashes": (None, None),
                 "color": "#bcbdbf",
                 "linewidth": 0.5,
                 "latmax": 88, "zorder": 10}

    lat_0 = np.median(l1b.time_orbit.latitude)
    lon_0 = np.median(l1b.time_orbit.longitude)

#    if np.nanmean(l1b.time_orbit.latitude) < 0:
#        lat_0 *= -1.0

    plt.figure("Sentinel3A Orbit Quickview", figsize=(12, 12),
               facecolor="#4b4b4d")
    m = Basemap(projection='ortho', lon_0=lon_0, lat_0=lat_0, resolution='i')
    m.fillcontinents(color='#4b4b4d', lake_color='#4b4b4d', zorder=20)
    # draw parallels and meridians.
    m.drawparallels(np.arange(-90., 120., 10.), **grid_keyw)
    m.drawmeridians(np.arange(0., 420., 15.), **grid_keyw)

    x, y = m(l1b.time_orbit.longitude, l1b.time_orbit.latitude)
    m.scatter(x, y, c=l1b.surface_type.flag, s=10,
              cmap=plt.get_cmap("cool"), edgecolors="none", zorder=201)
    # m.plot(x, y, color="#003e6e", linewidth=2.0)

    plt.title("")

    m.warpimage(temp_filename, zorder=200, alpha=0.80)
    os.remove(temp_filename)
    plt.tight_layout()
    plt.show(block=False)


def s3_l1b_fullorbit_plot(l1b):

    # AWI eisblau #00ace5
    # AWI tiefblau #003e6e
    # AWI grau 1 #4b4b4d
    # AWI grau 2 #bcbdbf

    grid_keyw = {"dashes": (None, None),
                 "color": "#bcbdbf",
                 "linewidth": 0.5,
                 "latmax": 88}

    plt.figure("Sentinel-3A Orbit Quickview", facecolor="white")
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='i')
    m.fillcontinents(color='#00ace5', lake_color='#00ace5')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-90., 120., 10.), **grid_keyw)
    m.drawmeridians(np.arange(0., 420., 15.), **grid_keyw)

    x, y = m(l1b.time_orbit.longitude, l1b.time_orbit.latitude)
    m.plot(x, y, color="#003e6e", linewidth=2.0)

    plt.title("")
    plt.show(block=True)


def s3_l1b_corrections_plot(l1b):

    # from matplotlib.ticker import MultipleLocator

    n = len(l1b.correction.parameter_list)
    f, ax = plt.subplots(n, sharex=True, facecolor="white", figsize=(10, 16))
    for i in np.arange(n):
        correction, name = l1b.correction.get_parameter_by_index(i)
        ax[i].plot(correction, lw=2, color="#00ace5")
        ax[i].set_title(name)
        # ax[i].yaxis.set_minor_locator(MultipleLocator(0.02))
        ax[i].yaxis.grid(True, which='minor')
        ax[i].yaxis.set_tick_params(direction='out')
        ax[i].yaxis.set_ticks_position('left')
        ax[i].xaxis.set_ticks([])

        spines_to_remove = ["top", "right", "bottom"]
        for spine in spines_to_remove:
            ax[i].spines[spine].set_visible(False)

    plt.tight_layout()
    plt.show(block=False)


def s3_l1b_classifier_plot(l1b):

    # from matplotlib.ticker import MultipleLocator

    n = len(l1b.classifier.parameter_list)
    f, ax = plt.subplots(n, sharex=True, facecolor="white", figsize=(10, 16))
    for i in np.arange(n):
        name = l1b.classifier.parameter_list[i]
        correction = l1b.classifier.get_parameter(name)
        ax[i].plot(correction, lw=2, color="#00ace5")
        ax[i].set_title(name)
        # ax[i].yaxis.set_minor_locator(MultipleLocator(0.02))
        ax[i].yaxis.grid(True, which='minor')
        ax[i].yaxis.set_tick_params(direction='out')
        ax[i].yaxis.set_ticks_position('left')
        ax[i].xaxis.set_ticks([])

        spines_to_remove = ["top", "right", "bottom"]
        for spine in spines_to_remove:
            ax[i].spines[spine].set_visible(False)

    plt.tight_layout()
    plt.show(block=False)


def s3_l1b_waveform_plot(l1b):

    plt.figure("Sentinel-3A Waveform Quickview",
               facecolor="white", figsize=(12, 6))
    ax = plt.gca()

    image, elevation_range = align_echo_power(
        l1b.waveform.power, l1b.waveform.range, l1b.time_orbit.altitude)
    image_extent = (0, len(l1b.time_orbit.altitude),
                    np.amin(elevation_range), np.amax(elevation_range))

    im = ax.imshow(
        np.log10(image.transpose()), cmap=plt.get_cmap("magma"),
        interpolation='none', origin='lower', extent=image_extent,
        aspect=12)

#    im = ax.imshow(l1b.waveform.power.transpose(),
#                   cmap=plt.get_cmap("magma"), aspect=12,
#                   vmin=0, vmax=1000)

    ax.yaxis.set_tick_params(direction='out')
    ax.yaxis.set_ticks_position('left')
    ax.set_ylabel("Elevation (m)")
    ax.xaxis.set_ticks([])
    ax.spines["left"].set_position(("data", -20))

    spines_to_remove = ["top", "right", "bottom"]
    for spine in spines_to_remove:
        ax.spines[spine].set_visible(False)

    cb = plt.colorbar(im, orientation="vertical", shrink=0.5)
    cb.set_label("Echo Power (log)")
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

    aligned_power = np.ndarray(shape=(n_records, len(elevation_range)),
                               dtype=np.float32)
    aligned_power *= np.nan

    for i in np.arange(n_records):
        f = interp1d(elevation[i, :].flatten(), power[i, :].flatten(),
                     bounds_error=False, fill_value=np.nan)
        aligned_power[i, :] = f(elevation_range)

    return aligned_power, elevation_range


def get_s3l1bdata_argparser():
    """ Handle command line arguments """
    parser = argparse.ArgumentParser()
    # Mission id string: cryosat2, envisat, ...
    # Start month as list: [yyyy, mm]
    parser.add_argument(
        '-month', action='store', dest='month',
        nargs='+', type=str, required=True,
        help='start date as year and month (-t0 yyyy mm)')
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

    sentinel3_to_l1bdata()
