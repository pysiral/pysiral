# -*- coding: utf-8 -*-
"""
Created on Mon Mar 07 15:48:22 2016

@author: Stefan

Sandbox script to parse CryoSat-2 Baseline-B files
and export time-orbit and waveform groups to netcdf
file

"""

from pysiral.config import ConfigInfo
from pysiral.l1bdata import L1bConstructor
from pysiral.iotools import L1bNCfile
from pysiral.roi import LatitudeLongitudeBox

import os
import glob
from netCDF4 import Dataset

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np


def cryosat2_baselineb_waveform_parameter():
    # Get Configuration Information
    config = ConfigInfo()

    # Get an L1B SAR file
    l1b_directory = r"D:\awi\altim\data\altimetry\cryosat2\baseline-b\SIR_SAR_L1"
    l1b_directory = os.path.join(l1b_directory, "2013", "12")
    l1b_files = glob.glob(os.path.join(l1b_directory, "*.DBL"))

    # netcdf export
    export_datagroups = ["time_orbit", "waveform", "classifier"]
    output_folder = r"D:\awi\product\altimetry\pysiral\tests\l1b_ncfiles"

    # Region of interest settings
    roi_definition = LatitudeLongitudeBox()
    roi_definition.set_options(lat_range=[-88., -50.], lon_range=[-30., 40.])

    # loop over iput files
    for l1b_filename in l1b_files:
        # Parse L1b files
        l1b = L1bConstructor(config)
        l1b.mission = "cryosat2"
        l1b.filename = l1b_filename
        l1b.construct()
        # Check region of interest
        l1b_antarctic_subset = l1b.extract_region_of_interest(roi_definition)
        # Validate overlap
        if l1b_antarctic_subset.info.n_records == -1:
            continue
        # Write to l1b netcdf container file
        ncfile = L1bNCfile()
        ncfile.l1b = l1b_antarctic_subset
        ncfile.datagroups = export_datagroups
        ncfile.output_folder = output_folder
        ncfile.export()

    # Test plot of netcdf content
    filename = os.path.join(
        r"D:\awi\product\altimetry\pysiral\tests\l1b_ncfiles",
        r"CS_OFFL_SIR_SAR_1B_20131201T012220_20131201T012916_B001.nc")
    nc = Dataset(filename, "r")
    wfm_power = nc.groups["waveform"].variables["power"][:]
    wfm_range = nc.groups["waveform"].variables["range"][:]
    altitude = nc.groups["time_orbit"].variables["altitude"][:]
    longitude = nc.groups["time_orbit"].variables["longitude"][:]
    latitude = nc.groups["time_orbit"].variables["latitude"][:]
    nc.close()

    # Test plot
    cryosat2_l1b_waveform_plot(wfm_power, wfm_range, altitude)
    l1b_orbit_plot(longitude, latitude)


def cryosat2_l1b_waveform_plot(wfm_power, wfm_range, altitude):

    plt.figure("Waveform Quickview", facecolor="white", figsize=(12, 6))
    ax = plt.gca()

    image, elevation_range = align_echo_power(
        wfm_power, wfm_range, altitude)
    image_extent = (0, len(altitude),
                    np.amin(elevation_range), np.amax(elevation_range))

    im = ax.imshow(
        np.log(image).transpose(), cmap=plt.get_cmap("magma"),
        interpolation='none', origin='lower', extent=image_extent,
        aspect=35)

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


def l1b_orbit_plot(longitude, latitude):

    # AWI eisblau #00ace5
    # AWI tiefblau #003e6e
    # AWI grau 1 #4b4b4d
    # AWI grau 2 #bcbdbf

    grid_keyw = {"dashes": (None, None),
                 "color": "#bcbdbf",
                 "linewidth": 0.5,
                 "latmax": 88}

    lat_0 = np.median(latitude)
    lon_0 = np.median(longitude)

#    if np.nanmean(l1b.time_orbit.latitude) < 0:
#        lat_0 *= -1.0

    plt.figure("Orbit Quickview", facecolor="white")
    m = Basemap(projection='ortho', lon_0=lon_0, lat_0=lat_0, resolution='i')
    m.fillcontinents(color='#00ace5', lake_color='#00ace5')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-90., 120., 10.), **grid_keyw)
    m.drawmeridians(np.arange(0., 420., 15.), **grid_keyw)

    x, y = m(longitude, latitude)
    m.plot(x, y, color="#003e6e", linewidth=2.0)

    plt.title("")
    plt.show(block=False)


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


if __name__ == "__main__":
    cryosat2_baselineb_waveform_parameter()
