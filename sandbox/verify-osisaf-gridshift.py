# -*- coding: utf-8 -*-
"""
Created on Fri Nov 04 16:31:40 2016

@author: shendric
"""

from pysiral.config import ConfigInfo
from pysiral.iotools import ReadNC
from pysiral.l1bdata import L1bdataNCFile

import matplotlib.pyplot as plt

import numpy as np

import scipy.ndimage as ndimage
from pyproj import Proj

l1b_file = r"testdata\osisaf-grid-issue\l1bdata_030dev_cryosat2_north_20160301T011804_20160301T012807.nc"
ncfile_shifted = r"testdata\osisaf-grid-issue\ice_conc_nh_polstere-100_multi_201609191200.nc"
ncfile_correct = r"testdata\osisaf-grid-issue\ice_conc_nh_polstere-100_multi_201609191200_corrected_grid.nc"


def grid_parameter_differences():

    grids_shifted = ReadNC(ncfile_shifted)
    grids_correct = ReadNC(ncfile_correct)

    for parameter in ["ice_conc", "lon", "lat"]:
        plt.figure(parameter)
        a, b = getattr(grids_shifted, parameter), getattr(grids_correct, parameter)
        print parameter, a.shape, b.shape
        if len(a.shape) == 3:
            a, b = a[0, :, :], b[0, :, :]
        plt.imshow(a-b)

    for parameter in ["xc", "yc"]:
        plt.figure(parameter)
        a, b = getattr(grids_shifted, parameter), getattr(grids_correct, parameter)
        plt.plot(a-b)

    plt.show()

def verify_osisaf_gridshift():

    orbit = L1bdataNCFile(l1b_file)
    orbit.parse()

    sic_shifted, sic_shifted_track = get_along_track_sic(ncfile_shifted, orbit, True)
    sic_correct, sic_correct_track = get_along_track_sic(ncfile_correct, orbit, False)

    plt.figure("sic grid shift effect")
    plt.imshow(np.roll(sic_shifted, -2, axis=1) - sic_correct, cmap=plt.get_cmap("bwr"),
               vmin=-10, vmax=10, interpolation="none")

#    plt.figure()
#    plt.imshow(sic_shifted, alpha=0.25, interpolation="none")
#    plt.imshow(sic_correct, alpha=0.25, interpolation="none")
#
#    plt.figure("along-track sic difference")
#    plt.plot(sic_shifted_track-sic_correct_track)
#    plt.ylim(-10, 10)

    plt.show()


def get_along_track_sic(ncfile, orbit, is_shifted):

    # Read data
    data = ReadNC(ncfile)

    # Get the necessary parameter
    config = ConfigInfo()
    projection = config.auxdata.sic.osisaf.options.north.projection
    dim = config.auxdata.sic.osisaf.options.north.dimension

    # prjection
    p = Proj(**projection)

    x, y = p(data.lon, data.lat)
    l2x, l2y = p(orbit.time_orbit.longitude, orbit.time_orbit.latitude)

    if is_shifted:
        x -= 20000.

    # Convert track projection coordinates to image coordinates
    # x: 0 < n_lines; y: 0 < n_cols
    x_min = x[dim.n_lines-1, 0]
    y_min = y[dim.n_lines-1, 0]
    ix, iy = (l2x-x_min)/dim.dx, (l2y-y_min)/dim.dy

    # Extract along track data from grid
    sic_orbit = ndimage.map_coordinates(data.ice_conc[0, :, :], [iy, ix],
                                        order=0)

    return data.ice_conc[0, :, :], sic_orbit


if __name__ == "__main__":
    # grid_parameter_differences()
    verify_osisaf_gridshift()
