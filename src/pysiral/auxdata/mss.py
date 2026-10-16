# -*- coding: utf-8 -*-
"""

Important Note:

    All mss data handlers must be subclasses of pysiral.auxdata.AuxdataBaseClass in order to work
    for the Level-2 Processor. If the auxiliary class is based on a static dataset, this should be parsed
    in `__init__`.

    Please review the variables and properties in the parent class, as well as the correspodning config and
    support classes for grid track interpolation in the pysiral.auxdata module for additional guidance.

    The only other hard requirements is the presence of on specific method in order to be a valid subclass of
    AuxdataBaseClass:


        get_l2_track_vars(l2)

            This method will be called during the Level-2 processor. The argument is the Level-2 data object and
            the purpose of the method is to compute the auxilary variable(s) and associated uncertainty. These
            variable need to be registered using the `register_auxvar(id, name, value, uncertainty)` method of
            the base class. All MSS subclasses need to register at minimum the following variable:

            mean sea surface (wgs84 elevation in meter):
                id: mss
                name: mean_sea_surface

            e.g., this code line is mandatory for `get_l2_track_vars` (uncertainty can be None):

                # Register Variables
                self.register_auxvar("mss", "mean_sea_surface", value, uncertainty)

"""

import numpy as np
import xarray as xr
import scipy.ndimage as ndimage

from loguru import logger
from time import time

from pysiral.auxdata import AuxdataBaseClass
from pysiral.core.iotools import ReadNC


class DTU1MinGrid(AuxdataBaseClass):
    """
    Parsing Routine for DTU 1 minute global mean sea surface height files
    """

    def __init__(self, *args, **kwargs):

        super(DTU1MinGrid, self).__init__(*args, **kwargs)

        # Read as standard netcdf
        dtu_grid = ReadNC(self.cfg.filename)

        # Cut to ROI regions (latitude only)
        # -> no need for world mss
        lat_range = self.cfg.options.latitude_range

        # Get the indices for the latitude subset
        latitude_indices = np.where(np.logical_and(dtu_grid.lat >= lat_range[0], dtu_grid.lat <= lat_range[1]))[0]

        # Crop data to subset
        self.elevation = dtu_grid.mss[latitude_indices, :]
        self.longitude = dtu_grid.lon
        self.latitude = dtu_grid.lat[latitude_indices]

        # Convert elevations to WGS84
        delta_h1 = egm2top_delta_h(self.latitude)
        delta_h = egm2wgs_delta_h(self.latitude)
        for i in np.arange(len(self.latitude)):
            self.elevation[i, :] += (delta_h[i]-delta_h1[i])

    def get_l2_track_vars(self, l2):

        # Use fast image interpolation (since DTU is on regular grid)
        # Longitudes must be 0 -> 360

        longitude = np.array(l2.track.longitude)
        latitude = np.array(l2.track.latitude)

        negative_lons = np.where(longitude < 0)[0]
        longitude[negative_lons] = longitude[negative_lons] + 360.

        # Calculate image coordinates of mss grid "image"
        mss_lon_min = self.longitude[0]
        mss_lon_step = self.longitude[1] - self.longitude[0]
        mss_lat_min = self.latitude[0]
        mss_lat_step = self.latitude[1] - self.latitude[0]
        ix = (longitude - mss_lon_min)/mss_lon_step
        iy = (latitude - mss_lat_min)/mss_lat_step

        # Extract and return the elevation along the track
        mss_track_elevation = ndimage.map_coordinates(self.elevation, [iy, ix])

        # Register auxdata variable
        self.register_auxvar("mss", "mean_sea_surface", mss_track_elevation, None)


class DTU25Grid(AuxdataBaseClass):
    """
    Parsing Routine for DTU 1 km global mean sea surface height files
    """

    def __init__(self, *args, **kwargs):

        super(DTU25Grid, self).__init__(*args, **kwargs)

        # Open lazily: the (compressed) file holds two global float64 grids (~10 GB in memory),
        # of which only the mss in the latitude range of interest is needed
        t0 = time()
        with xr.open_dataset(self.cfg.filename) as dtu_grid:

            # Cut to ROI regions (latitude only)
            # -> no need for world mss
            lat_range = self.cfg.options.latitude_range

            # Get the indices for the latitude subset
            latitude_indices = np.where(
                np.logical_and(
                    dtu_grid.latitude.values >= lat_range[0],
                    dtu_grid.latitude.values <= lat_range[1]
                )
            )[0]

            # Read only the subset (latitude is monotonic -> contiguous row range)
            rows = slice(latitude_indices[0], latitude_indices[-1] + 1)
            self.elevation = dtu_grid.mss[rows, :].values
            self.longitude = dtu_grid.longitude.values  # DTU MSS25 longitude is 0-360
            self.latitude = dtu_grid.latitude.values[rows]
        t1 = time()
        logger.debug(f"Loaded DTU25 dataset in : {t1 - t0:.3f} seconds")

    def get_l2_track_vars(self, l2):

        # Use fast image interpolation (since DTU is on regular grid)
        # Longitudes must be 0 -> 360

        longitude = np.array(l2.track.longitude)
        latitude = np.array(l2.track.latitude)

        negative_lons = np.where(longitude < 0)[0]
        longitude[negative_lons] = longitude[negative_lons] + 360.

        # Calculate image coordinates of mss grid "image"
        mss_lon_min = self.longitude[0]
        mss_lon_step = self.longitude[1] - self.longitude[0]
        mss_lat_min = self.latitude[0]
        mss_lat_step = self.latitude[1] - self.latitude[0]
        ix = (longitude - mss_lon_min)/mss_lon_step
        iy = (latitude - mss_lat_min)/mss_lat_step

        # Extract and return the elevation along the track
        mss_track_elevation = interpolate_cubic_tiled(self.elevation, iy, ix)

        # Register auxdata variable
        self.register_auxvar("mss", "mean_sea_surface", mss_track_elevation, None)


def interpolate_cubic_tiled(grid, iy, ix, tile_size=256, pad=12):
    """
    Cubic spline interpolation of a regular grid at fractional image coordinates, with the
    second (longitude) axis periodic.

    Equivalent to `ndimage.map_coordinates(grid, [iy, ix])`, but the spline prefilter (which
    creates a full-size copy of the input) is only applied to tiles of the grid that contain
    points, each with a margin of `pad` cells. The influence of the prefilter decays by a
    factor of ~0.27 per cell, so the margin makes the result independent of the tiling.
    At the first and last row of the grid the spline uses mirror boundary conditions.

    :param grid: 2D grid (rows: latitude, columns: periodic longitude)
    :param iy: fractional row coordinates
    :param ix: fractional column coordinates
    :param tile_size: tile size in grid cells
    :param pad: margin around each tile in grid cells

    :return: interpolated values (NaN for non-finite coordinates)
    """
    n_rows, n_cols = grid.shape
    values = np.full(iy.shape, np.nan)
    finite = np.flatnonzero(np.isfinite(iy) & np.isfinite(ix))
    if finite.size == 0:
        return values
    iy, ix = iy[finite], np.mod(ix[finite], n_cols)

    # Group the points by tile
    tile_row = np.floor(iy / tile_size).astype(int)
    tile_col = np.floor(ix / tile_size).astype(int)
    tile_id = tile_row * (n_cols // tile_size + 1) + tile_col
    order = np.argsort(tile_id, kind="stable")
    tile_starts = np.flatnonzero(np.r_[True, np.diff(tile_id[order]) != 0])

    for start, stop in zip(tile_starts, np.r_[tile_starts[1:], order.size]):
        points = order[start:stop]

        # Window: tile + margin (rows clipped to the grid, columns wrapped)
        row0 = max(tile_row[points[0]] * tile_size - pad, 0)
        row1 = min((tile_row[points[0]] + 1) * tile_size + pad + 1, n_rows)
        col0 = tile_col[points[0]] * tile_size - pad
        columns = np.arange(col0, (tile_col[points[0]] + 1) * tile_size + pad + 1) % n_cols
        window = np.take(grid[row0:row1], columns, axis=1)

        values[finite[points]] = ndimage.map_coordinates(
            window, [iy[points] - row0, ix[points] - col0], mode="mirror"
        )

    return values


def egm2wgs_delta_h(phi):
    aegm = 6378136.460000
    begm = 6356751.806631
    awgs = 6378137.000000
    bwgs = 6356752.314245
    return compute_delta_h(aegm, begm, awgs, bwgs, phi)


def egm2top_delta_h(phi):
    aegm = 6378136.460000
    begm = 6356751.806631
    atop = 6378136.300000
    btop = 6356751.600563
    return compute_delta_h(aegm, begm, atop, btop, phi)


def compute_delta_h(a1, b1, a2, b2, phi):
    dtor = np.pi/180.0
    delta_a = a2 - a1
    delta_b = b2 - b1
    phir = phi * dtor
    sinsqphi = np.sin(phir)
    sinsqphi = sinsqphi * sinsqphi
    cossqphi = 1.0 - sinsqphi
    return -1.0*(delta_a * cossqphi + delta_b * sinsqphi)
