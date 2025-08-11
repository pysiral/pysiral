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
import scipy.ndimage as ndimage

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
