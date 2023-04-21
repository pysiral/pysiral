# -*- coding: utf-8 -*-
"""

module for Mean Dynamic Topography datasets.

Important Note:

    All mdt data handlers must be subclasses of pysiral.auxdata.AuxdataBaseClass in order to work
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
            the base class. All MDT subclasses need to register at minimum the following variable:

            mean dynamic topography (relative to MSS):
                id: mdt
                name: mean_dynamic_topography

            e.g., this code line is mandatory for `get_l2_track_vars` (uncertainty can be None):

                # Register Variables
                self.register_auxvar("mdt", "mean_dynamic_topography", value, uncertainty)

"""

import numpy as np
import scipy.ndimage as ndimage
from xarray import open_dataset

from pysiral.auxdata import AuxdataBaseClass

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"


class DTUMDTGrid(AuxdataBaseClass):
    """
    Parsing Routine for DTU 1/8th degree mean dynamic topography grid
    """

    def __init__(self, *args, **kwargs):
        """
        Initialize the DTU MDT auxiliary data class by reading the expected
        latitude range subset into memory.
        :param args: AuxdataBaseClass arguments
        :param kwargs: AuxdataBaseClass keyword arguments
        """

        super(DTUMDTGrid, self).__init__(*args, **kwargs)

        # Read as standard netcdf
        dtu_grid = open_dataset(self.cfg.filename)

        # Cut to ROI regions (latitude only)
        # -> no need for world mdt
        lat_range = self.cfg.options.latitude_range

        # Get the indices for the latitude subset
        latitude_indices = np.where(np.logical_and(dtu_grid.lat.values >= lat_range[0],
                                                   dtu_grid.lat.values <= lat_range[1]))[0]

        # Crop data to subset
        if 'variables' in self.cfg.options:
            self.mean_dynamic_topography = dtu_grid.__getattr__(self.cfg.options['variables'][0]).values[latitude_indices, :]
            self.mean_dynamic_topography[np.isnan(self.mean_dynamic_topography)] = 0.0
            if self.cfg.options['variables'][1] is not None:
                self.mean_dynamic_topography_uncertainty = dtu_grid.__getattr__(self.cfg.options['variables'][1]).values[latitude_indices, :]
            else:
                self.mean_dynamic_topography_uncertainty = np.zeros_like(self.mean_dynamic_topography)
        else:
            self.mean_dynamic_topography = dtu_grid.MDT.values[latitude_indices, :]
            self.mean_dynamic_topography[np.isnan(self.mean_dynamic_topography)] = 0.0
            self.mean_dynamic_topography_uncertainty = dtu_grid.err.values[latitude_indices, :]
            self.mean_dynamic_topography_uncertainty[np.isnan(self.mean_dynamic_topography_uncertainty)] = 0.0

        self.longitude = dtu_grid.lon.values
        self.latitude = dtu_grid.lat.values[latitude_indices]

    def get_l2_track_vars(self, l2):

        # Use fast image interpolation (since DTU is on regular grid)
        # Longitudes must be 0 -> 360

        longitude = np.array(l2.track.longitude)
        latitude = np.array(l2.track.latitude)

        negative_lons = np.where(longitude < 0)[0]
        longitude[negative_lons] = longitude[negative_lons] + 360.

        # Calculate image coordinates of mss grid "image"
        mdt_lon_min = self.longitude[0]
        mdt_lon_step = self.longitude[1] - self.longitude[0]
        mdt_lat_min = self.latitude[0]
        mdt_lat_step = self.latitude[1] - self.latitude[0]
        ix = (longitude - mdt_lon_min)/mdt_lon_step
        iy = (latitude - mdt_lat_min)/mdt_lat_step

        # Extract and return the elevation along the track
        mdt = ndimage.map_coordinates(self.mean_dynamic_topography, [iy, ix])
        mdt_err = ndimage.map_coordinates(self.mean_dynamic_topography_uncertainty, [iy, ix])

        # Register auxdata variable
        self.register_auxvar("mdt", "mean_dynamic_topography", mdt, mdt_err)
