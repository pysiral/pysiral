# -*- coding: utf-8 -*-
"""

Important Note:

    All region data handlers must be subclasses of pysiral.auxdata.AuxdataBaseClass in order to work
    for the Level-2 Processor. If the auxiliary class is based on a static dataset, this should be parsed
    in `__init__`.

    Please review the variables and properties in the parent class, as well as the corresponding config and
    support classes for grid track interpolation in the pysiral.auxdata module for additional guidance.

    The only other hard requirements is the presence of on specific method in order to be a valid subclass of
    AuxdataBaseClass:


        get_l2_track_vars(l2)

            This method will be called during the Level-2 processor. The argument is the Level-2 data object and
            the purpose of the method is to compute the auxilary variable(s) and associated uncertainty. These
            variable need to be registered using the `register_auxvar(id, name, value, uncertainty)` method of
            the base class. All MSS subclasses need to register at minimum the following variable:

            region code (integer):
                id: reg_code
                name: region_code

            e.g., this code line is mandatory for `get_l2_track_vars` (uncertainty is None):

                # Register Variables
                self.register_auxvar("reg_code", "region_code", value, None)

"""

import numpy as np
from loguru import logger
from typing import Union
from xarray import open_dataset

from pysiral.l2data import Level2Data
from pysiral.auxdata import AuxdataBaseClass, GridTrackInterpol


class NSIDCRegionMask(AuxdataBaseClass):
    """ Provides region codes from NSIDC style region grids """

    def __init__(self, *args, **kwargs) -> None:

        super(NSIDCRegionMask, self).__init__(*args, **kwargs)

        # The region mask is static, parse the file during init
        self.nc = open_dataset(self.cfg.filename)

    def get_l2_track_vars(self, l2: "Level2Data") -> None:
        """
        API method to map gridded region id on the trajectory
        :param l2:
        :return:
        """

        # Extract from grid
        griddef = self.cfg.options[l2.hemisphere]
        grid_lons, grid_lats = self.nc.longitude.values, self.nc.latitude.values
        grid2track = GridTrackInterpol(l2.track.longitude, l2.track.latitude, grid_lons, grid_lats, griddef)
        region_code = grid2track.get_from_grid_variable(self.nc.region_id.values, flipud=False)

        # Register the variable
        self.register_auxvar("reg_code", "region_code", region_code, None)


class AntarcticSeas(AuxdataBaseClass):
    """
    Provides region codes for the *southern hemisphere only* based on longitude ranges
    from Antarctic Seas (source: Stefanie Arndt, AWI, pers comm).
    """

    def __init__(self, *args, **kwargs) -> None:
        """
        Init the class
        :param args:
        :param kwargs:
        """
        super(AntarcticSeas, self).__init__(*args, **kwargs)

    def get_l2_track_vars(self, l2: "Level2Data") -> None:
        """
        API method, will return the region cose based on longitude and latitude values.
        The location of the sea, the corresponding code and names are defined in the
        auxiliary data catalog `auxdata_def.yaml` in the pysiral resources directory.

        The expected structure of the options is:

        options:
            ice_free_ocean_code: 1
            region_def:
                # To be read as: [region code, region label, lon_min, lon_max, lat_limit]
                - [1, "Indian Ocean", 20.0,  90.0, -50.]
                - ...

        :param l2:
        :return:
        """

        # Get a default region code value (ice free ocean)
        # -> Land will be explicitly set
        # -> Definition of seas should be complete
        default_code = self.cfg.options.get("ice_free_ocean_code", -1)
        if default_code == -1:
            logger.warning("- Missing options: ice_free_ocean_code")
        region_code = np.full(l2.n_records, default_code).astype(int)

        # Verify latitudes are actually in the southern hemisphere
        if (l2.latitude > 0).any():
            msg = "- Northern hemisphere positions detected, region code will be erronous"
            logger.error(msg)

        # Identify and set the Antarctic seas
        region_def = self.cfg.options.get("region_def", [])
        if not region_def:
            logger.error("- No definition of polar seas")

        # Loop over all regions
        for region_definition in region_def:

            # Expand information
            code, name, lon_min, lon_max, lat_limit = region_definition

            # Compute angle of between region longitudes
            area_angle_dist = self.get_angular_distance(lon_max, lon_min)

            # Compute angles between data and two limits
            lon_dist1 = self.get_angular_distance(lon_max, l2.longitude)
            lon_dist2 = self.get_angular_distance(lon_min, l2.longitude)

            # Data is in region if sum of angles to both limits equals the angular
            # distance between the region limits and is below latitude threshold
            in_lon_range = (lon_dist1 + lon_dist2 - area_angle_dist) < 1.-0e-9
            in_lat_range = l2.latitude < lat_limit
            in_region = np.logical_and(in_lon_range, in_lat_range)
            region_code[in_region] = code

        # Register the variable
        self.register_auxvar("reg_code", "region_code", region_code, None)

    @staticmethod
    def get_angular_distance(lon1: Union[float, int, np.array],
                             lon2: Union[float, int, np.array]) -> Union[float, int, np.array]:
        """
        Return the angular distance between two longitude values that
        follow come in units of degrees east and degrees west
        :param lon1:
        :param lon2:
        :return: angular distance
        """

        # Compute angle
        phi = np.abs(lon2 - lon1) % 360

        # Choose smaller of the two possible angles
        # NOTE: The statement below is equivalent to
        #           360. - phi if phi > 180. else phi
        #       but works for numbers and array alike
        is_larger_angle = phi > 180.
        return 360. * is_larger_angle + (-1) ** is_larger_angle * phi
