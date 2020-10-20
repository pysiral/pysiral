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

from pysiral.auxdata import AuxdataBaseClass, GridTrackInterpol
from pysiral.iotools import ReadNC


class NSIDCRegionMask(AuxdataBaseClass):
    """ Provides region codes from NSIDC style region grids """

    def __init__(self, *args, **kwargs):

        super(NSIDCRegionMask, self).__init__(*args, **kwargs)

        # The region mask is static, parse the file during init
        self._data = ReadNC(self.cfg.filename)

    def get_l2_track_vars(self, l2):

        # Extract from grid
        griddef = self.cfg.options[l2.hemisphere]
        grid_lons, grid_lats = self._data.longitude, self._data.latitude
        grid2track = GridTrackInterpol(l2.track.longitude, l2.track.latitude, grid_lons, grid_lats, griddef)
        region_code = grid2track.get_from_grid_variable(self._data.region_id, flipud=False)

        # Register the variable
        self.register_auxvar("reg_code", "region_code", region_code, None)


class CCIAntarcticSeas(AuxdataBaseClass):
    """
    Provides region codes for the *southern hemisphere only* based on longitude ranges
    from Antarctic Seas (source: Stefanie Arndt, AWI, pers comm).
    """

    def __init__(self, *args, **kwargs):
        """
        Init the class
        :param args:
        :param kwargs:
        """
        super(CCIAntarcticSeas, self).__init__(*args, **kwargs)

    def get_l2_track_vars(self, l2):
        """
        API method, will return the region cose based on longitude and latitude values.
        The location of the sea, the corresponding code and names are defined in the
        auxiliary data catalog `auxdata_def.yaml` in the pysiral resources directory.

        The expected structure of the options is:

        options:
            inland_sea_or_lake_code: 0
            ice_free_ocean_code: 1
            land_code: 20
            region_def:
                # To be read as: [region code, region label, lon_min, lon_max, lat_limit]
                - [1, "Indian Ocean", 20.0,  90.0, -50.]
                - ...

        :param l2:
        :return:
        """

        # Get a default region code value (ice free ocean)
        # -> Land will be explicitely set
        # -> Definition of seas should be complete
        default_code = self.cfg.options.get("ice_free_ocean_code", -1)
        if not default_code:
            logger.warning("Missing options: ice_free_ocean_code")
        region_code = np.full(l2.n_records, default_code).astype(int)

        # Step 1: Transfer Land/inland ice/lake flags from surface type flag to region code
        from_surface_types = self.cfg.options.get("from_surface_types", [])
        if not from_surface_types:
            logger.warning("No code for land/land ice")
        for surface_type, code in from_surface_types:
            flag = l2.surface_type.get_by_name(surface_type)
            code = self.cfg.options.get(surface_type, -1)
            region_code[flag.indices] = code

        # Step 2: Identify and set the Antarctic seas
        region_def = self.cfg.options.get("region_def", [])
        if not from_surface_types:
            logger.error("No definition of polar seas")

        # Loop over all regions
        for region_definition in region_def:

            # Expand information
            code, name, lon_min, lon_max, lat_limit = region_definition

            # Compute angle of between region longitudes
            region_angle_range_deg = self.get_lon_angular_distance(lon_max, lon_min)

            # Compute angles between data and two limits

            # Data is in region if sum of angles to both limits equals the angular
            # distance between the region limits
            

    @staticmethod
    def get_lon_angular_distance(lon1, lon2):
        """
        Return the angular distance between two longitude values that
        follow come in units of degrees east and degrees west
        :param lon1:
        :param lon2:
        :return: angular distance
        """






        # Register the variable
        self.register_auxvar("reg_code", "region_code", region_code, None)


class BlankRegionMask(AuxdataBaseClass):
    """ A dummy region code class """

    def __init__(self, *args, **kwargs):
        super(BlankRegionMask, self).__init__(*args, **kwargs)
        pass

    def get_l2_track_vars(self, l2):
        # Just use a dummy array
        region_code = np.full(l2.track.longitude.shape, -1)
        # Register the variable
        self.register_auxvar("reg_code", "region_code", region_code, None)
