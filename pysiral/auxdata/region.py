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
