# -*- coding: utf-8 -*-
"""
Created on Sun Jun 11 19:24:04 2017

@author: Stefan
"""


from typing import Any, Tuple, Union

import numpy as np
import numpy.typing as npt
import scipy.ndimage as ndimage
from pyproj import Proj
from pyresample import geometry

from pysiral.core.config import get_yaml_config
from pysiral.core.legacy_classes import AttrDict, DefaultLoggingClass


class GridDefinition(DefaultLoggingClass):
    """ A container class for geospatial grids.  The main components are
    a) the projection (based on pyproj) and b) the extent and grid size.
    The class is currently designed for stereographic/equal area projection
    types. """

    def __init__(self, preset=None):
        super(GridDefinition, self).__init__(self.__class__.__name__)
        self._preset = preset
        self._metadata = {"grid_id": "n/a", "grid_tag": "n/a",
                          "hemisphere": "n/a", "resolution_tag": "n/a",
                          "name": "n/a", "netcdf_grid_description": "n/a"}
        self._griddef_filename = None
        self._proj = None
        self._proj_dict = {}
        self._extent_dict = {}

    def set_from_griddef_file(self, filename):
        """ Initialize the object with a grid definition (.yaml) file.
        Examples can be found in pysiral/settings/griddef """
        config = get_yaml_config(filename)
        for key in self._metadata.keys():
            try:
                self._metadata[key] = config[key]
            except KeyError:
                continue
        self.set_projection(**config.projection)
        self.set_extent(**config.extent)

    def set_projection(self, **kwargs):
        self._proj_dict = kwargs
        self._set_proj()

    def proj(self, longitude, latitude, **kwargs):
        projx, projy = self._proj(longitude, latitude, **kwargs)
        return projx, projy

    def grid_indices(self, longitude, latitude):
        """ Computes the grid indices the given lon/lat pairs would be sorted
        into (no clipping) """
        projx, projy = self.proj(longitude, latitude)
        extent = self.extent
        xi = np.floor((projx + extent.xsize/2.0)/extent.dx)
        yj = np.floor((projy + extent.ysize/2.0)/extent.dy)
        return xi, yj

    def get_grid_coordinates(self, mode="center"):
        """ Returns longitude/latitude points for each grid cell
        Note: mode keyword only for future use. center coordinates are
        returned by default """
#        x0, y0 = self.extent.xoff, self.extent.yoff
#        xsize, ysize = self.extent.xsize, self.extent.ysize
#        numx, numy = self.extent.numx, self.extent.numy
#        xmin, xmax = x0-(xsize/2.), x0+(xsize/2.)
#        ymin, ymax = y0-ysize/2., y0+ysize/2.
#        x = np.linspace(xmin, xmax, num=numx)
#        y = np.linspace(ymin, ymax, num=numy)
        xx, yy = np.meshgrid(self.xc, self.yc)
        lon, lat = self.proj(xx, yy, inverse=True)
        return lon, lat

    def set_extent(self, **kwargs):
        self._extent_dict = kwargs

    def _set_proj(self):
        self._proj = Proj(**self._proj_dict)

    @property
    def hemisphere(self):
        return self._metadata["hemisphere"]

    @property
    def grid_id(self):
        return self._metadata["grid_id"]

    @property
    def grid_tag(self):
        return self._metadata["grid_tag"]

    @property
    def grid_name(self):
        return self._metadata["name"]

    @property
    def resolution_tag(self):
        return self._metadata["resolution_tag"]

    @property
    def proj_dict(self):
        return dict(self._proj_dict)

    @property
    def extent(self):
        return AttrDict(self._extent_dict)

    @property
    def area_extent(self):
        xmin, ymin = -1*self.extent.xsize/2.0, -1*self.extent.ysize/2.0
        xmax, ymax = self.extent.xsize/2.0, self.extent.ysize/2.0
        return [xmin, ymin, xmax, ymax]

    @property
    def resolution(self):
        return self.extent.dx

    @property
    def pyresample_area_def(self):
        """ Returns a pyresample.geometry.AreaDefinition instance """

        area_def = None

        if self._proj is not None:
            # construct area definition
            area_def = geometry.AreaDefinition(
                    self.grid_id, self.grid_name, self.grid_id,
                    self.proj_dict, self.extent.numx, self.extent.numy,
                    self.area_extent)

        return area_def

    @property
    def xc(self):
        x0, numx, xsize, dx = (self.extent.xoff, self.extent.numx,
                               self.extent.xsize, self.extent.dx)
        xmin, xmax = x0-(xsize/2.)+dx/2., x0+(xsize/2.)-dx/2.
        return np.linspace(xmin, xmax, num=numx)

    @property
    def yc(self):
        y0, numy, ysize, dy = (self.extent.yoff, self.extent.numy,
                               self.extent.ysize, self.extent.dy)
        ymin, ymax = y0-(ysize/2.)+dy/2., y0+(ysize/2.)-dy/2.
        return np.linspace(ymin, ymax, num=numy)

    @property
    def xc_km(self):
        return self.xc/1000.

    @property
    def yc_km(self):
        return self.yc/1000.

    @property
    def netcdf_vardef(self):
        return self._metadata["netcdf_grid_description"]


class GridTrajectoryExtract(object):
    """
    Implements fast extraction of gridded data along a track using Image Interpolation.
    This class computes the track coordinates in image coordinates upon initialization
    and then allows to extract the multiple variables.

    Requirements are the longitude, latitude values of the trajectory and the projection
    and extent of the grid.

    Usage
    -----

    .. code-block:: python

        grid2track = GridTrajectoryExtract(track_longitude, track_latitude, grid_def)
        track_var_01 = grid2track.get_from_grid_variable(grid_var_01, outside_value=outside_value_01)
        track_var_02 = grid2track.get_from_grid_variable(grid_var_02, outside_value=outside_value_02)
        ...

    """

    def __init__(self,
                 trajectory_longitude: npt.NDArray,
                 trajectory_latitude: npt.NDArray,
                 griddef: Union[dict, AttrDict]
                 ) -> None:
        """
        Computes image coordinates ([0...1], [0...1]) of a trajectory with respect to a
        grid with known projection and extent. The image coordinates are stored in the
        instance and can be used for the extraction of several variables on the same
        grid.

        Grid Definition
        ---------------

        The grid definition contains the projection info as a dictiopnary that can be passed
        to pyproj (`pyproj.Proj(**grid_definition.projection)`) and the dimension dictionary
        that contains information, from which the grid extent can be computed, e.g.:

           projection: {'proj': stere, 'ellps': WGS84, 'lon_0': 0, 'lat_0': -90, 'lat_ts': -70,
           'a': 6378273, 'b': 6356889.44891}

           dimension: {'n_cols': 632, 'n_lines': 664, 'dx': 12500, 'dy': 12500}

           griddef = {'projection': projection, 'dimension': dimension}

        :param trajectory_longitude: longitude of the trajectory
        :param trajectory_latitude: latitude of the trajectory
        :param griddef: grid definitions (see above)

        :raises None:

        :return: None

        """

        # Save the arguments
        self.trajectory_longitude = trajectory_longitude
        self.trajectory_latitude = trajectory_latitude
        self.griddef = AttrDict(**griddef)

        # Compute image coordinates
        self.p = Proj(**self.griddef.projection)
        self.ix, self.iy = self._get_track_image_coordinates()

    def _get_track_image_coordinates(self) -> Tuple[npt.NDArray, npt.NDArray]:
        """
        Computes and returns the image coordinates, by converting the
        trajectory lon/lat values to projection coordinates and
        scaling them tho the grid extent

        :raises: None

        :return: image x coordindate, image y coordinate
        """
        tr_x, tr_y = self.p(self.trajectory_longitude, self.trajectory_latitude)
        dim = self.griddef.dimension
        x_min, y_max = -0.5 * dim.dx * (dim.n_cols - 1), 0.5 * dim.dy * (dim.n_lines - 1)
        return (tr_x - x_min) / dim.dx, (y_max - tr_y) / dim.dy

    def get_from_grid_variable(self,
                               grid_var: npt.NDArray,
                               order: int = 0,
                               flipud: bool = False,
                               outside_value: Any = np.nan
                               ) -> npt.NDArray:
        """
        Returns the grid variable along the trajectory using interpolation
        with the specified order. If required the grid var can be flipped
        upside down before extraction (flipud = True).

        :param grid_var: The gridded variable
        :param order: Order of the image interpolation (see scipy.ndimage.map_coordinates)
        :param flipud: Flag if the grid variable should be flipped before extraction.
        :param outside_value: Value to be used for part of the trajectory outside the grid

        :raises: None

        :return: The grid variable extracted and interpolation for the trajectory location
        """
        grid_var = np.flipud(grid_var) if flipud else grid_var
        return ndimage.map_coordinates(grid_var,
                                       [self.iy, self.ix],
                                       order=order,
                                       mode="constant",
                                       cval=outside_value
                                       )

# TODO: Add automatic debug map
# import matplotlib.pyplot as plt
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature
#
# fig = plt.figure(dpi=150)
# cmap_kwargs = {
#     "cmap": plt.get_cmap("plasma"),
#     "vmin": 0,
#     "vmax": 20
# }
# proj = ccrs.LambertAzimuthalEqualArea(central_latitude=90)
# ax = fig.add_subplot(1, 1, 1, projection=proj)
# ax.add_feature(cfeature.COASTLINE)
# ax.add_feature(cfeature.OCEAN)
# ax.add_feature(cfeature.LAND)
# ax.imshow(self.nc.region_id.values,
#           transform=proj,
#           extent=[-5400000, 5400000, -5400000, 5400000],
#           origin="lower",
#           **cmap_kwargs)
# ax.scatter(
#     l2.longitude,
#     l2.latitude,
#     c=region_code,
#     edgecolors="white",
#     transform=ccrs.PlateCarree(),
#     **cmap_kwargs
# )
# ax.set_extent([-180, 180, 45, 90], ccrs.PlateCarree())
# plt.show()
