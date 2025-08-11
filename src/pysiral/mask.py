# -*- coding: utf-8 -*-
"""
All about mask files

Created on Thu Sep 28 14:00:52 2017

@author: shendric
"""

import contextlib
import struct
from collections import OrderedDict
from pathlib import Path
from typing import Tuple

import numpy as np
import numpy.typing as npt
import xarray as xr
from loguru import logger
from netCDF4 import Dataset
from pyresample import geometry, image, kd_tree

from src.pysiral import psrlcfg
from core.class_template import DefaultLoggingClass
from core.errorhandler import ErrorStatus
from core.flags import SURFACE_TYPE_DICT
from core.iotools import ReadNC
from grid import GridDefinition, GridTrajectoryExtract
from l1data import Level1bData
from l1preproc.procitems import L1PProcItem


def MaskSourceFile(mask_name, mask_cfg):
    """ Wrapper method for different mask source file classes """

    error = ErrorStatus(caller_id="MaskSourceFile")

    try:
        mask_dir = psrlcfg.local_machine.auxdata_repository.mask[mask_name]
    except KeyError:
        mask_dir = None
        msg = f"path to mask {mask_name} not in local_machine_def.yaml"
        error.add_error("missing-lmd-def", msg)
        error.raise_on_error()

    # Return the Dataset class
    try:
        return globals()[mask_cfg.pyclass_name](mask_dir, mask_name, mask_cfg)
    except KeyError:
        msg = f"pysiral.mask.{str(mask_cfg.pyclass_name)} not implemented"
        error.add_error("missing-mask-class", msg)
        error.raise_on_error()


class MaskSourceBase(DefaultLoggingClass):
    """ Parent class for various source masks. Main functionality is to
    create gridded mask netCDF for level-3 grid definitions """

    def __init__(self, mask_dir, mask_name, cfg):
        super(MaskSourceBase, self).__init__(self.__class__.__name__)
        self._cfg = cfg
        self._mask_dir = mask_dir
        self._mask_name = mask_name
        self._mask = None
        self._area_def = None
        self._post_flipud = False

    def set_mask(self, mask, area_def):
        """ Set grid definition for the mask source grid using pyresample.
        The argument area_def needs to have the attributes needed as arguments
        for pyresample.geometry.AreaDefinition or be of the pyresampe types
        (geometry.AreaDefinition, geometry.GridDefinition) """

        # Set the Mask
        self._mask = mask

        # Set the area definition
        pyresample_instances = (geometry.AreaDefinition, geometry.GridDefinition)
        if isinstance(area_def, pyresample_instances):
            self._area_def = area_def
        else:
            self._area_def = geometry.AreaDefinition(
                    area_def.area_id, area_def.name, area_def.proj_id,
                    dict(area_def.proj_dict), area_def.x_size, area_def.y_size,
                    area_def.area_extent)

    def export_l3_mask(self, griddef, nc_filepath=None):
        """ Create a gridded mask product in pysiral compliant filenaming.
        The argument griddef is needs to be a pysiral.grid.GridDefinition
        instance """

        # Get the area definition for the grid
        if not isinstance(griddef, GridDefinition):
            msg = "griddef needs to be of type pysiral.grid.GridDefinition"
            self.error.add_error("value-error", msg)

        # Resample the mask
        if self.cfg.pyresample_method == "ImageContainerNearest":
            resample = image.ImageContainerNearest(
                    self.source_mask, self.source_area_def,
                    **self.cfg.pyresample_keyw)
            resample_result = resample.resample(griddef.pyresample_area_def)
            target_mask = resample_result.image_data

        elif self.cfg.pyresample_method == "resample_gauss":

            result, stddev, count = kd_tree.resample_gauss(
                    self.source_area_def, self.source_mask,
                    griddef.pyresample_area_def, with_uncert=True,
                    **self.cfg.pyresample_keyw)
            target_mask = result

        else:
            msg = f"Unrecognized opt pyresample_method: {str(self.cfg.pyresample_method)} need to be (" \
                  f"ImageContainerNearest, resample_gauss) "

            self.error.add_error("invalid-pr-method", msg)
            self.error.add_error()

        # pyresample may use masked arrays -> set nan's to missing data
        with contextlib.suppress(AttributeError):
            target_mask[np.where(target_mask.mask)] = np.nan

        if "post_processing" in self.cfg:
            pp_method = getattr(self, self.cfg.post_processing)
            target_mask = pp_method(target_mask, griddef)

        # Write the mask to a netCDF file
        # (the filename will be automatically generated if not specifically
        # passed to this method
        if nc_filepath is None:
            nc_filename = f"{self.mask_name}_{griddef.grid_id}.nc"
            nc_filepath = Path(self.mask_dir) / nc_filename
        logger.info(f"Export mask file: {nc_filepath}")
        self._write_netcdf(nc_filepath, griddef, target_mask)

    def _write_netcdf(self, nc_filepath, griddef, mask):
        """ Write a netCDF file with the mask in the target
        grid projections"""

        # Get metadata
        shape = np.shape(mask)
        dimdict = OrderedDict([("x", shape[0]), ("y", shape[1])])

        # Get longitude/latitude from grid definition
        lons, lats = griddef.get_grid_coordinates()

        # Open the file
        try:
            rootgrp = Dataset(nc_filepath, "w")
        except RuntimeError:
            rootgrp = None
            msg = f"Unable to create netCDF file: {nc_filepath}"
            self.error.add_error("nc-runtime", msg)
            self.error.raise_on_error()

        # Write Global Attributes
        rootgrp.setncattr("title", "Mask file for pysiral Level3 Processor")
        rootgrp.setncattr("mask_id", self.mask_name)
        rootgrp.setncattr("comment", self.cfg.comment)
        rootgrp.setncattr("description", self.cfg.label)
        rootgrp.setncattr("grid_id", griddef.grid_id)

        # Write dimensions
        dims = dimdict.keys()
        for key in dims:
            rootgrp.createDimension(key, dimdict[key])

        # Write Variables
        dim = tuple(dims[:len(mask.shape)])
        dtype_str = mask.dtype.str
        varmask = rootgrp.createVariable("mask", dtype_str, dim, zlib=True)
        varmask[:] = mask

        dtype_str = lons.dtype.str
        varlon = rootgrp.createVariable("longitude", dtype_str, dim, zlib=True)
        setattr(varlon, "long_name", "longitude of grid cell center")
        setattr(varlon, "standard_name", "longitude")
        setattr(varlon, "units", "degrees")
        setattr(varlon, "scale_factor", 1.0)
        setattr(varlon, "add_offset", 0.0)
        varlon[:] = lons

        varlat = rootgrp.createVariable("latitude", dtype_str, dim, zlib=True)
        setattr(varlat, "long_name", "latitude of grid cell center")
        setattr(varlat, "standard_name", "latitude")
        setattr(varlat, "units", "degrees")
        setattr(varlat, "scale_factor", 1.0)
        setattr(varlat, "add_offset", 0.0)
        varlat[:] = lats

        # Close the file
        rootgrp.close()

    @property
    def cfg(self):
        return self._cfg

    @property
    def mask_name(self):
        return str(self._mask_name)

    @property
    def mask_dir(self):
        return str(self._mask_dir)

    @property
    def source_mask(self):
        return self._mask

    @property
    def source_area_def(self):
        return self._area_def


class MaskLandSea2Min(MaskSourceBase):
    """ A land/sea mask based on a binary file on a 2 minute grid.
    Content of orignial mask: (0: sea, 1: lakes, 2: land ice: 3: land)
    There seems to be a few issues with the land ice mask in some places,
    therefore we limit the mask to 0: sea, 1: mixed, 2: non-sea (land) """

    def __init__(self, mask_dir, mask_name, cfg):
        super(MaskLandSea2Min, self).__init__(mask_dir, mask_name, cfg)
        self.construct_source_mask()

    def construct_source_mask(self):
        """ Read the binary file and set the mask """

        # Settings for the binary file
        xdim, ydim = 10800, 5400
        mask_struct_fmt = "<%0.fB" % (xdim * ydim)
        n_bytes_header = 1392

        # Read the content of the landmask in a string
        with open(str(self.mask_filepath), "rb") as fh:
            # Skip header
            fh.seek(n_bytes_header)
            content = fh.read(xdim*ydim)

        # decode string & order to array
        mask_val = np.array(struct.unpack(mask_struct_fmt, content))
        mask = mask_val.reshape((xdim, ydim))
        mask = mask.transpose()

        # Convert to only land/sea flag
        # Note: mask must be a byte data type since netCDF does not handle
        #       variables of type bool very well
        mask = np.int8(mask > 0)

        # Compute longitude/latitude grids
        lons_1d = np.linspace(0., 360., xdim)
        lats_1d = np.linspace(-90, 90, ydim)
        lons, lats = np.meshgrid(lons_1d, lats_1d)

        # Create geometry definitions
        area_def = geometry.GridDefinition(lons=lons, lats=lats)

        # Set the mask
        self.set_mask(mask, area_def)

    @staticmethod
    def pp_classify(resampled_mask, *args):
        """ Post-processing method after resampling to target grid
        The resampled mask contains a land fraction (datatype float), which
        needs to be simplified to the flags (0: ocean, 1: mixed, 2: land) """

        # Input array is range [0:1], scale to [0:2]
        pp_mask = np.copy(resampled_mask*2)

        # Get indices that are neither 0 or 2 and set them to mixed flag
        is_not_sea = np.logical_not(np.isclose(pp_mask, 0.))
        is_not_land = np.logical_not(np.isclose(pp_mask, 2.))
        is_mixed = np.logical_and(is_not_sea, is_not_land)
        pp_mask[np.where(is_mixed)] = 1.

        # Return as byte array (also needs to be flipped)
        return np.flipud(pp_mask.astype(np.int8))

    @property
    def mask_filepath(self):
        return Path(self.mask_dir) / self.cfg.filename


class MaskW99Valid(MaskSourceBase):
    """ A valid mask for the Warren climatology  """

    def __init__(self, mask_dir, mask_name, cfg):
        super(MaskW99Valid, self).__init__(mask_dir, mask_name, cfg)

        # Read the data and transfer the ice_mask (1: valid, 0: invalid)
        content = ReadNC(self.mask_filepath)
        mask = content.ice_mask
        mask = np.flipud(mask)

        # Set the mask (pyresample area definition from config file)
        self.set_mask(mask, self.cfg.area_def)

    @staticmethod
    def pp_limit_lat(resampled_mask, griddef):
        """ There are some artefacts in the source mask that need to be
        filtered out based on a simple latitude threshold filter.
        We also set all NaN values to 0 and fix a small problem at the
        north pole """

        # Get longitude/latitude valies for target grid
        lons, lats = griddef.get_grid_coordinates()

        # Set mask to false for all grid cells south of 65N
        resampled_mask[np.where(lats <= 65.)] = 0

        # Fix north pole issue
        resampled_mask[np.where(lats >= 89.)] = 1

        # Set all NaN's to 0
        resampled_mask[np.where(np.isnan(resampled_mask))] = 0

        # Done
        return resampled_mask

    @property
    def mask_filepath(self):
        return Path(self.mask_dir) / self.cfg.filename


class L3Mask(DefaultLoggingClass):
    """ Container for Level-3 mask compliant netCDF files
    (see output of pysiral.mask.MaskSourceBase.export_l3_mask) """

    def __init__(self, mask_name, grid_id, flipud=False):
        """ Mask container for Level3Processor. Arguments are the
        name (id) of the mask (e.g. warren99_is_valid) and the id of the
        grid (e.g. nh25kmEASE2) """

        super(L3Mask, self).__init__(self.__class__.__name__)
        self.error = ErrorStatus()

        # Save input
        self._mask_name = mask_name
        self._grid_id = grid_id
        self._flipud = flipud

        # Read the mask
        self._read_mask_netcdf()

    def _read_mask_netcdf(self):
        """ Read the mask """
        if self.mask_filepath is not None:
            self._nc = ReadNC(self.mask_filepath)

    @property
    def mask(self):
        mask = self._nc.mask
        if self._flipud:
            mask = np.flipud(mask)
        return mask

    @property
    def lat(self):
        lat = self._nc.latitude
        if self._flipud:
            lat = np.flipud(lat)
        return lat

    @property
    def lon(self):
        lon = self._nc.longitude
        if self._flipud:
            lon = np.flipud(lon)
        return lon

    @property
    def mask_name(self):
        return str(self._mask_name)

    @property
    def grid_id(self):
        return str(self._grid_id)

    @property
    def mask_filepath(self):

        # Get the path to the mask file
        # (needs to be in local_machine_def.yaml)
        mask_dir = psrlcfg.local_machine.auxdata_repository.mask
        try:
            mask_dir = mask_dir[self.mask_name]
        except KeyError:
            msg = "cannot find mask entry [%s] in local_machine_def.yaml"
            self.error.add_error("lmd-error", msg % self.mask_name)
            return None

        mask_filename = f"{self.mask_name}_{self.grid_id}.nc"
        filepath = Path(mask_dir) / mask_filename

        if not filepath.is_file():
            msg = f"cannot find mask file: {filepath}"
            self.error.add_error("io-error", msg)
            return None

        return filepath


class L1PHighResolutionLandMask(L1PProcItem):
    """
    Level-1 processor item providing access to a high resolution
    land mask and distance to land fields
    """

    def __init__(self, **cfg):
        """
        Initialize the class. This step includes parsing the static mask
        and keeping it in memory
        :param cfg:
        """
        super(L1PHighResolutionLandMask, self).__init__(**cfg)

        # Get the local file path
        type_ = self.cfg.get("local_machine_def_auxclass")
        tag = self.cfg.get("local_machine_def_tag")
        filename = self.cfg.get("filename")
        lookup_directory = psrlcfg.local_machine.auxdata_repository[type_][tag]

        # Set the file for each hemisphere type
        mask_filepath = Path(lookup_directory) / str(filename)
        with xr.open_dataset(mask_filepath) as nc:
            self.grid_def = self.cfg.get("grid_def")
            self.land_ocean_flag_grid = nc.land_ocean_flag.values
            self.distance_to_coast_grid = nc.distance_to_coast.values

    def apply(self, l1: Level1bData) -> None:
        """
        Extract land/ocean flag and distance to coast along the l1p trajectory if a mask exists
        for the corresponding hemisphere of the l1p data object.

        The parameters are stored in the classifier data group among the original surface type
        value, which is then update for the mask coverage

        :param l1:
        :return: None
        """

        # Get the land/ocean flag from the grid
        land_ocean_flag, distance_to_coast = self.get_trajectory(l1.time_orbit.longitude, l1.time_orbit.latitude)

        # TODO: Update when global land/ocean flag becomes available
        # NOTE: Currently only a grid for the northern hemisphere is available, thus
        #       there will be regions not covered by the gridded land/ocean flag.
        #       That's why the built/in land/ocean flag will only be updated for
        #       region with coverage and be left untouched everywhere else.
        #       The arrays are added to the l1 classifier container regardless
        #       for consistency.

        # --- Update the L1 data container ---
        # 1. Save both extracted variables to classifier data groups
        l1.classifier.add(land_ocean_flag, "hr_land_ocean_flag")
        l1.classifier.add(distance_to_coast, "distance_to_coast")

        # 2. Save original ESA surface type variable to classifier data group
        #    and update the surface type flag in the surface type data group
        l1.classifier.add(l1.surface_type.flag, "orig_land_ocean_flag")

        # 3. Update the surface type instance
        valid_mask_indices = land_ocean_flag != self.cfg.get("dummy_val")["land_ocean_flag"]
        flag_update = np.full(l1.n_records, SURFACE_TYPE_DICT["invalid"])
        flag_update[land_ocean_flag == 1] = SURFACE_TYPE_DICT["land"]
        flag_update[land_ocean_flag == 0] = SURFACE_TYPE_DICT["ocean"]
        updated_surface_type_flag = l1.surface_type.flag.copy()
        updated_surface_type_flag[valid_mask_indices] = flag_update[valid_mask_indices]
        l1.surface_type.set_flag(updated_surface_type_flag)

    def get_trajectory(self, longitude: npt.NDArray, latitude: npt.NDArray) -> Tuple[npt.NDArray, npt.NDArray]:
        """
        Get extract of the land/ocean flag and the distance to coast value for an
        array of (longitude, latitude) positions.

        If longitude, latitude is outside the grid, a pre-defined dummy value will
        be returned.

        :param longitude: Longitude values in degrees
        :param latitude: latitude values in degrees

        :raises None:

        :return: land ocean flag & distance to coast values for longitude, latitude positions
        """
        grid2track = GridTrajectoryExtract(longitude, latitude, self.grid_def)
        land_ocean_flag = grid2track.get_from_grid_variable(
            self.land_ocean_flag_grid,
            outside_value=self.dummy_val["land_ocean_flag"]
        )

        distance_to_coast = grid2track.get_from_grid_variable(
            self.distance_to_coast_grid,
            outside_value=self.dummy_val["distance_to_coast"]
        )
        return land_ocean_flag, distance_to_coast
