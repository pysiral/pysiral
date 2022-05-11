# -*- coding: utf-8 -*-
"""
All about mask files

Created on Thu Sep 28 14:00:52 2017

@author: shendric
"""
import pyproj

from pysiral import psrlcfg
from pysiral.errorhandler import ErrorStatus
from pysiral.grid import GridDefinition
from pysiral._class_template import DefaultLoggingClass
from pysiral.iotools import ReadNC
from pysiral.l1bdata import Level1bData
from pysiral.l1preproc.procitems import L1PProcItem

from collections import OrderedDict
from netCDF4 import Dataset
from loguru import logger

from pyresample import image, geometry, kd_tree
import numpy as np
import struct
import xarray as xr
from pyproj import Proj
from pathlib import Path


def MaskSourceFile(mask_name, mask_cfg):
    """ Wrapper method for different mask source file classes """

    error = ErrorStatus(caller_id="MaskSourceFile")

    try:
        mask_dir = psrlcfg.local_machine.auxdata_repository.mask[mask_name]
    except KeyError:
        mask_dir = None
        msg = "path to mask %s not in local_machine_def.yaml" % mask_name
        error.add_error("missing-lmd-def", msg)
        error.raise_on_error()

    # Return the Dataset class
    try:
        return globals()[mask_cfg.pyclass_name](mask_dir, mask_name, mask_cfg)
    except KeyError:
        msg = "pysiral.mask.%s not implemented" % str(mask_cfg.pyclass_name)
        error.add_error("missing-mask-class", msg)
        error.raise_on_error()


class MaskSourceBase(DefaultLoggingClass):
    """ Parent class for various source masks. Main functionality is to
    create gridded mask netCDF for for level-3 grid definitions """

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
            msg = "Unrecognized opt pyresample_method: %s need to be %s" % (
                    str(self.cfg.pyresample_method),
                    "(ImageContainerNearest, resample_gauss)")
            self.error.add_error("invalid-pr-method", msg)
            self.error.add_error()

        # pyresample may use masked arrays -> set nan's to missing data
        try:
            target_mask[np.where(target_mask.mask)] = np.nan
        except AttributeError:
            pass

        if "post_processing" in self.cfg:
            pp_method = getattr(self, self.cfg.post_processing)
            target_mask = pp_method(target_mask, griddef)

        # Write the mask to a netCDF file
        # (the filename will be automatically generated if not specifically
        # passed to this method
        if nc_filepath is None:
            nc_filename = "%s_%s.nc" % (self.mask_name, griddef.grid_id)
            nc_filepath = Path(self.mask_dir) / nc_filename
        logger.info("Export mask file: %s" % nc_filepath)
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
            msg = "Unable to create netCDF file: %s" % nc_filepath
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
        dim = tuple(dims[0:len(mask.shape)])
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

    def pp_classify(self, resampled_mask, griddef):
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

    def pp_limit_lat(self, resampled_mask, griddef):
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

        mask_filename = "%s_%s.nc" % (self.mask_name, self.grid_id)
        filepath = Path(mask_dir) / mask_filename

        if not filepath.is_file():
            msg = "cannot find mask file: %s" % filepath
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

        self.masks = {}
        self._init_masks()

    def _init_masks(self) -> None:
        """
        Store the masks in memory
        :return:
        """

        # Get the local file path
        type_ = self.cfg.get("local_machine_def_auxclass")
        tag = self.cfg.get("local_machine_def_tag")
        lookup_directory = psrlcfg.local_machine.auxdata_repository[type_][tag]

        # Set the file for each hemisphere type
        hemispheres = self.cfg.get("hemispheres", {})
        for hemisphere in hemispheres:
            hemisphere_cfg = self.cfg["hemispheres"][hemisphere]
            mask_filepath = Path(lookup_directory) / hemisphere_cfg["filename"]
            nc = xr.open_dataset(mask_filepath)
            self.masks[hemisphere] = {
                "projection": pyproj.Proj(nc.geospatial_bounds_crs),
                "grid_def": hemisphere_cfg["grid_def"],
                "land_ocean_flag": nc.land_ocean_flag.values,
                "distance_to_coast": nc.distance_to_coast.values
            }
            nc = None

    def apply(self, l1: Level1bData) -> None:
        """
        Extract land/ocean flag and distance to coast along the
        l1p trajectory if a mask exists for the corresponding
        hemisphere of the l1p data object.

        :param l1: 
        :return: 
        """

        # Determine the hemisphere
        if l1.info.hemisphere not in self.masks:
            logger.info(f"{self.__class__.__name__}: No mask for hemisphere {l1.info.hemisphere}")
            return

        # Get the mask val
        mask = self.masks[l1.info.hemisphere]

        # Compute the track position in image coordinates
        prj_x, prj_y = mask["projection"](l1.time_orbit.longitude, l1.time_orbit.latitude)

        breakpoint()
