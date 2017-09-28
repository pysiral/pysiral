# -*- coding: utf-8 -*-
"""
All about mask files

Created on Thu Sep 28 14:00:52 2017

@author: shendric
"""

from pysiral.config import ConfigInfo
from pysiral.errorhandler import ErrorStatus
from pysiral.grid import GridDefinition
from pysiral.logging import DefaultLoggingClass
from pysiral.iotools import ReadNC

from collections import OrderedDict
from netCDF4 import Dataset

from pyresample import image, geometry
import numpy as np

import os


def MaskSourceFile(mask_name, mask_cfg):
    """ Wrapper method for different mask source file classes """

    error = ErrorStatus(caller_id="MaskSourceFile")

    # Get the full mask filename
    pysiral_cfg = ConfigInfo()
    try:
        mask_dir = pysiral_cfg.local_machine.auxdata_repository.mask[mask_name]
    except KeyError:
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
        self.error = ErrorStatus()
        self._cfg = cfg
        self._mask_dir = mask_dir
        self._mask_name = mask_name
        self._mask = None
        self._area_def = None
        self._post_flipud = False

    def set_mask(self, mask, area_def):
        """ Set grid definition for the mask source grid using pyresample.
        The argument area_def needs to have the attributes needed as arguments
        for pyresample.geometry.AreaDefinition """

        # Set the Mask
        self._mask = mask

        # Set the area definition
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
        resample = image.ImageContainerNearest(
            self.source_mask, self.source_area_def,
            radius_of_influence=self.cfg.resample_radius_of_influence,
            fill_value=None)
        resample_result = resample.resample(griddef.pyresample_area_def)

        # pyresample uses masked arrays -> set nan's to missing data
        target_mask = resample_result.image_data
        target_mask[np.where(target_mask.mask)] = np.nan

        # Write the mask to a netCDF file
        # (the filename will be automatically generated if not specifically
        # passed to this method
        if nc_filepath is None:
            nc_filename = "%s_%s.nc" % (self.mask_name, griddef.grid_id)
            nc_filepath = os.path.join(self.mask_dir, nc_filename)
        self.log.info("Export mask file: %s" % nc_filepath)
        self._write_netcdf(nc_filepath, griddef, target_mask)

#        import matplotlib.pyplot as plt
#
#        plt.figure("source", dpi=300)
#        plt.imshow(self.source_mask, interpolation="none")
#
#        plt.figure("target", dpi=300)
#        plt.imshow(target_mask, interpolation="none")
#        plt.show()
#        stop

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
    """ A land/sea mask based on a binary file on a 2 minute grid """

    def __init__(self, mask_dir, mask_name, cfg):
        super(MaskLandSea2Min, self).__init__(mask_dir, mask_name, cfg)
        self.construct_source_mask()

    def construct_source_mask(self):
        stop


class MaskW99Valid(MaskSourceBase):
    """ A valid mask for the Warren climatology  """

    def __init__(self, mask_dir, mask_name, cfg):
        super(MaskW99Valid, self).__init__(mask_dir, mask_name, cfg)

        # Read the data and transfer the ice_mask (1: valid, 0: invalid)
        mask_source_filename = os.path.join(mask_dir, cfg.filename)
        content = ReadNC(mask_source_filename)
        mask = content.ice_mask
        mask = np.flipud(mask)

        # Set the mask (pyresample area definition from config file)
        self.set_mask(mask, self.cfg.area_def)
