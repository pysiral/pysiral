# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 16:30:24 2015

@author: Stefan
"""

from pysiral.errorhandler import ErrorStatus
from pysiral.output import PysiralOutputFilenaming
from pysiral.path import filename_from_path
from pysiral.iotools import ReadNC

import numpy as np
from geopy.distance import great_circle
from collections import OrderedDict


class Level2Data(object):

    _L2_DATA_ITEMS = ["mss", "ssa", "elev",  "afrb", "frb", "range", "sic",
                      "sitype", "snow_depth",  "snow_dens",  "ice_dens",
                      "sit"]

    _HEMISPHERE_CODES = {"north": "nh", "south": "sh"}

    _PARAMETER_CATALOG = {
        "timestamp": "timestamp",
        "longitude": "longitude",
        "latitude": "latitude",
        "surface_type": "surface_type_flag",
        "elevation": "elev",
        "mean_sea_surface": "mss",
        "sea_surface_anomaly": "ssa",
        "radar_freeboard": "afrb",
        "freeboard": "frb",
        "sea_ice_type": "sitype",
        "snow_depth": "snow_depth",
        "snow_density": "snow_dens",
        "ice_density": "ice_dens",
        "sea_ice_thickness": "sit",
        "sea_ice_concentration": "sic"}

    def __init__(self, l1b):
        # Copy necessary fields form l1b
        self.error = ErrorStatus()
        self._n_records = l1b.n_records
        self.info = l1b.info
        self.track = l1b.time_orbit
        # Create Level2 Data Groups
        self._create_l2_data_items()

    def set_surface_type(self, surface_type):
        self.surface_type = surface_type

    def set_parameter(self, target, value, uncertainty=None, bias=None):
        """ Convienience method to safely add a parameter with optional
        uncertainty and/or bias to the level-2 data structure """

        # Sanity checks
        is_valid = self._check_if_valid_parameter(target)
        is_correct_size = self._check_valid_size(value)
        if not is_valid or not is_correct_size:
            return

        # Set values, uncertainty bias
        parameter = getattr(self, target)
        parameter.set_value(value)
        if uncertainty is not None:
            uncertainty_value = self._get_as_array(uncertainty)
            parameter.set_uncertainty(uncertainty_value)
        if bias is not None:
            bias_value = self._get_as_array(bias)
            parameter.set_bias(bias, bias_value)
        setattr(self, target, parameter)

    def update_retracked_range(self, retracker):
        # Update only for indices (surface type) supplied by retracker class
        # XXX: should get an overhaul
        ii = retracker.indices
        self.range[ii] = retracker.range[ii]
        self.range.uncertainty[ii] = retracker.uncertainty[ii]
        self.elev[ii] = self.track.altitude[ii] - retracker.range[ii]
        self.elev.uncertainty[ii] = retracker.uncertainty[ii]

    def get_parameter_by_name(self, parameter_name):
        """ Method to retrieve a level-2 parameter """

        if "_uncertainty" in parameter_name:
            parameter_name = parameter_name.replace("_uncertainty", "")
            source = self._PARAMETER_CATALOG[parameter_name]
            parameter = getattr(self, source)
            return parameter.uncertainty
        elif "_bias" in parameter_name:
            parameter_name = parameter_name.replace("_bias", "")
            source = self._PARAMETER_CATALOG[parameter_name]
            parameter = getattr(self, source)
            return parameter.bias
        else:
            source = self._PARAMETER_CATALOG[parameter_name]
            parameter = getattr(self, source)
            return parameter

    def _create_l2_data_items(self):
        for item in self._L2_DATA_ITEMS:
            setattr(self, item, L2ElevationArray(shape=(self._n_records)))

    def _check_if_valid_parameter(self, parameter_name):
        """ Performs a test if parameter name is a valid level-2 parameter
        name. Adds error if result negative and returns flag (valid: True,
        invalid: False) """
        if parameter_name not in self._L2_DATA_ITEMS:
            msg = "Invalid level-2 parameter: %s" % parameter_name
            self.error.add_error("l2-invalid-parameter", msg)
            return False
        else:
            return True

    def _check_valid_size(self, array, name=""):
        """ Test if array has the correct size shape=(n_records). Adds error
        if not and returns flag (valid: True, invalid: False) """
        condition = array.ndim == 1 and len(array) == self._n_records
        if condition:
            return True
        else:
            self.error.add_error("Invalid array added to level-2 class")
            return False

    def _get_as_array(self, value, dtype=np.float32):
        """ Create an output array from values that is of length n_records.
        Value can be scalar or array of length n_records. If value is any other
        length or dimension, an error will be added and a nan array of length
        n_records will be returned

        Arguments:
            value (integer, float or )

        Note: This method is mostly used to allow scalar uncertainty and
              bias values. It also makes sure that uncertainty and bias
              are of the same shape than the value, which is not guaranteed
              in L2ElevationArray. If a wrong uncertainty, bias shape is
              passed, the result will be nan uncertainties/biases throughout
              the processing chain and the start of NaN occurences can be used
              to trace the origin of the error.
        """

        # Check if value is either float or integer
        is_numeric = np.asarray(value).dtype.kind in "if"
        if not is_numeric:
            return np.full(self.arrshape, np.nan)

        # Check if value is scalar or array
        if np.isscalar(value):
            return np.full(self.arrshape, value).astype(dtype)

        # if array, check if correct size
        else:
            is_np_array = np.isinstance(value, (np.ndarray, np.array))
            is_correct_size = self._check_valid_size(value)
            if is_np_array and is_correct_size:
                return value
            else:
                return np.full(self.arrshape, np.nan)

    @property
    def arrshape(self):
        return (self.n_records)

    @property
    def n_records(self):
        return self._n_records

    @property
    def hemisphere(self):
        return self.info.subset_region_name

    @property
    def hemisphere_code(self):
        return self._HEMISPHERE_CODES[self.hemisphere]

    @property
    def footprint_spacing(self):
        spacing = great_circle(
            (self.track.latitude[1], self.track.longitude[1]),
            (self.track.latitude[0], self.track.longitude[0])).meters
        return spacing

    @property
    def dimdict(self):
        """ Returns dictionary with dimensions"""
        dimdict = OrderedDict([("n_records", self.n_records)])
        return dimdict

    @property
    def timestamp(self):
        return self.track.timestamp

    @property
    def longitude(self):
        return self.track.longitude

    @property
    def latitude(self):
        return self.track.latitude

    @property
    def surface_type_flag(self):
        return self.surface_type.flag


class L2ElevationArray(np.ndarray):
    """
    Recipe from:
    http://docs.scipy.org/doc/numpy/user/basics.subclassing.html
    XXX: not yet full slicing capability! -> __getitem__ trouble
         always use cls[list] and cls.uncertainty[list]
         cls[list].uncertainty will fail
    """

    def __new__(subtype, shape, dtype=float, buffer=None, offset=0,
                strides=None, order=None, info=None):
        obj = np.ndarray.__new__(
            subtype, shape, dtype, buffer, offset, strides, order)*np.nan
        obj.uncertainty = np.zeros(shape=shape, dtype=float)
        obj.bias = np.ones(shape=shape, dtype=float)*0.1
        obj.source_class = "n/a"
        obj.source_files = "n/a"
        obj.long_name = "n/a"
        obj.unit = "n/a"
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.uncertainty = getattr(obj, 'uncertainty', None)
        self.bias = getattr(obj, 'bias', None)
        self.source_class = getattr(obj, 'source_class', None)
        self.source_files = getattr(obj, 'source_files', None)
        self.long_name = getattr(obj, 'long_name', None)
        self.unit = getattr(obj, 'unit', None)

    def __getslice__(self, i, j):
        r = np.ndarray.__getslice__(self, i, j)
        r.uncertainty = r.uncertainty[i:j]
        r.bias = r.bias[i:j]
        return r

    def set_value(self, value):
#        uncertainty = self.uncertainty
#        bias = self.bias
        self[:] = value[:]
#        setattr(self, "uncertainty", uncertainty)
#        setattr(self, "bias", bias)

    def set_uncertainty(self, uncertainty):
        self.uncertainty = uncertainty

    def set_bias(self, bias):
        self.bias = bias

    def set_nan_indices(self, indices):
        value = self[:]
        value[indices] = np.nan
        self.set_value(value)
        self.uncertainty[indices] = np.nan
        self.bias[indices] = np.nan


class AttributeList(object):

    def __init__(self):
        pass

    def set_attribute(self, name, value):
        setattr(self, name, value)


class L2iNCFileImport(object):
    # TODO: Needs proper implementation

    def __init__(self, filename):
        from pysiral.output import NCDateNumDef
        self.filename = filename
        self._n_records = 0
        self.time_def = NCDateNumDef()
        self.info = AttributeList()
        self.attribute_list = []
        self._parse()

    def _parse(self):
        from pysiral.path import file_basename
        from netCDF4 import num2date

        content = ReadNC(self.filename)

        for attribute_name in content.attributes:
            self.attribute_list.append(attribute_name)
            self.info.set_attribute(attribute_name,
                                    getattr(content, attribute_name))

        for parameter_name in content.parameters:
            setattr(self, parameter_name, getattr(content, parameter_name))

        self._n_records = len(self.longitude)

        # Get mission id from filename
        l2i_filename = filename_from_path(self.filename)
        filenaming = PysiralOutputFilenaming()
        filenaming.parse_filename(l2i_filename)
        self.mission = filenaming.mission_id

        self.timestamp = num2date(self.timestamp, self.time_def.units,
                                  self.time_def.calendar)

    def project(self, griddef):
        from pyproj import Proj
        p = Proj(**griddef.projection)
        self.projx, self.projy = p(self.longitude, self.latitude)
        # Convert projection coordinates to grid indices
        extent = griddef.extent
        self.xi = np.floor((self.projx + extent.xsize/2.0)/extent.dx)
        self.yj = np.floor((self.projy + extent.ysize/2.0)/extent.dy)

    @property
    def n_records(self):
        return self._n_records
