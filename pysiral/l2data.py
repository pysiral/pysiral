# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 16:30:24 2015

@author: Stefan
"""

from pysiral.config import (PYSIRAL_VERSION, PYSIRAL_VERSION_FILENAME,
                            SENSOR_NAME_DICT, MISSION_NAME_DICT)
from pysiral.errorhandler import ErrorStatus
from pysiral.output import PysiralOutputFilenaming
from pysiral.path import filename_from_path
from pysiral.iotools import ReadNC
from pysiral.logging import DefaultLoggingClass
from pysiral.l1bdata import L1bMetaData, L1bTimeOrbit

import numpy as np
from datetime import datetime
from geopy.distance import great_circle
from collections import OrderedDict
import re


class Level2Data(object):

    _L2_DATA_ITEMS = ["mss", "ssa", "elev",  "afrb", "frb", "range", "sic",
                      "sitype", "snow_depth",  "snow_dens",  "ice_dens",
                      "sit", "radar_mode"]

    _HEMISPHERE_CODES = {"north": "nh", "south": "sh"}

    _PARAMETER_CATALOG = {
        "timestamp": "timestamp",
        "time": "time",
        "longitude": "longitude",
        "latitude": "latitude",
        "surface_type": "surface_type_flag",
        "radar_mode": "radar_mode",
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
        "sea_ice_concentration": "sic",
        "radar_mode": "radar_mode"}

    _PROPERTY_CATALOG = {
        "sea_surface_height": "ssh"}

    def __init__(self, metadata, time_orbit, period=None):

        # Copy necessary fields form l1b
        self.error = ErrorStatus()
        self._n_records = metadata.n_records
        self.info = metadata
        self.track = time_orbit
        self.period = period

        # Metadata
        self._auxdata_source_dict = {}
        self._source_primary_filename = "unkown"
        self._l2_algorithm_id = "unkown"
        self._l2_version_tag = "unkown"

        # Define time of dataset creation as the time of object initialization
        # to avoid slightly different timestamps for repated calls of
        # datatime.now()
        self._creation_time = datetime.now()

        # Other Class properties
        self._is_evenly_spaced = time_orbit.is_evenly_spaced

        # Create Level2 Data Groups
        self._create_l2_data_items()

    def set_surface_type(self, surface_type):
        self.surface_type = surface_type

    def set_radar_mode(self, radar_mode):
        self.radar_mode = radar_mode

    def set_parameter(self, target, value, uncertainty=None, bias=None):
        """ Convienience method to safely add a parameter with optional
        uncertainty and/or bias to the level-2 data structure """

        # Sanity checks
        is_valid = self._check_if_valid_parameter(target)

        # Check if the full name has been passed
        if not is_valid and target in self._PARAMETER_CATALOG.keys():
            target = self._PARAMETER_CATALOG[target]
            is_valid = True

        # Next check: Needs to be of correct shape
        is_correct_size = self._check_valid_size(value)
        if not is_valid or not is_correct_size:
            msg = "Invalid parameter name: %s (See self._L2_DATA_ITEMS)"
            msg = msg % str(target)
            self.error.add_error("l2-invalid-parameter_name", msg)
            self.error.raise_on_error()

        # Test if parameter exists
        # (older l2i files might not have all parameters)
        try:
            parameter = getattr(self, target)
        except AttributeError:
            return

        # Set values, uncertainty bias
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
        self.elev[ii] = self.altitude[ii] - retracker.range[ii]
        self.elev.uncertainty[ii] = retracker.uncertainty[ii]

    def set_metadata(self, auxdata_source_dict=None,
                     source_primary_filename=None,
                     l2_algorithm_id=None, l2_version_tag=None):
        if auxdata_source_dict is not None:
            self._auxdata_source_dict = auxdata_source_dict
        if source_primary_filename is not None:
            self._source_primary_filename = source_primary_filename
        if l2_algorithm_id is not None:
            self._l2_algorithm_id = l2_algorithm_id
        if l2_version_tag is not None:
            self._l2_version_tag = l2_version_tag

    def get_parameter_by_name(self, parameter_name):
        """ Method to retrieve a level-2 parameter """

        # Combine parameter and property catalogs
        catalog = self._PARAMETER_CATALOG
        catalog.update(self._PROPERTY_CATALOG)

        if "_uncertainty" in parameter_name:
            parameter_name = parameter_name.replace("_uncertainty", "")
            source = catalog[parameter_name]
            parameter = getattr(self, source)
            return parameter.uncertainty
        elif "_bias" in parameter_name:
            parameter_name = parameter_name.replace("_bias", "")
            source = catalog[parameter_name]
            parameter = getattr(self, source)
            return parameter.bias
        else:
            source = catalog[parameter_name]
            parameter = getattr(self, source)
            return parameter

    def get_attribute(self, attribute_name, *args):
        """ Return a string for a given attribute name. This method is
        required for the output data handler """

        try:
            attr_getter = getattr(self, "_get_attr_"+attribute_name)
            attribute = attr_getter(*args)
            return attribute
        except AttributeError:
            return "unkown"

    def _create_l2_data_items(self):
        for item in self._L2_DATA_ITEMS:
            setattr(self, item, L2ElevationArray(shape=(self.n_records)))

    def _check_if_valid_parameter(self, parameter_name):
        """ Performs a test if parameter name is a valid level-2 parameter
        name. Adds error if result negative and returns flag (valid: True,
        invalid: False) """
        if parameter_name not in self._L2_DATA_ITEMS:
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
            is_np_array = isinstance(value, (np.ndarray, np.array))
            is_correct_size = self._check_valid_size(value)
            if is_np_array and is_correct_size:
                return value
            else:
                return np.full(self.arrshape, np.nan)

    def _get_attr_pysiral_version(self, target):
        versions = {"filename": PYSIRAL_VERSION_FILENAME,
                    "default": PYSIRAL_VERSION}
        return versions[target]

    def _get_attr_mission_id(self, *args):
        # XXX: Deprecated
        return self.info.mission

    def _get_attr_source_mission_id(self, *args):
        mission_id = self.info.mission
        if args[0] == "uppercase":
            mission_id = mission_id.upper()
        return mission_id

    def _get_attr_source_mission_name(self, *args):
        mission_name = MISSION_NAME_DICT[self.info.mission]
        if args[0] == "uppercase":
            mission_name = mission_name.upper()
        return mission_name

    def _get_attr_source_mission_sensor(self, *args):
        mission_sensor = SENSOR_NAME_DICT[self.info.mission]
        if args[0] == "uppercase":
            mission_sensor = mission_sensor.upper()
        return mission_sensor

    def _get_attr_source_hemisphere(self, *args):
        return self.hemisphere

    def _get_attr_hemisphere(self, *args):
        # XXX: Deprecated
        return self.hemisphere

    def _get_attr_hemisphere_code(self, *args):
        hemisphere_code = self.hemisphere_code
        if args[0] == "uppercase":
            hemisphere_code = hemisphere_code.upper()
        return hemisphere_code

    def _get_attr_startdt(self, dtfmt):
        # XXX: Deprecated
        return self.info.start_time.strftime(dtfmt)

    def _get_attr_stopdt(self, dtfmt):
        # XXX: Deprecated
        return self.info.stop_time.strftime(dtfmt)

    def _get_attr_geospatial_lat_min(self, *args):
        return self._gett_attr_geospatial_str(np.nanmin(self.latitude))

    def _get_attr_geospatial_lat_max(self, *args):
        return self._gett_attr_geospatial_str(np.nanmax(self.latitude))

    def _get_attr_geospatial_lon_min(self, *args):
        return self._gett_attr_geospatial_str(np.nanmin(self.longitude))

    def _get_attr_geospatial_lon_max(self, *args):
        return self._gett_attr_geospatial_str(np.nanmax(self.longitude))

    def _gett_attr_geospatial_str(self, value):
        return "%.4f" % value

    def _get_attr_source_auxdata_sic(self, *args):
        value = self._auxdata_source_dict.get("sic", "unkown")
        if value == "unkown":
            value = self.info.source_auxdata_sic
        return value

    def _get_attr_source_auxdata_sitype(self, *args):
        value = self._auxdata_source_dict.get("sitype", "unkown")
        if value == "unkown":
            value = self.info.source_auxdata_sitype
        return value

    def _get_attr_source_auxdata_mss(self, *args):
        value = self._auxdata_source_dict.get("mss", "unkown")
        if value == "unkown":
            value = self.info.source_auxdata_mss
        return value

    def _get_attr_source_auxdata_snow(self, *args):
        value = self._auxdata_source_dict.get("snow", "unkown")
        if value == "unkown":
            value = self.info.source_auxdata_snow
        return value

    def _get_attr_source_sic(self, *args):
        # XXX: Deprecated
        return self._auxdata_source_dict.get("sic", "unkown")

    def _get_attr_source_sitype(self, *args):
        # XXX: Deprecated
        return self._auxdata_source_dict.get("sitype", "unkown")

    def _get_attr_source_mss(self, *args):
        # XXX: Deprecated
        return self._auxdata_source_dict.get("mss", "unkown")

    def _get_attr_source_snow(self, *args):
        # XXX: Deprecated
        return self._auxdata_source_dict.get("snow", "unkown")

    def _get_attr_source_primary(self, *args):
        return self._source_primary_filename

    def _get_attr_l2_algorithm_id(self, *args):
        return self._l2_algorithm_id

    def _get_attr_l2_version_tag(self, *args):
        return self._l2_version_tag

    def _get_attr_utcnow(self, *args):
        return self._creation_time.isoformat()

    def _get_attr_time_coverage_start(self, *args):
        datetime = self.period.start
        if re.match("%", args[0]):
            time_string = datetime.strftime(args[0])
        else:
            time_string = datetime.isoformat()
        return time_string

    def _get_attr_time_coverage_end(self, *args):
        datetime = self.period.stop
        if re.match("%", args[0]):
            time_string = datetime.strftime(args[0])
        else:
            time_string = datetime.isoformat()
        return time_string

    def _get_attr_time_coverage_duration(self, *args):
        return self.period.duration_isoformat

    def _get_attr_time_resolution(self, *args):
        tdelta = self.timestamp[-1]-self.timestamp[0]
        seconds = tdelta.total_seconds() + 1e-6 * tdelta.microseconds
        resolution = seconds/self.n_records
        return "%.2f seconds" % resolution

    def _get_attr_source_timeliness(self, *args):
        """ Return the timeliness of the l1b source data. Set default to
        NTC for backwark compability """
        try:
            timeliness = self.info.timeliness
        except AttributeError:
            timeliness = "NTC"
        if timeliness is None:
            timeliness = "NTC"
        if args[0] == "lowercase":
            timeliness = timeliness.lower()
        return timeliness

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
            (self.latitude[1], self.longitude[1]),
            (self.latitude[0], self.longitude[0])).meters

        if np.isclose(spacing, 0.0):
            spacing = great_circle(
                (self.latitude[-2], self.longitude[-2]),
                (self.latitude[-1], self.longitude[-1])).meters

        return spacing

    @property
    def dimdict(self):
        """ Returns dictionary with dimensions"""
        dimdict = OrderedDict([("n_records", self.n_records)])
        return dimdict

    @property
    def timestamp(self):
        try:
            time = self.track.time
        except AttributeError:
            time = self.track.timestamp
        return time

    @property
    def longitude(self):
        return self.track.longitude

    @property
    def latitude(self):
        return self.track.latitude

    @property
    def altitude(self):
        return self.track.altitude

    @property
    def surface_type_flag(self):
        return self.surface_type.flag

    @property
    def ssh(self):
        ssh = L2ElevationArray(shape=self._n_records)
        ssh.set_value(self.mss+self.ssa)
        ssh.set_uncertainty(self.ssa.uncertainty)
        return ssh


class Level2iMetadata(L1bMetaData):
    """ Container for Level-2 intermediate meta data (Essentially
    mimicks the L1bdata equivalent since the data location
    are idential. This also allows to directly use the l1b.info
    object directly) """

    def __init__(self):
        super(Level2iMetadata, self).__init__()


class Level2iTimeOrbit(L1bTimeOrbit):
    """ Container for Level-2 intermediate time orbit group (Essentially
    mimicks the L1bdata equivalent since the data location
    are idential. This also allows to directly use the l1b.time_orbit
    oject directly) """

    def __init__(self, **kwargs):
        """ Accepts `is_evenly_spaced` keyword (default: True).
        This based on the assumption that l2i data is evenly spaced (all
        along-track data points). For l2p data which excludes nan's
        this must be explicetely set to false """

        super(Level2iTimeOrbit, self).__init__(None, **kwargs)

    def from_l2i_stack(self, l2i_stack, index_list=None):
        """ Creates a TimeOrbit group object from l2i import. This is
        necessary when the Level2Data object shall be constructed from an
        l2i netcdf product.
        The index list can be used for subsetting (e.g. only use positions
        with valid freeboard, etc) """

        # Extract parameters from l2i stack
        try:
            time = l2i_stack["timestamp"]
        except KeyError:
            time = l2i_stack["time"]
        longitude = l2i_stack["longitude"]
        latitude = l2i_stack["latitude"]

        # Subset (if necessary)
        if index_list is not None:
            time = time[index_list]
            longitude = longitude[index_list]
            latitude = latitude[index_list]

        # Get dummy altitude
        dummy_altitude = np.full(longitude.shape, np.nan)

        # Set the timestamp
        self.time = time

        # Set the position
        self.set_position(longitude, latitude, dummy_altitude)

    def from_l2i_nc_import(self, l2i):
        """ Creates a TimeOrbit group object from l2i import. This is
        necessary when the Level2Data object shall be constructed from an
        l2i netcdf product """
        # Set the timestamp
        self.time = l2i.timestamp
        # Set the position
        dummy_altitude = np.full(l2i.longitude.shape, np.nan)
        self.set_position(self, l2i.longitude, l2i.latitude, dummy_altitude)


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


class Level2PContainer(DefaultLoggingClass):

    def __init__(self, period):
        super(Level2PContainer, self).__init__(self.__class__.__name__)
        self.error = ErrorStatus()
        self._period = period
        self._l2i_stack = []

    def append_l2i(self, l2i):
        self._l2i_stack.append(l2i)

    def get_merged_l2(self):
        """ Returns a Level2Data object with data from all l2i objects """

        # Merge the parameter
        data = self._get_merged_data(valid_mask="freeboard")

        # There are rare occasion, where no valid freeboard data is found for an entire day
        if len(data["longitude"]) == 0:
            return None

        # Set up a timeorbit group
        timeorbit = Level2iTimeOrbit()
        timeorbit.from_l2i_stack(data)

        # Use the first l2i object in stack to retrieve metadata
        l2i = self._l2i_stack[0]

        # Set up a metadata container
        metadata = Level2iMetadata()
        metadata.set_attribute("n_records", len(timeorbit.time))
        metadata.set_attribute("start_time", timeorbit.time[0])
        metadata.set_attribute("stop_time", timeorbit.time[-1])

        # XXX: Very ugly, but required due to a non-standard use of
        #      region_subset_set (originally idea to crop regions in
        #      Level-2 Processor)
        region_name = "north" if np.nanmean(data["latitude"]) > 0 else "south"
        metadata.subset_region_name = region_name

        # Retrieve the following constant attributes from the first
        # l2i object in the stack
        info = self.l2i_stack[0].info

        # Old notation (for backward compability)
        try:
            mission_id = info.mission_id
            # Transfer auxdata information
            metadata.source_auxdata_sic = l2i.info.source_sic
            metadata.source_auxdata_snow = l2i.info.source_snow
            metadata.source_auxdata_sitype = l2i.info.source_sitype
            metadata.source_auxdata_mss = l2i.info.source_mss

        # New (fall 2017) pysiral product notaion
        except AttributeError:
            mission_id = info.source_mission_id
            # Transfer auxdata information
            metadata.source_auxdata_sic = l2i.info.source_auxdata_sic
            metadata.source_auxdata_snow = l2i.info.source_auxdata_snow
            metadata.source_auxdata_sitype = l2i.info.source_auxdata_sitype
            metadata.source_auxdata_mss = l2i.info.source_auxdata_mss

        try:
            metadata.timeliness = l2i.info.source_timeliness
        except AttributeError:
            pass

        metadata.set_attribute("mission", mission_id)
        mission_sensor = SENSOR_NAME_DICT[mission_id]
        metadata.set_attribute("mission_sensor", mission_sensor)

        # Construct level-2 object
        l2 = Level2Data(metadata, timeorbit, period=self._period)

        #  Transfer the level-2 data items

        # 1. Get the list of parameters
        # (assumuning all l2i files share the same)
        parameter_list_all = l2i.parameter_list

        # 2. Exclude variables that end with `_uncertainty`
        parameter_list = [p for p in parameter_list_all
                          if not re.search("_uncertainty", p)]

        # 3. Remove parameters from the timeorbit group, surface type &
        # orbit id. This will be added to level 2 object by other means
        # or do not make sense (surface type for valid freeboard will
        # always be sea ice)
        for parameter_name in ["timestamp", "time", "longitude", "latitude",
                               "surface_type"]:
            try:
                parameter_list.remove(parameter_name)
            except ValueError:
                pass

        # 4. Set parameters
        for parameter_name in parameter_list:

            # Get the parameter
            value = data[parameter_name]

            # Test if uncertainty exists
            uncertainty_name = parameter_name+"_uncertainty"
            if uncertainty_name in parameter_list_all:
                uncertainty = data[uncertainty_name]
            else:
                uncertainty = np.full(value.shape, 0.0)

            # Add to l2 object
            l2.set_parameter(parameter_name, value, uncertainty=uncertainty)

        return l2

    def _get_merged_data(self, valid_mask=None):
        """ Returns a dict with merged data groups for all parameters
        in the l2i file (assumed to be identical for all files in the stack
        """
        parameter_list = self.l2i_stack[0].parameter_list
        data = self._get_empty_data_group(parameter_list)
        for l2i in self.l2i_stack:
            if valid_mask is not None:
                valid_mask_parameter = getattr(l2i, valid_mask)
                is_valid = np.where(np.isfinite(valid_mask_parameter))[0]
            else:
                is_valid = np.arange(l2i.n_records)
            for parameter in parameter_list:
                stack_data = getattr(l2i, parameter)
                stack_data = stack_data[is_valid]
                data[parameter] = np.append(data[parameter], stack_data)
        return data

    def _get_empty_data_group(self, parameter_list):
        data = {}
        for parameter_name in parameter_list:
            data[parameter_name] = np.array([], dtype=np.float32)
        return data

    @property
    def l2i_stack(self):
        return self._l2i_stack

    @property
    def n_l2i_objects(self):
        return len(self.l2i_stack)

    @property
    def period(self):
        return self._period


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
        self.parameter_list = []
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
            self.parameter_list.append(parameter_name)
            setattr(self, parameter_name, getattr(content, parameter_name))

        self._n_records = len(self.longitude)

        # Get mission id from filename
#        l2i_filename = filename_from_path(self.filename)
#        filenaming = PysiralOutputFilenaming()
#        filenaming.parse_filename(l2i_filename)
#        self.mission = filenaming.mission_id

        # Get timestamp (can be either time or timestamp in l2i files)
        if hasattr(self, "time"):
            time = self.time
            time_parameter_name = "time"
        else:
            time = self.timestamp
            time_parameter_name = "timestamp"
        self._time_parameter_name = time_parameter_name
        dt = num2date(time, self.time_def.units, self.time_def.calendar)
        setattr(self, "time", dt)
        self.timestamp = self.time

    def transfer_nan_mask(self, source, targets):
        source_parameter = getattr(self, source)
        nan_indices = np.where(np.isnan(source_parameter))
        for target in targets:
            parameter = getattr(self, target)
            parameter[nan_indices] = np.nan
            setattr(self, target, parameter)

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

    @property
    def mission(self):
        if not hasattr(self, "info"):
            return None
        try:
            return self.info.mission_id
        except AttributeError:
            pass
        try:
            return self.info.source_mission_id
        except AttributeError:
            pass
        return None

    @property
    def timeliness(self):
        if not hasattr(self, "info"):
            return None
        try:
            return self.info.source_timeliness
        except AttributeError:
            return "NTC"
