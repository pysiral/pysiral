# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 16:30:24 2015

@author: Stefan
"""


import contextlib
from pysiral import psrlcfg
from pysiral.errorhandler import ErrorStatus
from pysiral.iotools import ReadNC
from pysiral.core import DefaultLoggingClass
from pysiral.l1bdata import L1bMetaData, L1bTimeOrbit

import numpy as np
from datetime import datetime
from loguru import logger
from geopy.distance import great_circle
from collections import OrderedDict
import uuid
import re


class Level2Data(object):

    # TODO: Using class variables may break things during multi-processing
    _L2_DATA_ITEMS = ["range", "sla", "sla_raw", "dot", "elev", "afrb", "frb", "sit", "radar_mode"]

    _HEMISPHERE_CODES = {"north": "nh", "south": "sh"}

    # These are only the standard Level-2 parameters
    # NOTE: Auxiliary parameter are handled differently
    _PARAMETER_CATALOG = {
        "time": "time",
        "longitude": "longitude",
        "latitude": "latitude",
        "surface_type": "surface_type_flag",
        "radar_mode": "radar_mode",
        "flag": "flag",
        "elevation": "elev",
        "sea_level_anomaly": "sla",
        "sea_level_anomaly_raw": "sla_raw",
        "dynamic_ocean_topography": "dot",
        "radar_freeboard": "afrb",
        "sea_ice_freeboard": "frb",
        "sea_ice_thickness": "sit"}

    _PROPERTY_CATALOG = {"sea_surface_height": "ssh"}

    def __init__(self, metadata, time_orbit, period=None):

        # Copy necessary fields form l1b
        self.error = ErrorStatus()
        self._n_records = metadata.n_records
        self.info = metadata
        self.track = time_orbit
        self.period = period

        self._auto_auxvar_num = 0

        # A dictionary similar to the parameter catalog
        # To be filled during the set auxdata method
        self._auxiliary_catalog = {}

        # A dictionary containing information for auxiliary data with
        # more then one dimension
        self._multidim_auxiliary_catalog = {}

        # Metadata
        self._auxdata_source_dict = {}
        self._source_primary_filename = "unkown"
        self._l2_algorithm_id = "unkown"
        self._l2_version_tag = "unkown"
        self._doi = ""
        self._data_record_type = None

        # Define time of dataset creation as the time of object initialization
        # to avoid slightly different timestamps for repated calls of datetime.now()
        self._creation_time = datetime.now()

        # Other Class properties
        self._is_evenly_spaced = time_orbit.is_evenly_spaced

        # Algorithm error flag
        self.flag = np.zeros(self.n_records, dtype=np.uint16)
        self.radar_mode = np.zeros(self.n_records, dtype=np.uint16)

        # Create Level2 Data Groups
        self._create_l2_data_items()

    def set_parameter(self, target, value, uncertainty=None):
        """ Convienience method to safely add a parameter with optional
        uncertainty and/or bias to the level-2 data structure """

        # Sanity checks
        # TODO: This will need to be overhauled
        is_l2_default = self._check_if_valid_parameter(target)

        # Check if the full name has been passed
        if not is_l2_default and target in self.parameter_catalog.keys():
            target = self.parameter_catalog[target]
        else:
            # TODO: Need to figure something out for the auxvar id (not known if reinstated from l2i)
            par_name = self.auto_auxvar_id
            self.set_auxiliary_parameter(par_name, target, value, uncertainty)
            return

        # Next check: Needs to be of correct shape
        is_correct_size = self._check_valid_size(value)
        if not is_correct_size:
            msg = "Invalid parameter dimension: %s (See self._L2_DATA_ITEMS)"
            msg %= str(target)
            self.error.add_error("l2-invalid-parameter_name", msg)
            self.error.raise_on_error()

        # Test if parameter exists
        # (older l2i files might not have all parameters)
        try:
            parameter = getattr(self, target)
        except AttributeError:
            return

        # Set values, uncertainty bias
        if isinstance(parameter, L2DataArray):
            parameter.set_value(value)
            if uncertainty is not None:
                uncertainty_value = self._get_as_array(uncertainty)
                parameter.set_uncertainty(uncertainty_value)
        else:
            parameter = value
        setattr(self, target, parameter)

    def set_auxiliary_parameter(self, var_id, var_name, value, uncertainty=None):
        """ Adds an auxiliary parameter to the data object"""

        # Use L2Elevation Array
        # TODO: This is to cumbersome, replace by xarray at due time
        param = L2DataArray(shape=self.n_records)
        # Allow value to be None
        # NOTE: In this case an empty value will be generated
        if value is None:
            value = np.full(self.n_records, np.nan)
        try:
            param.set_value(value)
        except ValueError:
            logger.error(f"Could not set auxiliary parameter: {var_name}")
        if uncertainty is not None:
            param.set_uncertainty(uncertainty)
        setattr(self, var_id, param)

        # Register auxiliary parameter (this allows to find the parameter
        # by its long name
        self._auxiliary_catalog[var_name] = var_id

    def set_multidim_auxiliary_parameter(self, var_id, var_name, value, dim_dict):
        """
        Adds an auxiliary parameter to the data object different dimensions than
        the standard (l2.n_records) default data array.
        :param var_id: (str) The target id for the variable
        :param var_name: (str) The long name of the variable
        :param value: (np.ndarray) The auxiliary parameter
        :param dim_dict: (dictionary) The dimenstion dictionary
        :return:
        """
        setattr(self, var_id, value)
        # Register auxiliary parameter (this allows to find the parameter by its long name
        self._auxiliary_catalog[var_name] = var_id
        self._multidim_auxiliary_catalog[var_name] = dim_dict

    def get_multidim_auxdata_dimdict(self, var_name):
        return self._multidim_auxiliary_catalog.get(var_name, None)

    def set_data_record_type(self, data_record_type):
        self._data_record_type = data_record_type

    def update_retracked_range(self, retracker):
        # Update only for indices (surface type) supplied by retracker class
        # TODO: should get an overhaul
        ii = retracker.indices
        self.range[ii] = retracker.range[ii]
        self.range.uncertainty[ii] = retracker.uncertainty[ii]
        self.elev[ii] = self.altitude[ii] - retracker.range[ii]
        self.elev.uncertainty[ii] = retracker.uncertainty[ii]

        # Register potential auxiliary data
        for var_id, var_name, value, uncertainty in retracker.auxdata_output:

            # --- Check if output variable already exists ---

            # Create if new
            if var_name not in self.auxvar_names:
                self.set_auxiliary_parameter(var_id, var_name, value, uncertainty)

            # Transfer values for indices if already exists
            else:
                auxdata = getattr(self, var_id)
                auxdata[ii] = value[ii]
                if uncertainty is not None:
                    auxdata.uncertainty[ii] = uncertainty[ii]
                setattr(self, var_id, auxdata)

    # TODO: Metadata needs to become a CF compliant data class
    def set_metadata(self, auxdata_source_dict=None, source_primary_filename=None, l2_algorithm_id=None,
                     l2_version_tag=None):
        if auxdata_source_dict is not None:
            self._auxdata_source_dict = auxdata_source_dict
        if source_primary_filename is not None:
            self._source_primary_filename = source_primary_filename
        if l2_algorithm_id is not None:
            self._l2_algorithm_id = l2_algorithm_id
        if l2_version_tag is not None:
            self._l2_version_tag = l2_version_tag

    def set_doi(self, doi):
        self._doi = doi

    def get_parameter_by_name(self, parameter_name, raise_on_error=True):
        """ Method to retrieve a level-2 parameter """

        # Combine parameter and property catalogs
        catalog = self.parameter_catalog
        catalog.update(self.property_catalog)
        catalog.update(self._auxiliary_catalog)

        if "_uncertainty" in parameter_name:
            parameter_name = parameter_name.replace("_uncertainty", "")
            source = catalog[parameter_name]
            parameter = getattr(self, source)
            return parameter.uncertainty

        else:
            try:
                source = catalog[parameter_name]
                parameter = getattr(self, source)
            except KeyError:
                msg = f"Variable name `{parameter_name}` is not in the catalog of this l2 object"

                self.error.add_error("l2data-missing-variable", msg)
                if raise_on_error:
                    self.error.raise_on_error()
                parameter = np.full(self.n_records, np.nan)
            return parameter

    def get_attribute(self, attribute_name, *args):
        """ Return a string for a given attribute name. This method is
        required for the output data handler """

        try:
            attr_getter = getattr(self, f"_get_attr_{attribute_name}")
            return attr_getter(*args)
        except AttributeError:
            return "unkown"

    def _create_l2_data_items(self):
        for item in self._L2_DATA_ITEMS:
            setattr(self, item, L2DataArray(shape=self.n_records))

    def _check_if_valid_parameter(self, parameter_name):
        """ Performs a test if parameter name is a valid level-2 parameter
        name. Adds error if result negative and returns flag (valid: True,
        invalid: False) """
        return parameter_name in self._L2_DATA_ITEMS

    def _check_valid_size(self, array, **kwargs):
        """ Test if array has the correct size shape=(n_records). Adds error
        if not and returns flag (valid: True, invalid: False) """
        condition = array.ndim == 1 and len(array) == self._n_records
        if condition:
            return True
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
        is_np_array = isinstance(value, (np.ndarray, np.array))
        is_correct_size = self._check_valid_size(value)
        if is_np_array and is_correct_size:
            return value
        else:
            return np.full(self.arrshape, np.nan)

    # TODO: All this needs to go to the metadata class with standardized attributes
    @staticmethod
    def _get_attr_pysiral_version(*args):
        return psrlcfg.version

    def _get_attr_mission_id(self, *args):
        return self.info.mission

    def _get_attr_source_mission_id(self, *args):
        mission_id = self.info.mission
        if args[0] == "uppercase":
            mission_id = mission_id.upper()
        if args[0] == "select":
            for entry in args[1:]:
                mission_id_code, label = entry.split(":")
                if mission_id == mission_id_code:
                    return label
            return f"Error (mission id {mission_id} not in select statement)"
        return mission_id

    def _get_attr_source_mission_name(self, *args):
        mission_name = psrlcfg.platforms.get_name(self.info.mission)
        if args[0] == "uppercase":
            mission_name = mission_name.upper()
        return mission_name

    def _get_attr_source_mission_sensor(self, *args):
        mission_sensor = psrlcfg.platforms.get_sensor(self.info.mission)
        if args[0] == "uppercase":
            mission_sensor = mission_sensor.upper()
        return mission_sensor

    def _get_attr_source_mission_sensor_fn(self, *args):
        """ Same as source mission sensor, only a sanitized version for filenames """
        mission_sensor = psrlcfg.platforms.get_sensor(self.info.mission)
        for character in ["-"]:
            mission_sensor = mission_sensor.replace(character, "")
        if args[0] == "uppercase":
            mission_sensor = mission_sensor.upper()
        return mission_sensor

    def _get_attr_source_hemisphere(self, *args):
        if args[0] == "select":
            choices = {"north": args[1], "south": args[2]}
            return choices.get(self.hemisphere, "n/a")
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
        return self.info.start_time.strftime(dtfmt)

    def _get_attr_stopdt(self, dtfmt):
        return self.info.stop_time.strftime(dtfmt)

    def _get_attr_geospatial_lat_min(self, *args):
        return self._gett_attr_geospatial_str(np.nanmin(self.latitude))

    def _get_attr_geospatial_lat_max(self, *args):
        return self._gett_attr_geospatial_str(np.nanmax(self.latitude))

    def _get_attr_geospatial_lon_min(self, *args):
        return self._gett_attr_geospatial_str(np.nanmin(self.longitude))

    def _get_attr_geospatial_lon_max(self, *args):
        return self._gett_attr_geospatial_str(np.nanmax(self.longitude))

    @staticmethod
    def _gett_attr_geospatial_str(value):
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
        dt = self._creation_time
        return dt.strftime(args[0]) if re.match("%", args[0]) else dt.isoformat()

    def _get_attr_time_coverage_start(self, *args):
        # Cryo-TEMPO change from start of invocation timeperiod to start of L2 object coverage
        dt = self.info.start_time
        return dt.strftime(args[0]) if re.match("%", args[0]) else dt.isoformat()

    def _get_attr_time_coverage_end(self, *args):
        # Cryo-TEMPO change from end of invocation timeperiod to end of L2 object coverage
        dt = self.info.stop_time
        return dt.strftime(args[0]) if re.match("%", args[0]) else dt.isoformat()

    def _get_attr_period_coverage_start(self, *args):
        dt = self.period.tcs.dt
        return dt.strftime(args[0]) if re.match("%", args[0]) else dt.isoformat()

    def _get_attr_period_coverage_end(self, *args):
        dt = self.period.tce.dt
        if re.match("%", args[0]):
            time_string = dt.strftime(args[0])
        else:
            time_string = dt.isoformat()
        return time_string

    def _get_attr_time_coverage_duration(self, *args):
        return self.period.duration.isoformat

    def _get_attr_time_resolution(self, *args):
        tdelta = self.time[-1] - self.time[0]
        seconds = tdelta.total_seconds() + 1e-6 * tdelta.microseconds
        resolution = seconds / self.n_records
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

    def _get_attr_cycle(self, *args):
        """ Return the cycle number of the l1b source data. Set default to
        -1 """
        try:
            cycle = str(self.info.cycle)
        except AttributeError:
            cycle = "-1"
        if cycle is None:
            cycle = "-1"
        return cycle

    def _get_attr_fncycle(self, *args):
        """ Return the cycle number of the l1b source data. Set default to
        -1 """
        try:
            fncycle = '{:0>2}'.format(str(self.info.cycle))
        except AttributeError:
            fncycle = "-1"
        if fncycle is None:
            fncycle = "-1"
        return fncycle

    def _get_attr_orbit(self, *args):
        """ Return the orbit number of the l1b source data. Set default to
        -1 """
        try:
            orbit = str(self.info.orbit)
        except AttributeError:
            orbit = "-1"
        if orbit is None:
            orbit = "-1"
        return orbit

    def _get_attr_rel_orbit(self, *args):
        """ Return the orbit number of the l1b source data. Set default to
        -1 """
        try:
            rel_orbit = str(self.info.rel_orbit)
        except AttributeError:
            rel_orbit = "-1"
        if rel_orbit is None:
            rel_orbit = "-1"
        return rel_orbit

    def _get_attr_fnrel_orbit(self, *args):
        """ Return the orbit number of the l1b source data. Set default to
        -1 """
        try:
            fnrel_orbit = '{:0>5}'.format(str(self.info.rel_orbit))
        except AttributeError:
            fnrel_orbit = "000-1"
        if fnrel_orbit is None:
            fnrel_orbit = "000-1"
        return fnrel_orbit

    @staticmethod
    def cycle_to_subcycle(cy, tr):
        NB_TRACKS = 10688  # cycle of 368.24 days
        NB_SUB_TRACKS = 840  # sub-cycle of 28.94 days

        absolute_tr = (cy - 1) * NB_TRACKS + tr
        sub_cy = ((absolute_tr - 1) // NB_SUB_TRACKS) + 1
        sub_tr = ((absolute_tr - 1) % NB_SUB_TRACKS) + 1
        return sub_cy, sub_tr

    # Add CNES sub-cycle. Need to check what to do after orbit change in Jul 2020
    def _get_attr_cnes_subcycle(self, *args):
        """ Return the orbit number of the l1b source data. Set default to
        -1 """
        try:
            cycle = int(self.info.cycle)
            rel_orbit = int(self.info.rel_orbit)
            cnes_subcycle, cnes_track = self.cycle_to_subcycle(cycle, rel_orbit)
            cnes_subcycle = str(cnes_subcycle)
        except AttributeError:
            cnes_subcycle = "-1"
        if cnes_subcycle is None:
            cnes_subcycle = "-1"
        return cnes_subcycle

    def _get_attr_cnes_track(self, *args):
        """ Return the orbit number of the l1b source data. Set default to
        -1 """
        try:
            cycle = int(self.info.cycle)
            rel_orbit = int(self.info.rel_orbit)
            cnes_subcycle, cnes_track = self.cycle_to_subcycle(cycle, rel_orbit)
            cnes_track = str(cnes_track)
        except AttributeError:
            cnes_track = "-1"
        if cnes_track is None:
            cnes_track = "-1"
        return cnes_track

    @staticmethod
    def _get_attr_uuid(*args):
        """ Provide an uuid code (for tracking id's) """
        return str(uuid.uuid4())

    def _get_attr_doi(self, *args):
        return self._doi

    @property
    def parameter_catalog(self):
        return dict(self._PARAMETER_CATALOG)

    @property
    def full_variable_catalog(self):
        full_variable_catalog = self.parameter_catalog
        full_variable_catalog.update(self.auxiliary_catalog)
        return full_variable_catalog

    @property
    def property_catalog(self):
        return dict(self._PROPERTY_CATALOG)

    @property
    def auxvar_names(self):
        return sorted(self._auxiliary_catalog.keys())

    @property
    def auto_auxvar_id(self):
        name = "auxvar%02g" % self._auto_auxvar_num
        self._auto_auxvar_num += 1
        return name

    @property
    def auxiliary_catalog(self):
        return dict(self._auxiliary_catalog)

    @property
    def arrshape(self):
        return self.n_records

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
        if self.n_records < 2:
            return np.nan

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
        return OrderedDict([("time", self.n_records)])

    @property
    def time(self):
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
        ssh = L2DataArray(shape=self._n_records)
        ssh.set_value(self.mss + self.sla)
        ssh.set_uncertainty(self.sla.uncertainty)
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

        # Properties
        self.time = None

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


# TODO: Is this going to stay, or better replaced by xarray variables?
class L2DataArray(np.ndarray):
    """
    Recipe from:
    http://docs.scipy.org/doc/numpy/user/basics.subclassing.html
    XXX: not yet full slicing capability! -> __getitem__ trouble
         always use cls[list] and cls.uncertainty[list]
         cls[list].uncertainty will fail
    """

    def __new__(cls, shape, dtype=float, buffer=None, offset=0, strides=None, order=None, info=None):
        obj = np.ndarray.__new__(cls, shape, dtype, buffer, offset, strides, order) * np.nan
        obj.uncertainty = np.zeros(shape=shape, dtype=float)
        obj.source_class = ""
        obj.source_files = ""
        obj.long_name = ""
        obj.unit = ""
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.uncertainty = getattr(obj, 'uncertainty', None)
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
        self[:] = value[:]

    def set_uncertainty(self, uncertainty):
        self.uncertainty = uncertainty

    def set_nan_indices(self, indices):
        value = self[:]
        value[indices] = np.nan
        self.set_value(value)
        self.uncertainty[indices] = np.nan


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
        data = self._get_merged_data(valid_mask="sea_ice_freeboard")

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

        # Old notation (for backward compatibility)
        # TODO: This will soon be obsolete
        mission_id = None
        if hasattr(info, "mission_id"):
            mission_id = info.mission_id
            metadata.source_auxdata_sic = l2i.info.source_sic
            metadata.source_auxdata_snow = l2i.info.source_snow
            metadata.source_auxdata_sitype = l2i.info.source_sitype
            metadata.source_auxdata_mss = l2i.info.source_mss

        # New (fall 2017) pysiral product notation
        if hasattr(info, "source_mission_id"):
            mission_id = info.source_mission_id
            # Transfer auxdata information
            metadata.source_auxdata_sic = l2i.info.source_auxdata_sic
            metadata.source_auxdata_snow = l2i.info.source_auxdata_snow
            metadata.source_auxdata_sitype = l2i.info.source_auxdata_sitype
            metadata.source_auxdata_mss = l2i.info.source_auxdata_mss

        # Conversion of l2i to CF/ACDD conventions (Fall 2021)
        if hasattr(info, "platform"):
            mission_id = psrlcfg.platforms.get_platform_id(info.platform)

        if mission_id is None:
            self.error.add_error("unknown-platform", "Cannot determine platform name from source l2i stack")
            self.error.raise_on_error()

        if hasattr(l2i.info, "source_timeliness"):
            metadata.timeliness = l2i.info.source_timeliness

        if hasattr(l2i.info, "data_record_type"):
            metadata.timeliness = l2i.info.data_record_type

        metadata.set_attribute("mission", mission_id)
        mission_sensor = psrlcfg.platforms.get_sensor(mission_id)
        metadata.set_attribute("mission_sensor", mission_sensor)

        # Construct level-2 object
        l2 = Level2Data(metadata, timeorbit, period=self._period)

        #  Transfer the level-2 data items

        # 1. Get the list of parameters
        # (assuming all l2i files share the same)
        parameter_list_all = l2i.parameter_list

        # 2. Exclude variables that end with `_uncertainty`
        parameter_list = [p for p in parameter_list_all if not re.search("_uncertainty", p)]

        # 3. Remove parameters from the timeorbit group, surface type &
        # orbit id. This will be added to level 2 object by other means
        # or do not make sense (surface type for valid freeboard will
        # always be sea ice)
        for parameter_name in ["timestamp", "time", "longitude", "latitude", "surface_type"]:
            with contextlib.suppress(ValueError):
                parameter_list.remove(parameter_name)
        # 4. Set parameters
        for parameter_name in parameter_list:

            # Get the parameter
            value = data[parameter_name]

            # Test if uncertainty exists
            uncertainty_name = parameter_name + "_uncertainty"
            if uncertainty_name in parameter_list_all:
                uncertainty = data[uncertainty_name]
            else:
                uncertainty = None

            # Add to l2 object
            l2.set_parameter(parameter_name, value, uncertainty=uncertainty)

        return l2

    def _get_merged_data(self, valid_mask: str = None) -> dict:
        """
        Returns a dict with merged data groups for all parameters
        in the l2i file (assumed to be identical for all files in the stack)

        :param valid_mask: The name of the parameter that defines the
            mask of valid l2i data points

        :return: Dictionary with all mergered l2i parameters
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

                # TODO: This needs better handling
                # NOTE: Some variables are dimensions and not data variables.
                #       These should not be concatenated
                if stack_data.size == l2i.n_records:
                    stack_data = stack_data[is_valid]
                    data[parameter] = np.append(data[parameter], stack_data)
                else:
                    data[parameter] = stack_data
        return data

    @staticmethod
    def _get_empty_data_group(parameter_list):
        return {parameter_name: np.array([], dtype=np.float32) for parameter_name in parameter_list}

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
        from cftime import num2pydate

        content = ReadNC(self.filename)

        for attribute_name in content.attributes:
            self.attribute_list.append(attribute_name)
            self.info.set_attribute(attribute_name,
                                    getattr(content, attribute_name))

        for parameter_name in content.parameters:
            self.parameter_list.append(parameter_name)
            setattr(self, parameter_name, getattr(content, parameter_name))

        self._n_records = len(self.longitude)

        # Get timestamp (can be either time or timestamp in l2i files)
        time_parameter_name = "time" if hasattr(self, "time") else "timestamp"
        time = self.time
        self._time_parameter_name = time_parameter_name
        dt = num2pydate(time, content.time_def.units, content.time_def.calendar)
        setattr(self, "time", dt)
        self.time = self.time

    def transfer_nan_mask(self, source, targets):
        source_parameter = getattr(self, source)
        nan_indices = np.where(np.isnan(source_parameter))
        for target in targets:
            parameter = getattr(self, target)
            parameter[nan_indices] = np.nan
            setattr(self, target, parameter)

    def mask_variables(self, indices, targets, value=np.nan):
        for target in targets:
            parameter = getattr(self, target)
            with contextlib.suppress(ValueError):
                parameter[indices] = value
            setattr(self, target, parameter)

    # TODO: Does this need to be here?
    def project(self, griddef):
        from pyproj import Proj
        p = Proj(**griddef.projection)
        self.projx, self.projy = p(self.longitude, self.latitude)
        # Convert projection coordinates to grid indices
        extent = griddef.extent
        self.xi = np.floor((self.projx + extent.xsize / 2.0) / extent.dx)
        self.yj = np.floor((self.projy + extent.ysize / 2.0) / extent.dy)

    @property
    def n_records(self):
        return self._n_records

    @property
    def mission(self):

        if not hasattr(self, "info"):
            return None

        if hasattr(self.info, "mission_id"):
            return self.info.mission_id

        if hasattr(self.info, "source_mission_id"):
            return self.info.source_mission_id

        if hasattr(self.info, "platform"):
            platform_name = self.info.platform
            return psrlcfg.platforms.get_platform_id(platform_name)

        return None

    @property
    def timeliness(self):

        if not hasattr(self, "info"):
            return "unknown"

        if hasattr(self.info, "source_mission_id"):
            return self.info.source_mission_id

        if hasattr(self.info, "data_record_type"):
            return self.info.data_record_type

        return "unknown"
