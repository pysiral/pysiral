# -*- coding: utf-8 -*-
"""
Created on Tue Jul 07 14:10:34 2015

@author: Stefan

L1bdata is a data container that unifies radar altimeter L1b orbit data
from different missions. It allows subsetting and merging of adjacent
orbit segments. L1bdata can be stored as a netCDF file, thus allowing
faster access to pre-processed subsets of RA orbit data for L2 processing.

The scope of the L1bdata container comprises:

---------
Metadata:
---------

- descriptors of RA source data
- period and geographical location
- processing history (subsetted, merged)
- software version

-------------
Waveform Data
-------------

- waveform echo power
  dimension: (n_records, n_range_bins)
- range for each range bin to the satellite in meters
  dimension: (n_records, n_range_bins)
- radar mode flag for each waveform:
    0: LRM
    1: SAR
    2: SIN
  (this is necessary for merging CryoSat-2 SAR and SIN adjacent orbit segments)
- summarizing flag from source data
    0: invalid
    1: valid
- optional: Additional named flags


----------------------
Time-Orbit Information
----------------------

- timestamp in UTC
- longitude, latitude (of satellite/nadir point)
- altitude (of satellite above WGS84 reference ellipsoid)

All parameter are of dimension (n_records).


-----------------
Range Corrections
-----------------

A list of range corrections (usually from RA source data files). The list
of correction is not predefined, but usally contains range corrections for:

- dry troposphere
- wet troposphere
- ionosphere
- inverse barometric / dynamic atmosphere
- ocean tide
- solid earth tide
- long period tide
- pole tide
- tidal loading

All parameter are of dimension (n_records) and of unit meter


----------
Classifier
----------

A list of optional named parameters that can be used for waveform
classification in the L2 processor. (e.g. stack parameter from the
CryoSat-2 l1b files)

All parameter are of dimension (n_records)


------------
Surface Type
------------



"""

from pysiral.surface_type import SurfaceType
from pysiral.output import NCDateNumDef
from pysiral.config import RadarModes

from netCDF4 import Dataset, num2date
from collections import OrderedDict
import numpy as np
import copy
import os


class Level1bData(object):
    """
    Unified L1b Data Class
    """
    data_groups = ["time_orbit", "correction", "classifier",
                   "waveform", "surface_type"]

    def __init__(self):
        self.info = L1bMetaData()
        self.waveform = L1bWaveforms(self.info)
        self.time_orbit = L1bTimeOrbit(self.info)
        self.correction = L1bRangeCorrections(self.info)
        self.classifier = L1bClassifiers(self.info)
        self.surface_type = SurfaceType()

    def append(self, l1b_annex):
        """ Appends another l1b object to this one """

        # Append data in each datagroup
        for data_group in self.data_groups:
            this_data_group = getattr(self, data_group)
            annex_data_group = getattr(l1b_annex, data_group)
            this_data_group.append(annex_data_group)

        # Update the statistics
        self.info.set_attribute("is_merged_orbit", True)
        self.info.set_attribute("n_records", len(self.time_orbit.timestamp))
        mission_data_source = ";".join([self.info.mission_data_source,
                                        l1b_annex.info.mission_data_source])
        self.info.set_attribute("mission_data_source", mission_data_source)
        self.update_l1b_metadata()

    def trim_to_subset(self, subset_list):
        """ Create a subset from an indix list """

        # Trim all datagroups
        for data_group in self.data_groups:
            content = getattr(self, data_group)
            content.set_subset(subset_list)

        # Update metadata
        self.info.set_attribute("is_orbit_subset", True)
        self.info.set_attribute("n_records", len(subset_list))
        self.update_l1b_metadata()

    def apply_range_correction(self, correction):
        """  Apply range correction """
        range_delta = self.correction.get_parameter_by_name(correction)
        if range_delta is None:
            # TODO: raise warning
            return
        self.waveform.add_range_delta(range_delta)

    def extract_subset(self, subset_list):
        """ Same as trim_to_subset, except returns a new l1bdata instance """
        if len(subset_list) > 0:
            l1b = copy.deepcopy(self)
            l1b.trim_to_subset(subset_list)
        else:
            return None
        return l1b

    def extract_region_of_interest(self, roi):
        """ Extracts data for a given region of interest definition """
        subset_list = roi.get_roi_list(self.time_orbit.longitude,
                                       self.time_orbit.latitude)
        if len(subset_list) > 0:
            l1b = copy.copy(self)
            l1b.trim_to_subset(subset_list)
        else:
            l1b = Level1bData()
        return l1b

    def update_l1b_metadata(self):
        self.update_data_limit_attributes()
        self.update_waveform_statistics()
        self.update_surface_type_statistics()
        self.update_region_name()

    def update_data_limit_attributes(self):
        """
        Set latitude/longitude and timestamp limits in the metadata container
        """

        info = self.info

        # time orbit group infos
        info.set_attribute("lat_min", np.nanmin(self.time_orbit.latitude))
        info.set_attribute("lat_max", np.nanmax(self.time_orbit.latitude))
        info.set_attribute("lon_min", np.nanmin(self.time_orbit.longitude))
        info.set_attribute("lon_max", np.nanmax(self.time_orbit.longitude))
        info.set_attribute("start_time", self.time_orbit.timestamp[0])
        info.set_attribute("stop_time", self.time_orbit.timestamp[-1])

    def update_waveform_statistics(self):
        """ Compute waveform metadata attributes """

        from pysiral.config import RadarModes

        # waveform property infos (lrm, sar, sarin)
        radar_modes = RadarModes()
        radar_mode = self.waveform.radar_mode

        # Check if radar mode is none
        # (e.g. if only header information is parsed at this stage)
        if radar_mode is None:
            return

        # Compute the percentage of each radar mode in the l1b object
        nrecs_fl = float(self.n_records)
        for flag in range(radar_modes.num):
            is_this_radar_mode = np.where(radar_mode == flag)[0]
            radar_mode_percent = 100.*float(len(is_this_radar_mode))/nrecs_fl
            attribute_name = "%s_mode_percent" % radar_modes.name(flag)
            self.info.set_attribute(attribute_name, radar_mode_percent)

    def update_surface_type_statistics(self):
        """ Re-calculate the open ocean percent """
        n_ocean_records = self.surface_type.get_by_name("ocean").num
        open_ocean_percent = 100.*float(n_ocean_records)/float(self.n_records)
        self.info.set_attribute("open_ocean_percent", open_ocean_percent)

    def update_region_name(self):
        """ Estimate the region (north/south/global) for metatdata class """

        lat_range = np.array([self.info.lat_min, self.info.lat_max])

        if np.amin(lat_range) > 0 and np.amax(lat_range) > 0:
            region_name = "north"
        elif np.amin(lat_range) < 0 and np.amax(lat_range) < 0:
            region_name = "south"
        else:
            region_name = "global"
        self.info.set_attribute("region_name", region_name)

    def reduce_waveform_bin_count(self, target_count, maxloc=0.4):
        """
        Reduce the bin count of waveform power and range arrays.
        (e.g. for merging CryoSat-2 SAR [256 bins] and SIN [1024 bins])

        Creates a subset and updates the l1b.waveform container

        Arguments
        ---------
        target_count (int)
            target number of waveform bins
            (needs to be smaller than full waveform bin count)

        Keywords
        --------
        maxloc (float, default=0.4)
            preferred location of the maximum of the waveform in the subset
        """
        # Extract original waveform
        orig_power, orig_range = self.waveform.power, self.waveform.range
        n_records, n_bins = orig_power.shape
        # Get the bin with the waveform maximum
        max_index = np.argmax(orig_power, axis=1)
        # Compute number of leading and trailing bins
        lead_bins = int(maxloc*target_count)
        trail_bins = target_count-lead_bins
        # Get the start/stop indeces for each waveform
        start, stop = max_index - lead_bins, max_index + trail_bins
        # Create new arrays
        rebin_shape = (n_records, target_count)
        power = np.ndarray(shape=rebin_shape, dtype=orig_power.dtype)
        range = np.ndarray(shape=rebin_shape, dtype=orig_range.dtype)
        # Validity check
        overflow = np.where(stop > n_bins)[0]
        if len(overflow) > 0:
            offset = n_bins - stop[overflow]
            stop[overflow] += offset
            start[overflow] += offset

        underflow = np.where(start < 0)[0]
        if len(underflow) > 0:
            offset = start[underflow]
            stop[underflow] -= offset
            start[underflow] -= offset
        # Extract the waveform with reduced bin count
        for i in np.arange(n_records):
            power[i, :] = orig_power[i, start[i]:stop[i]]
            range[i, :] = orig_range[i, start[i]:stop[i]]
        # Push to waveform container
        self.waveform.set_waveform_data(power, range, self.radar_modes)

    def get_parameter_by_name(self, data_group, parameter_name):
        """ API method to retrieve any parameter from any data group """
        try:
            data_group = getattr(self, data_group)
            return getattr(data_group, parameter_name)
        except:
            return None

    def set_parameter_by_name(self, data_group_name, parameter_name, value):
        """ API method to set any parameter in any data group """
        # Sanity check
        if len(value) != self.n_records:
            raise ValueError(
                    "value for %s.%s has wrong shape: %s" % (
                            data_group_name,
                            parameter_name,
                            str(value.shape)))
        try:
            # Get data group
            data_group = getattr(self, data_group_name)
            # Update data group
            setattr(data_group, parameter_name, value)
            # Update l1b container
            setattr(self, data_group_name, data_group)
        except:
            raise ValueError("Could not set value for %s.%s" % (
                   data_group_name, parameter_name))

    @property
    def n_records(self):
        return self.info.n_records

    @property
    def radar_modes(self):
        radar_modes = RadarModes()
        radar_mode_flag_list = np.unique(self.waveform.radar_mode)
        radar_mode_list = []
        for radar_mode_flag in radar_mode_flag_list:
            radar_mode_list.append(radar_modes.name(radar_mode_flag))
        return ";".join(radar_mode_list)


class L1bConstructor(Level1bData):
    """
    Class to be used to construct a L1b data object from any mission
    L1b data files
    """

    # TODO: should be coming from config file
    _SUPPORTED_MISSION_LIST = ["cryosat2", "envisat", "ers1", "ers2",
                               "sentinel3a", "icesat"]

    def __init__(self, config, header_only=False):
        super(L1bConstructor, self).__init__()
        self._config = config
        self._mission = None
        self._mission_options = None
        self._filename = None
        self._header_only = header_only
        self.error_status = False

    @property
    def mission(self):
        return self._mission

    @mission.setter
    def mission(self, value):
        if value in self._SUPPORTED_MISSION_LIST:
            self._mission = value
        else:
            # XXX: An ErrorHandler is needed here
            raise ValueError("Unsupported mission type")
        # Get mission default options
        self._mission_options = self._config.get_mission_defaults(value)

    @property
    def filename(self):
        return self._filename

    @filename.setter
    def filename(self, value):
        if os.path.isfile(value):
            self._filename = value
        else:
            # XXX: An ErrorHandler is needed here
            raise IOError("Not a valid path")

    def set_mission_options(self, **kwargs):
        self._mission_options = kwargs

    def construct(self):
        """ Parse the file and construct the L1bData object """
        adapter = get_l1b_adapter(self._mission)(self._config)
        adapter.filename = self.filename
        adapter.construct_l1b(self)

    def get_header_info(self):
        adapter = get_l1b_adapter(self._mission)(self._config)
        adapter.filename = self.filename
        adapter.construct_l1b(self, header_only=True)


class L1bdataNCFile(Level1bData):

    def __init__(self, filename):

        super(L1bdataNCFile, self).__init__()
        self.filename = filename
        self.nc = None
        self.time_def = NCDateNumDef()
        self.ncattrs_ignore_list = ['_NCProperties']

    def parse(self):
        """ populated the L1b data container from the l1bdata netcdf file """
        self.nc = Dataset(self.filename, "r")
        self._import_metadata()
        self._import_timeorbit()
        self._import_waveforms()
        self._import_corrections()
        self._import_surface_type()
        self._import_classifier()
        self.nc.close()

    def _import_metadata(self):
        """
        transfers l1b metadata attributes
        (stored as global attributes in l1bdata netCDF files)
        """
        for attribute_name in self.nc.ncattrs():
            if attribute_name in self.ncattrs_ignore_list:
                continue
            attribute_value = getattr(self.nc, attribute_name)
            # Convert timestamps back to datetime objects
            if attribute_name in ["start_time", "stop_time"]:
                attribute_value = num2date(
                    attribute_value, self.time_def.units,
                    calendar=self.time_def.calendar)
            # Convert flags (integers back to bool)
            if attribute_name in ["is_orbit_subset", "is_merged_orbit"]:
                attribute_value = bool(attribute_value)
            self.info.set_attribute(attribute_name, attribute_value)

    def _import_timeorbit(self):
        """
        transfers l1b timeorbit group
        (timeorbit datagroup in l1bdata netCDF files)
        """
        # Get the datagroup
        datagroup = self.nc.groups["time_orbit"]
        # Set satellite position data (measurement is nadir)
        self.time_orbit.set_position(
            datagroup.variables["longitude"][:],
            datagroup.variables["latitude"][:],
            datagroup.variables["altitude"][:])
        # Convert the timestamp to datetimes
        self.time_orbit.timestamp = num2date(
             datagroup.variables["timestamp"][:],
             self.time_def.units,
             calendar=self.time_def.calendar)

    def _import_waveforms(self):
        """
        transfers l1b waveform group
        (waveform datagroup in l1bdata netCDF files)
        """
        # Get the datagroup
        datagroup = self.nc.groups["waveform"]
        # Set waveform (measurement is nadir)
        self.waveform.set_waveform_data(
            datagroup.variables["power"][:],
            datagroup.variables["range"][:],
            datagroup.variables["radar_mode"][:])
        # Set the valid flag
        is_valid = datagroup.variables["is_valid"][:].astype(bool)
        self.waveform.set_valid_flag(is_valid)

    def _import_corrections(self):
        """
        transfers l1b corrections group
        (waveform corrections in l1bdata netCDF files)
        """
        # Get the datagroup
        datagroup = self.nc.groups["correction"]
        # Loop over parameters
        for key in datagroup.variables.keys():
            variable = np.array(datagroup.variables[key][:])
            self.correction.set_parameter(key, variable)

    def _import_surface_type(self):
        """
        transfers l1b surface_type group
        (waveform corrections in l1bdata netCDF files)
        """
        # Get the datagroup
        datagroup = self.nc.groups["surface_type"]
        self.surface_type.set_flag(datagroup.variables["flag"][:])

    def _import_classifier(self):
        """
        transfers l1b corrections group
        (waveform corrections in l1bdata netCDF files)
        """
        # Get the datagroup
        datagroup = self.nc.groups["classifier"]
        # Loop over parameters
        for key in datagroup.variables.keys():
            variable = np.array(datagroup.variables[key][:])
            self.classifier.add(variable, key)


class L1bMetaData(object):
    """
    Container for L1B Metadata information
    (see property attribute_list for a list of attributes)
    """

    _attribute_list = [
        "pysiral_version", "mission", "mission_data_version",
        "mission_sensor", "mission_data_source", "n_records", "orbit",
        "cycle", "sar_mode_percent", "lrm_mode_percent", "sin_mode_percent",
        "is_orbit_subset", "is_merged_orbit", "start_time", "stop_time",
        "region_name", "lat_min", "lat_max", "lon_min", "lon_max",
        "open_ocean_percent"]

    def __init__(self):
        # Init all fields
        for field in self.attribute_list:
            setattr(self, field, None)
        # Set some fields to False (instead of none)
        self.orbit = 999999
        self.is_orbit_subset = False
        self.is_merged_orbit = False
        self.n_records = -1

    def __repr__(self):
        output = "pysiral.L1bdata object:\n"
        for field in self._attribute_list:
            output += "%22s: %s" % (field, getattr(self, field))
            output += "\n"
        return output

    @property
    def attribute_list(self):
        return self._attribute_list

    @property
    def attdict(self):
        """ Return attributes as dictionary (e.g. for netCDF export) """
        attdict = {}
        for field in self.attribute_list:
            attdict[field] = getattr(self, field)
        return attdict

    @property
    def hemisphere(self):
        hemisphere = "global"
        if self.lat_min > 0. and self.lat_max > 0.:
            hemisphere = "north"
        if self.lat_min < 0. and self.lat_max < 0.0:
            hemisphere = "south"
        return hemisphere

    @property
    def year(self):
        return "%04g" % self.start_time.year

    @property
    def month(self):
        return "%02g" % self.start_time.month

    def set_attribute(self, tag, value):
        if tag not in self.attribute_list:
            raise ValueError("Unknown attribute: ", tag)
        setattr(self, tag, value)

    def check_n_records(self, n_records):
        # First time a data set is set: Store number of records as reference
        if self.n_records == -1:
            self.n_records = n_records
        else:  # n_records exists: verify consistenty
            if n_records == self.n_records:  # all good
                pass
            else:  # raise Erro
                raise ValueError("n_records mismatch, len must be: ",
                                 str(self.n_records))


class L1bTimeOrbit(object):

    """ Container for Time and Orbit Information of L1b Data """
    def __init__(self, info, is_evenly_spaced=True):
        self._info = info  # Pointer to metadata container
        self._timestamp = None
        self._longitude = None
        self._latitude = None
        self._altitude = None
        self._is_evenly_spaced = is_evenly_spaced

    @property
    def longitude(self):
        return self._longitude

    @property
    def latitude(self):
        return self._latitude

    @property
    def altitude(self):
        return self._altitude

    @property
    def timestamp(self):
        return self._timestamp

    @timestamp.setter
    def timestamp(self, value):
        if self._info is not None:
            self._info.check_n_records(len(value))
        self._timestamp = value

    @property
    def parameter_list(self):
        return ["timestamp", "longitude", "latitude", "altitude"]

    @property
    def dimdict(self):
        """ Returns dictionary with dimensions"""
        dimdict = OrderedDict([("n_records", len(self._timestamp))])
        return dimdict

    @property
    def is_evenly_spaced(self):
        return self._is_evenly_spaced

    def set_position(self, longitude, latitude, altitude):
        # Check dimensions
        if self._info is not None:
            self._info.check_n_records(len(longitude))
            self._info.check_n_records(len(latitude))
            self._info.check_n_records(len(altitude))
        # All fine => set values
        self._longitude = longitude
        self._latitude = latitude
        self._altitude = altitude

    def append(self, annex):
        for parameter in self.parameter_list:
            this_data = getattr(self, "_"+parameter)
            annex_data = getattr(annex, parameter)
            this_data = np.append(this_data, annex_data)
            setattr(self,  "_"+parameter, this_data)

    def set_subset(self, subset_list):
        for parameter in self.parameter_list:
            data = getattr(self, "_"+parameter)
            data = data[subset_list]
            setattr(self,  "_"+parameter, data)

    def __getstate__(self):
        return self.__dict__

    def __setstate__(self, d):
        self.__dict__.update(d)


class L1bRangeCorrections(object):
    """ Container for Range Correction Information """

    def __init__(self, info):
        self._info = info  # Pointer to Metadata object
        self._parameter_list = []

    def set_parameter(self, tag, value):
        self._info.check_n_records(len(value))
        setattr(self, tag, value)
        self._parameter_list.append(tag)

    @property
    def parameter_list(self):
        return self._parameter_list

    @property
    def n_records(self):
        parameter, name = self.get_parameter_by_index(0)
        return len(parameter)

    @property
    def dimdict(self):
        """ Returns dictionary with dimensions"""
        dimdict = OrderedDict([("n_records", self.n_records)])
        return dimdict

    def get_parameter_by_index(self, index):
        name = self._parameter_list[index]
        return getattr(self, name), name

    def get_parameter_by_name(self, name):
        try:
            return getattr(self, name)
        except:
            return None

    def append(self, annex):
        for parameter in self.parameter_list:
            this_data = getattr(self, parameter)
            annex_data = getattr(annex, parameter)
            this_data = np.append(this_data, annex_data)
            setattr(self, parameter, this_data)

    def set_subset(self, subset_list):
        for parameter in self.parameter_list:
            data = getattr(self, parameter)
            data = data[subset_list]
            setattr(self, parameter, data)


class L1bClassifiers(object):
    """ Containier for parameters that can be used as classifiers """

    def __init__(self, info):
        self._info = info  # Pointer to Metadate object
        # Make a pre-selection of different classifier types
        self._list = {
            "surface_type": [],
            "warning": [],
            "error": []}

    def add(self, value, name, classifier_type="surface_type"):
        """ Add a parameter for a given classifier type """
        setattr(self, name, np.array(value))
        self._list[classifier_type].append(name)

    @property
    def parameter_list(self):
        parameter_list = []
        for key in self._list.keys():
            parameter_list.extend(self._list[key])
        return parameter_list

    @property
    def n_records(self):
        parameter_list = self.parameter_list
        if len(parameter_list) == 0:
            return 0
        else:
            return len(getattr(self, parameter_list[0]))

    @property
    def dimdict(self):
        """ Returns dictionary with dimensions"""
        dimdict = OrderedDict([("n_records", self.n_records)])
        return dimdict

    def has_parameter(self, parameter_name):
        return parameter_name in self.parameter_list

    def get_parameter(self, parameter_name):
            return getattr(self, parameter_name)

    def append(self, annex):
        for parameter in self.parameter_list:
            this_data = getattr(self, parameter)
            annex_data = getattr(annex, parameter)
            this_data = np.append(this_data, annex_data)
            setattr(self, parameter, this_data)

    def set_subset(self, subset_list):
        for parameter in self.parameter_list:
            data = getattr(self, parameter)
            data = data[subset_list]
            setattr(self, parameter, data)


class L1bWaveforms(object):
    """ Container for Echo Power Waveforms """

    _valid_radar_modes = ["lrm", "sar", "sin"]
    _parameter_list = ["power", "range", "radar_mode", "is_valid"]
    _attribute_list = ["echo_power_unit"]

    def __init__(self, info):
        self._info = info  # Pointer to Metadate object
        # Attributes
        self.echo_power_unit = None
        self.radar_mode_def = RadarModes()
        # Parameter
        self._power = None
        self._range = None
        self._radar_mode = None
        self._is_valid = None

    @property
    def power(self):
        return np.copy(self._power)

    @property
    def range(self):
        return np.copy(self._range)

    @property
    def radar_mode(self):
        return self._radar_mode

    @property
    def is_valid(self):
        return self._is_valid

    @property
    def parameter_list(self):
        return self._parameter_list

    @property
    def n_range_bins(self):
        return self._get_wfm_shape(1)

    @property
    def n_records(self):
        return self._get_wfm_shape(0)

    @property
    def radar_modes(self):
        if self._radar_mode is None:
            return "none"
        flags = np.unique(self._radar_mode)
        return [self.radar_mode_def.get_name(flag) for flag in flags]

    @property
    def dimdict(self):
        """ Returns dictionary with dimensions"""
        shape = np.shape(self._power)
        dimdict = OrderedDict([("n_records", shape[0]), ("n_bins", shape[1])])
        return dimdict

    def set_waveform_data(self, power, range, radar_mode):
        # Validate input
        if power.shape != range.shape:
            raise ValueError("power and range must be of same shape",
                             power.shape, range.shape)
        if len(power.shape) != 2:
            raise ValueError("power and range arrays must be of dimension" +
                             " (n_records, n_bins)")

            # Validate number of records
        self._info.check_n_records(power.shape[0])

        # Assign values
        self._power = power
        self._range = range

        # Create radar mode arrays
        if type(radar_mode) is str and radar_mode in self._valid_radar_modes:
            mode_flag = self.radar_mode_def.get_flag(radar_mode)
            self._radar_mode = np.repeat(
                mode_flag, self.n_records).astype(np.byte)
        elif len(radar_mode) == self._info.n_records and \
                radar_mode.dtype == "int8":
                self._radar_mode = radar_mode
        else:
            raise ValueError("Invalid radar_mode: ", radar_mode)

        # Set valid flag (assumed to be valid for all waveforms)
        # Flag can be set separately using the set_valid_flag method
        if self._is_valid is None:
            self._is_valid = np.ones(shape=(self.n_records), dtype=bool)

    def set_valid_flag(self, valid_flag):
        # Validate number of records
        self._info.check_n_records(len(valid_flag))
        self._is_valid = valid_flag

    def append(self, annex):
        self._power = np.concatenate((self._power, annex.power), axis=0)
        self._range = np.concatenate((self._range, annex.range), axis=0)
        self._radar_mode = np.append(self._radar_mode, annex.radar_mode)
        self._is_valid = np.append(self._is_valid, annex.is_valid)

    def set_subset(self, subset_list):
        self._power = self._power[subset_list, :]
        self._range = self._range[subset_list, :]
        self._radar_mode = self._radar_mode[subset_list]
        self._is_valid = self._is_valid[subset_list]

    def add_range_delta(self, range_delta):
        # XXX: Should range delta needs to be reshaped?
        range_delta_reshaped = np.repeat(range_delta, self.n_range_bins)
        range_delta_reshaped = range_delta_reshaped.reshape(
            self.n_records, self.n_range_bins)
        self._range += range_delta_reshaped

    def _get_wfm_shape(self, index):
        shape = np.shape(self._power)
        return shape[index]


def get_l1b_adapter(mission):
    """ Select and returns the correct IO Adapter for the specified mission """

    from pysiral.io_adapter import (
        L1bAdapterCryoSat, L1bAdapterEnvisat, L1bAdapterERS1,
        L1bAdapterERS2, L1bAdapterSentinel3A, L1bAdapterICESat)

    if mission == "cryosat2":
        return L1bAdapterCryoSat
    elif mission == "envisat":
        return L1bAdapterEnvisat
    elif mission == "ers1":
        return L1bAdapterERS1
    elif mission == "ers2":
        return L1bAdapterERS2
    elif mission == "sentinel3a":
        return L1bAdapterSentinel3A
    elif mission == "icesat":
        return L1bAdapterICESat
    else:
        raise ValueError("Unknown mission id: %s" % mission)
