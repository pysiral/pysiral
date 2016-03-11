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

from pysiral.io_adapter import (L1bAdapterCryoSat, L1bAdapterEnvisat)
from pysiral.surface_type import SurfaceType

from collections import OrderedDict
import numpy as np
import os


class Level1bData(object):
    """
    Unified L1b Data Class
    """
    def __init__(self):
        self.info = L1bMetaData()
        self.waveform = L1bWaveforms(self.info)
        self.time_orbit = L1bTimeOrbit(self.info)
        self.correction = L1bRangeCorrections(self.info)
        self.classifier = L1bClassifiers(self.info)
        self.surface_type = SurfaceType()

    def trim_to_subset(self, subset_list):
        """
        Create a subset from an indice list
        """
        data_groups = ["time_orbit", "correction", "classifier", "waveform"]
        for data_group in data_groups:
            content = getattr(self, data_group)
            content.set_subset(subset_list)
        # TODO: Updating Metadata?

    def apply_range_correction(self, correction):
        """
        Apply range correction
        """
        range_delta = self.correction.get_parameter_by_name(correction)
        if range_delta is None:
            # TODO: raise warning
            return
        self.waveform.add_range_delta(range_delta)

    @property
    def n_records(self):
        try:
            n_records = len(self.time_orbit.timestamp)
        except:
            n_records = 0
        return n_records


class L1bConstructor(Level1bData):
    """
    Class to be used to construct a L1b data object from any mission
    L1b data files
    """

    _SUPPORTED_MISSION_LIST = ["cryosat2", "envisat"]

    def __init__(self, config):
        super(L1bConstructor, self).__init__()
        self._config = config
        self._mission = None
        self._mission_options = None
        self._filename = None

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


class L1bMetaData(object):
    """ Container for L1B Metadata information """

    field_list = ["mission", "mission_data_version", "radar_mode",
                  "orbit", "start_time", "stop_time"]

    def __init__(self):
        self.mission = None
        self.mission_data_version = None
        self.radar_mode = None
        self.orbit = None
        self.start_time = None
        self.stop_time = None

    def __repr__(self):

        output = "pysiral.L1bdata object:\n"
        for field in self.field_list:
            output += "%22s: %s" % (field, getattr(self, field))
            output += "\n"
        return output

    @property
    def attdict(self):
        """ Return attributes as dictionary (e.g. for netCDF export) """
        attdict = {}
        for field in self.field_list:
            attdict[field] = getattr(self, field)
        return attdict


class L1bTimeOrbit(object):
    """ Container for Time and Orbit Information of L1b Data """
    def __init__(self, info):
        self._info = info  # Pointer to metadata container
        self._timestamp = None
        self._longitude = None
        self._latitude = None
        self._altitude = None

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

    def set_position(self, longitude, latitude, altitude):
        # XXX: This is developing stuff
        self._longitude = longitude
        self._latitude = latitude
        self._altitude = altitude

    def set_subset(self, subset_list):
        for parameter in self.parameter_list:
            data = getattr(self, "_"+parameter)
            data = data[subset_list]
            setattr(self,  "_"+parameter, data)


class L1bRangeCorrections(object):
    """ Container for Range Correction Information """

    def __init__(self, info):
        self._info = info  # Pointer to Metadate object
        self._parameter_list = []

    def set_parameter(self, tag, value):
        setattr(self, tag, value)
        self._parameter_list.append(tag)

    @property
    def parameter_list(self):
        return self._parameter_list

    @property
    def dimdict(self):
        """ Returns dictionary with dimensions"""
        dimdict = OrderedDict([("n_records",
                                len(self.get_parameter_by_index(0)))])
        return dimdict

    def get_parameter_by_index(self, index):
        name = self._parameter_list[index]
        return getattr(self, name), name

    def get_parameter_by_name(self, name):
        try:
            return getattr(self, name)
        except:
            return None

    def set_subset(self, subset_list):
        for parameter in self.parameter_list:
            data = getattr(self, parameter)
            data = data[subset_list]


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

    def set_subset(self, subset_list):
        for parameter in self.parameter_list:
            data = getattr(self, parameter)
            data = data[subset_list]


class L1bWaveforms(object):
    """ Container for Echo Power Waveforms """

    def __init__(self, info):
        self._info = info  # Pointer to Metadate object
        self._power = None
        self._range = None

    @property
    def power(self):
        return self._power

    @property
    def range(self):
        return self._range

    @property
    def parameter_list(self):
        return ["power", "range"]

    @property
    def n_range_bins(self):
        return self._get_wfm_shape(1)

    @property
    def n_records(self):
        return self._get_wfm_shape(0)

    @property
    def dimdict(self):
        """ Returns dictionary with dimensions"""
        shape = np.shape(self._power)
        dimdict = OrderedDict([("n_records", shape[0]), ("n_bins", shape[1])])
        return dimdict

    def add_waveforms(self, power, range):
        self._power = power
        self._range = range

    def set_subset(self, subset_list):
        self._power = self._power[subset_list, :]
        self._range = self._range[subset_list, :]

    def add_range_delta(self, range_delta):
        range_delta_reshaped = np.repeat(range_delta, self.n_range_bins)
        range_delta_reshaped = range_delta_reshaped.reshape(
            self.n_records, self.n_range_bins)
        self._range += range_delta_reshaped

    def _get_wfm_shape(self, index):
        shape = np.shape(self._power)
        return shape[index]


def get_l1b_adapter(mission):
    """ XXX: Early development state only """
    if mission == "cryosat2":
        return L1bAdapterCryoSat
    if mission == "envisat":
        return L1bAdapterEnvisat
