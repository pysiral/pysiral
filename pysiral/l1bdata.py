# -*- coding: utf-8 -*-
"""
Created on Tue Jul 07 14:10:34 2015

@author: Stefan
"""

from pysiral.io_adapter import L1bAdapterCryoSat

import numpy as np
import os


class Level1bData(object):
    """
    Unified L1b Data Class
    """
    def __init__(self):

        self.info = L1bMetaData()
        self.waveform = L1bWaveforms()
        self.time_orbit = L1bTimeOrbit()
        self.correction = L1bRangeCorrections()
        self.classifier = L1bClassifiers()

    def trim_to_subset(self, subset_list):
        """
        Create a subset from an indice list
        """
        data_groups = ["time_orbit", "correction", "classifier", "waveform"]
        for data_group in data_groups:
            content = getattr(self, data_group)
            content.set_subset(subset_list)
        # TODO: Updating Metadata?

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

    _SUPPORTED_MISSION_LIST = ["cryosat2"]

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
    def __init__(self):
        self.mission = None
        self.mission_data_version = None
        self.radar_mode = None
        self.orbit = None
        self.start_time = None
        self.stop_time = None

    def __repr__(self):
        field_list = ["mission", "mission_data_version", "radar_mode",
                      "orbit", "start_time", "stop_time"]
        output = "pysiral.L1bdata object:\n"
        for field in field_list:
            output += "%22s: %s" % (field, getattr(self, field))
            output += "\n"
        return output


class L1bTimeOrbit(object):
    """ Container for Time and Orbit Information of L1b Data """
    def __init__(self):
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
        self._timestamp = value

    @property
    def parameter_list(self):
        return ["timestamp", "longitude", "latitude", "altitude"]

    def set_position(self, longitude, latitude, altitude):
        # XXX: This is developing stuff
        self._longitude = longitude
        self._latitude = latitude
        self._altitude = altitude

    def set_subset(self, subset_list):
        for parameter in self.parameter_list:
            data = getattr(self, "_"+parameter)
            data = data[subset_list]


class L1bRangeCorrections(object):
    """ Container for Range Correction Information """

    def __init__(self):
        self._parameter_list = []

    def set_parameter(self, tag, value):
        setattr(self, tag, value)
        self._parameter_list.append(tag)

    @property
    def parameter_list(self):
        return self._parameter_list

    def get_parameter_by_index(self, index):
        name = self._parameter_list[index]
        return getattr(self, name), name

    def set_subset(self, subset_list):
        for parameter in self.parameter_list:
            data = getattr(self, parameter)
            data = data[subset_list]


class L1bClassifiers(object):
    """ Containier for parameters that can be used as classifiers """
    def __init__(self):
        # Make a pre-selection of different classifier types
        self._list = {
            "surface_type": [],
            "warning": [],
            "error": []}

    def add_parameter(self, name, value, classifier_type):
        """ Add a parameter for a given classifier type """
        setattr(self, name, np.array(value))
        self._list[classifier_type].append(name)

    @property
    def parameter_list(self):
        parameter_list = []
        for key in self._list.keys():
            parameter_list.extend(self._list[key])
        return parameter_list

    def set_subset(self, subset_list):
        for parameter in self.parameter_list:
            data = getattr(self, parameter)
            data = data[subset_list]


class L1bWaveforms(object):
    """ Container for Echo Power Waveforms """
    def __init__(self):
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

    def add_waveforms(self, power, range):
        self._power = power
        self._range = range

    def set_subset(self, subset_list):
        self._power = self._power[subset_list, :]
        self._range = self._range[subset_list, :]


def get_l1b_adapter(mission):
    """ XXX: Early development state only """
    if mission == "cryosat2":
        return L1bAdapterCryoSat
