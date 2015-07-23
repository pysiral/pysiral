# -*- coding: utf-8 -*-
"""
Created on Tue Jul 07 14:10:34 2015

@author: Stefan
"""

import os

from pysiral.helper import parse_datetime_str
from pysiral.cryosat2.l1bfile import CryoSatL1B
from pysiral.cryosat2.functions import (get_structarr_attr, tai2utc,
                                        get_tai_datetime_from_timestamp)


class L1bData(object):
    """
    Unified L1b Data Class
    """
    def __init__(self):

        self.info = L1bMetaData()
        self.waveform = None
        self.time_orbit = L1bTimeOrbit()
        self.correction = None
        self.classifier = None


class L1bConstructor(L1bData):
    """
    Class to be used to construct a L1b data object from any mission
    L1b data files
    """

    _SUPPORTED_MISSION_LIST = ["cryosat2"]

    def __init__(self):
        super(L1bConstructor, self).__init__()
        self._mission = None
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

    def construct(self):
        """ Parse the file and construct the L1bData object """
        adapter = get_l1b_adapter(self._mission)()
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

    def set_position(self, longitude, latitude):
        # XXX: This is developing stuff
        self._longitude = longitude
        self._latitude = latitude


    def __init__(self):


def get_l1b_adapter(mission):
    """ XXX: Early development state only """
    if mission == "cryosat2":
        return L1bAdapterCryoSat
