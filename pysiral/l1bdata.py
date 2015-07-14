# -*- coding: utf-8 -*-
"""
Created on Tue Jul 07 14:10:34 2015

@author: Stefan
"""

import os

from pysiral.cryosat2.l1bfile import CryoSatL1B
from pysiral.cryosat2.functions import get_structarr_attr


class L1bData(object):
    """
    Unified L1b Data Class
    """
    def __init__(self):

        self._waveform_group = None
        self.time_orbit = L1bTimeOrbit()
        self._geocorrection_group = None
        self._classifier_group = None


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


class L1bTimeOrbit(object):
    """ Container for Time and Orbit Information of L1b Data """
    def __init__(self):
        self._timestamp = None
        self._longitude = None
        self._latitude = None
        self._altitude = None
        self._orbit = None

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

    def set_position(self, longitude, latitude):
        # XXX: This is developing stuff
        self._longitude = longitude
        self._latitude = latitude


class L1bAdapterCryoSat(object):
    """ Converts a CryoSat2 L1b object into a L1bData object """
    def __init__(self):
        self.filename = None

    def construct_l1b(self, l1bdata):
        # Store the pointer to the L1bData object
        self.l1bdata = l1bdata
        # Read the CryoSat-2 L1b data file
        self._read_cryosat2l1b()
        # Transfer the time and orbit data
        self._transfer_timeorbit()
        # Transfer the waveform data
        self._transfer_waveform_collection()
        # Transfer the range corrections
        self._transfer_range_corrections()
        # Transfer any classifier data
        self._transfer_classifiers()

    def _read_cryosat2l1b(self):
        """ Read the L1b file and create a CryoSat-2 native L1b object """
        self.cs2l1b = CryoSatL1B()
        self.cs2l1b.filename = self.filename
        self.cs2l1b.parse()
        error_status = self.cs2l1b.get_status()
        if error_status:
            # TODO: Needs ErrorHandler
            raise IOError()
        self.cs2l1b.post_processing()

    def _transfer_timeorbit(self):
        # Transfer the orbit position
        longitude = get_structarr_attr(self.cs2l1b.time_orbit, "longitude")
        latitude = get_structarr_attr(self.cs2l1b.time_orbit, "latitude")
        self.l1bdata.time_orbit.set_position(longitude, latitude)
        # Transfer the timestamp

    def _transfer_waveform_collection(self):
        pass

    def _transfer_range_corrections(self):
        pass

    def _transfer_classifiers(self):
        pass


def get_l1b_adapter(mission):
    """ XXX: Early development state only """
    if mission == "cryosat2":
        return L1bAdapterCryoSat
