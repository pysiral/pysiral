import numpy as np

import cPickle as pickle
import logging

# SICCI library required for loading pickled SICCI (grid or vector objects)
try:
    import SICCI
except:
    pass

from pysiral.iotools import ReadNC
from pysiral.maptools import GeoPcolorGrid
from pysiral.path import (filename_from_path, file_basename)


class SICCIGrid(object):

    def __init__(self, filename, logging=True):
        self._filename = filename
        self.metadata = SICCIGridMetadata()
        self._logging = logging
        self._parameter_list = []
        self._get_metadata()

    def calculate_pcolor_grid(self, **grid_projection):
        self._log("Calculate pcolor grid ...")
        self.pcolor = GeoPcolorGrid(self.longitude, self.latitude)
        self.pcolor.calc_from_proj(**grid_projection)
        self._log("... done")

    def _register_parameter(self, param):
        self._parameter_list.append(param)

    @property
    def period_label(self):
        return get_month_name(self.metadata.month)+" "+str(self.metadata.year)

    @property
    def parameter_list(self):
        return self._parameter_list

    def _log(self, message):
        if self._logging:
            logging.info(message)


class SICCIGridMetadata(object):

    def __init__(self):
        self.year = None
        self.month = None
        self.week = None
        self.source = None


class SICCIGridFMIPickled(SICCIGrid):

    def __init__(self, filename, **kwargs):
        super(SICCIGridFMIPickled, self).__init__(filename, **kwargs)
        self.metadata.source = "SICCI-1 pickled"

    def parse(self):
        self._log(
            "Parsing sicci pickled grid file: %s" %
            filename_from_path(self._filename))
        G = pickle.load(open(self._filename))

        # Extract vector fb values.
        self.sea_ice_freeboard = np.ma.array(G.AverageFreeboard(False)[0])
        self.sea_ice_freeboard *= 1.0e-3  # in mm!
        self.sea_ice_freeboard.mask = np.isnan(self.sea_ice_freeboard)
        self.latitude = np.array(G.lat)  # in degrees
        self.longitude = np.array(G.lon)  # in degrees

        self._register_parameter("longitude")
        self._register_parameter("latitude")
        self._register_parameter("sea_ice_freeboard")

    def _get_metadata(self):
        basename = file_basename(self._filename)
        self.metadata.year = np.int(basename[8:12])
        self.metadata.month = np.int(basename[12:14])
        self._log("Year: %g, Month: %2.2g" % (
            self.metadata.year, self.metadata.month))


class SICCIGridCS2AWI(SICCIGrid):
    """ Grid Parameter based on cs2awi netCDF files """

    def __init__(self, filename, **kwargs):
        super(SICCIGridCS2AWI, self).__init__(filename, **kwargs)
        self.metadata.source = "cs2awi netcdf"

    def parse(self):
        self._log(
            "Parsing cs2awi grid file: %s" %
            filename_from_path(self._filename))
        nc = ReadNC(self._filename)
        for parameter in nc.parameters:
            data = np.ma.array(getattr(nc, parameter))
            data.mask = np.isnan(data)
            setattr(self, parameter, data)
            self._register_parameter(parameter)
            self._log("- found parameter: %s dims:%s" % (
                parameter, str(np.shape(data))))

    def _get_metadata(self):
        basename = file_basename(self._filename)
        # TODO: do proper parsing
        strarr = basename.split("_")
        self.metadata.year = np.int(strarr[2][0:4])
        self.metadata.month = np.int(strarr[2][4:6])
        self._log("Year: %g, Month: %2.2g" % (
            self.metadata.year, self.metadata.month))


def get_month_name(month_int):
    month_str = ["January", "February", "March", "April", "May", "June",
                 "July", "August", "September", "October", "November",
                 "December"]
    return month_str[month_int-1]
