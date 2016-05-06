# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 13:57:56 2016

@author: Stefan
"""
from pysiral.config import options_from_dictionary
from pysiral.filter import idl_smooth

from pyproj import Proj
import numpy as np


class SnowBaseClass(object):

    def __init__(self):
        self._options = None
        self._local_repository = None
        self._subfolders = []

    def set_options(self, **opt_dict):
        self._options = options_from_dictionary(**opt_dict)

    def set_local_repository(self, path):
        self._local_repository = path

    def set_filenaming(self, filenaming):
        self._filenaming = filenaming

    def set_subfolders(self, subfolder_list):
        self._subfolders = subfolder_list

    def get_along_track_snow(self, l2):
        snow_depth, snow_dens, msg = self._get_along_track_snow(l2)
        return snow_depth, snow_dens, msg


class Warren99(SnowBaseClass):

    sd_coefs = np.array([
        [28.01, 0.1270, -1.1833, -0.1164, -0.0051, 0.0243, 7.6, -0.06, 0.07, 4.6],
        [30.28, 0.1056, -0.5908, -0.0263, -0.0049, 0.0044, 7.9, -0.06, 0.08, 5.5],
        [33.89, 0.5486, -0.1996, 0.0280, 0.0216, -0.0176, 9.4, -0.04, 0.10, 6.2],
        [36.80, 0.4046, -0.4005, 0.0256, 0.0024, -0.0641, 9.4, -0.09, 0.09, 6.1],
        [36.93, 0.0214, -1.1795, -0.1076, -0.0244, -0.0142, 10.6, -0.21, 0.09, 6.3],
        [36.59, 0.7021, -1.4819, -0.1195, -0.0009, -0.0603, 14.1, -0.16, 0.12, 8.1],
        [11.02, 0.3008, -1.2591, -0.0811, -0.0043, -0.0959, 9.5, 0.02, 0.10, 6.7],
        [4.64, 0.3100, -0.6350, -0.0655, 0.0059, -0.0005, 4.6, -0.01, 0.05, 3.3],
        [15.81, 0.2119, -1.0292, -0.0868, -0.0177, -0.0723, 7.8, -0.03, 0.06, 3.8],
        [22.66, 0.3594, -1.3483, -0.1063, 0.0051, -0.0577, 8.0, -0.08, 0.06, 4.0],
        [25.57, 0.1496, -1.4643, -0.1409, -0.0079, -0.0258, 7.9, -0.05, 0.07, 4.3],
        [26.67, -0.1876, -1.4229, -0.1413, -0.0316, -0.0029, 8.2, -0.06, 0.07, 4.8]])

    swe_coefs = np.array([
        [8.37, -0.0270, -0.3400, -0.0319, -0.0056, -0.0005, 2.5, -0.005, 0.024, 1.6],
        [9.43, 0.0058, -0.1309, 0.0017, -0.0021, -0.0072, 2.6, -0.007, 0.028, 1.8],
        [10.74, 0.1618, 0.0276, 0.0213, 0.0076, -0.0125, 3.1, 0.007, 0.032, 2.1],
        [11.67, 0.0841, -0.1328, 0.0081, -0.0003, -0.0301, 3.2, -0.013, 0.032, 2.1],
        [11.80, -0.0043, -0.4284, -0.0380, -0.0071, -0.0063, 3.5, -0.047, 0.033, 2.2],
        [12.48, 0.2084, -0.5739, -0.0468, -0.0023, -0.0253, 4.9, -0.030, 0.044, 2.9],
        [4.01, 0.0970, -0.4930, -0.0333, -0.0026, -0.0343, 3.5, 0.008, 0.037, 2.4],
        [1.08, 0.0712, -0.1450, -0.0155, 0.0014, -0.0000, 1.1, -0.001, 0.012, 0.8],
        [3.84, 0.0393, -0.2107, -0.0182, -0.0053, -0.0190, 2.0, -0.003, 0.016, 1.0],
        [6.24, 0.1158, -0.2803, -0.0215, 0.0015, -0.0176, 2.3, -0.005, 0.021, 1.4],
        [7.54, 0.0567, -0.3201, -0.0284, -0.0032, -0.0129, 2.4, -0.000, 0.023, 1.5],
        [8.00, -0.0540, -0.3650, -0.0362, -0.0112, -0.0035, 2.5, -0.003, 0.024, 1.5]])

    earth_radius = 6371000.8
    water_density = 1024.0

    def __init__(self):
        super(Warren99, self).__init__()
        self._p = Proj(proj="stere", lat_0=90, lon_0=-90, lat_ts=70)

    def _get_along_track_snow(self, l2):
        # Validate hemisphere
        if l2.hemisphere == "south":
            snow_depth = np.zeros(l2.n_records)
            snow_density = np.zeros(l2.n_records)
            msg = "ERROR - Warren99 not valid for southern hemisphere"
            return snow_depth, snow_density, msg
        # Get orginial warren values
        snow_depth, snow_density = self._get_warren99_fit(l2)
        # Filter invalid values
        valid_min, valid_max = self._options.valid_snow_depth_range
        invalid = np.logical_or(snow_depth < valid_min, snow_depth > valid_max)
        invalid_records = np.where(invalid)[0]
        snow_depth[invalid_records] = np.nan
        snow_density[invalid_records] = np.nan
        # Apply ice_type (myi_fraction correction)
        scaling = l2.sitype * self._options.fyi_correction_factor + 0.5
        snow_depth *= scaling
        # Smooth snow depth (if applicable)
        if self._options.smooth_snow_depth:
            filter_width = self._options.smooth_filter_width_m
            # Convert filter width to index
            filter_width /= l2.footprint_spacing
            # Round to odd number
            filter_width = np.floor(filter_width) // 2 * 2 + 1
            snow_depth = idl_smooth(snow_depth, filter_width)
        return snow_depth, snow_density, ""

    def _get_warren99_fit(self, l2):
        # get projection coordinates
        month = l2.track.timestamp[0].month
        l2x, l2y = self._p(l2.track.longitude, l2.track.latitude)
        # convert to degrees of arc
        l2x = l2x/(self.earth_radius * np.pi/180.0)
        l2y = l2y/(self.earth_radius * np.pi/180.0)
        # Get W99 snow depth
        snow_depth = self._get_snow_depth(month, l2x, l2y)
        # Get W99 snow density
        snow_density = self._get_snow_density(snow_depth, month, l2x, l2y)
        return snow_depth, snow_density

    def _get_sd_coefs(self, month):
        return self.sd_coefs[month-1, 0:6]

    def _get_swe_coefs(self, month):
        return self.swe_coefs[month-1, 0:6]

    def _get_snow_depth(self, month, l2x, l2y):
        sd = self._get_sd_coefs(month)
        snow_depth = sd[0] + sd[1]*l2x + sd[2]*l2y + \
            sd[3]*l2x*l2y + sd[4]*l2x*l2x + sd[5]*l2y*l2y
        snow_depth *= 0.01
        return snow_depth

    def _get_snow_density(self, snow_depth, month, l2x, l2y):
        # get snow water equivalent coefs
        swe = self._get_swe_coefs(month)
        snow_water_equivalent = swe[0] + swe[1]*l2x + swe[2]*l2y + \
            swe[3]*l2x*l2y + swe[4]*l2x*l2x + swe[5]*l2y*l2y
        snow_water_equivalent *= 0.01
        # Convert sd and swe to snow density
        snow_density = snow_water_equivalent/snow_depth*self.water_density
        return snow_density


class FixedSnowDepthDensity(SnowBaseClass):
    """ Always returns zero snow depth """

    def __init__(self):
        super(FixedSnowDepthDensity, self).__init__()

    def _get_along_track_snow(self, l2):
        snow_depth = np.ones(shape=(l2.n_records), dtype=np.float32)
        snow_depth *= self._options.fixed_snow_depth
        snow_density = np.ones(shape=(l2.n_records), dtype=np.float32)
        snow_density *= self._options.fixed_snow_density
        return snow_depth, snow_density, ""


def get_l2_snow_handler(name):
    return globals()[name]()
