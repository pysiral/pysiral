# -*- coding: utf-8 -*-
"""
Created on Sat Aug 01 17:03:19 2015

@author: Stefan
"""
from pysiral.io_tools import ReadNC

from scipy.interpolate import interp1d
from treedict import TreeDict
import scipy.ndimage as ndimage
import numpy as np


class BaseMSS(object):

    def __init__(self):
        pass

    def set_filename(self, filename):
        self._filename = filename

    def set_roi(self, roi):
        self._roi = roi

    def parse(self):
        self._parse()

    def roi_latitude_range(self):
        if hasattr(self, "_roi"):
            return self._roi.get_latitude_range()
        else:
            return [-90.0, 90.0]


class DTU1MinGrid(BaseMSS):
    """
    Parsing Routine for DTU 1 minute global mean sea surface height files
    """
    def __init__(self):
        super(DTU1MinGrid, self).__init__()

    def _parse(self):
        dtu_grid = ReadNC(self._filename)
        # Cut to ROI regions (latitude only)
        # -> no need for world mss
        latitude_range = self.roi_latitude_range()
        latitude_indices = np.where(
            np.logical_and(dtu_grid.lat >= latitude_range[0],
                           dtu_grid.lat <= latitude_range[1]))[0]
        self.elevation = dtu_grid.mss[latitude_indices, :]
        self.longitude = dtu_grid.lon
        self.latitude = dtu_grid.lat[latitude_indices]
        # Convert elevations to WGS84
        delta_h1 = egm2top_delta_h(self.latitude)
        delta_h = egm2wgs_delta_h(self.latitude)
        for i in np.arange(len(self.latitude)):
            self.elevation[i, :] += (delta_h[i]-delta_h1[i])

    def get_track(self, longitude, latitude):
        # Use fast image interpolation (since DTU is on regular grid)
        # Longitudes must be 0 -> 360
        negative_lons = np.where(longitude < 0)[0]
        longitude[negative_lons] = longitude[negative_lons] + 360.
        # Calculate image coordinates of mss grid "image"
        mss_lon_min = self.longitude[0]
        mss_lon_step = self.longitude[1] - self.longitude[0]
        mss_lat_min = self.latitude[0]
        mss_lat_step = self.latitude[1] - self.latitude[0]
        ix = (longitude - mss_lon_min)/mss_lon_step
        iy = (latitude - mss_lat_min)/mss_lat_step
        # Extract and return the elevation along the track
        mss_track_elevation = ndimage.map_coordinates(self.elevation, [iy, ix])
        return mss_track_elevation


class SSAInterpolator(object):
    """
    Parent class for sea surface anomaly retrieval and interpolation
    """
    def __init__(self):
        pass

    def set_options(self, **opt_dict):
        self._options = TreeDict.fromdict(opt_dict, expand_nested=True)

    def interpolate(self, l2):
        self._value = np.ndarray(shape=(l2.n_records))
        self._uncertainty = np.ndarray(shape=(l2.n_records))
        self._interpolate(l2)

    @property
    def value(self):
        return self._value

    @property
    def uncertainty(self):
        return self._uncertainty


class SSASmoothedLinear(SSAInterpolator):
    """ Default CS2AWI Method """

    def __init__(self):
        super(SSASmoothedLinear, self).__init__()

    def _interpolate(self, l2):
        self._linear_interpolation_between_leads(l2)
        self._calculate_uncertainty()
        self._smooth_interpolated_ssa()

    def _linear_interpolation_between_leads(self, l2):
        import matplotlib.pyplot as plt
        # Use ocean and lead elevations
        ssh_tiepoints = l2.surface_type.lead.indices
        if self._options.use_ocean_wfm:
            ssh_tiepoints.append(l2.surface_type.ocean.indices)
            ssh_tiepoints = np.sort(ssh_tiepoints)
        # Interpolate
        ssa_raw = np.ndarray(shape=(l2.n_records))*np.nan
        ssa_raw[ssh_tiepoints] = l2.elev[ssh_tiepoints]
        # Calculate filter width
        filter_width = self._options.smooth_filter_width_m/300.0
        filter_width = np.floor(filter_width) // 2 * 2 + 1
        print filter_width

        x = np.arange(l2.n_records)
        plt.figure(facecolor="white")
        plt.scatter(x, ssa_raw)
        plt.scatter(x, smooth(ssa_raw, filter_width), color="red")
        plt.show()


        raise Exception()

#        x = np.arange(l2.n_records)
#        interpolator = interp1d(x[ssh], l2.elev[ssh],
#                                kind='linear', bounds_error=False)
#        self._value = interpolator(x)

    def _calculate_uncertainty(self):
        pass

    def _smooth_interpolated_ssa(self):
        pass


def egm2wgs_delta_h(phi):
    aegm = 6378136.460000
    begm = 6356751.806631
    awgs = 6378137.000000
    bwgs = 6356752.314245
    return compute_delta_h(aegm, begm, awgs, bwgs, phi)


def egm2top_delta_h(phi):
    aegm = 6378136.460000
    begm = 6356751.806631
    atop = 6378136.300000
    btop = 6356751.600563
    return compute_delta_h(aegm, begm, atop, btop, phi)


def compute_delta_h(a1, b1, a2, b2, phi):
    dtor = np.pi/180.0
    delta_a = a2 - a1
    delta_b = b2 - b1
    phir = phi * dtor
    sinsqphi = np.sin(phir)
    sinsqphi = sinsqphi * sinsqphi
    cossqphi = 1.0 - sinsqphi
    return -1.0*(delta_a * cossqphi + delta_b * sinsqphi)


def smooth(x, window):
    """ Numpy implementation of the IDL SMOOTH function """
    from scipy.ndimage.filters import uniform_filter
#    return convolve(x, np.ones(window)/window)
    return uniform_filter(x, size=window)
#    from astropy.convolution import convolve, Box1DKernel
#    kernel = Box1DKernel(window)
#    smoothed = convolve(x, kernel, boundary="extend", normalize_kernel=True)
#    return smoothed
