# -*- coding: utf-8 -*-
"""
Created on Sat Aug 01 17:03:19 2015

@author: Stefan
"""
from pysiral.iotools import ReadNC
from pysiral.filter import (fill_nan, idl_smooth)

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
        self._linear_smoothed_interpolation_between_tiepoints(l2)
        self._calculate_uncertainty(l2)

    def _linear_smoothed_interpolation_between_tiepoints(self, l2):
        """ Based in cs2awi code from Robert Ricker """

        # Use ocean and lead elevations
        self._ssh_tiepoints = l2.surface_type.lead.indices
        if self._options.use_ocean_wfm:
            self._ssh_tiepoints.append(l2.surface_type.ocean.indices)
            self._ssh_tiepoints = np.sort(self._ssh_tiepoints)

        # Get initial elevation at tiepoint locations
        mss_frb = l2.elev - l2.mss

        # Remove ssh tiepoints from the list if their elevation
        # corrected by the median offset of all tiepoints from the mss
        # exceeds a certain threshold
        if self._options.pre_filtering:
            index_dict = np.arange(l2.surface_type.lead.num)
            threshold = self._options.pre_filter_maximum_mss_median_offset
            median_mss_offset = np.median(mss_frb[self._ssh_tiepoints])
            offset = np.abs(mss_frb[self._ssh_tiepoints] - median_mss_offset)
            valid = np.where(offset < threshold)
            self._ssh_tiepoints = self._ssh_tiepoints[index_dict[valid]]

        self.ssa_raw = np.ndarray(shape=(l2.n_records))*np.nan
        self.ssa_raw[self._ssh_tiepoints] = mss_frb[self._ssh_tiepoints]
        non_tiepoints = np.where(np.isnan(self.ssa_raw))

        # Filtered raw values
        # Use custom implementation of IDL SMOOTH:
        # idl_smooth(x, w) equivalent to SMOOTH(x, w, /edge_truncate, /nan)
        ssa_filter1 = idl_smooth(self.ssa_raw, self.filter_width)

        # Leave only the original ssh tie points
        ssa_filter1[non_tiepoints] = np.nan

        # Fill nans with linear interpolation and contant values at borders
        # python: fill_nan(x) = IDL: FILL_NAN(x, /NEIGHBOUR)
        ssa_filter2 = fill_nan(ssa_filter1)

        # Final smoothing
        ssa = idl_smooth(ssa_filter2, self.filter_width)
        self._value = ssa

#        # TODO: Make example plot of individual filter steps
#        import matplotlib.pyplot as plt
#        x = np.arange(l2.n_records)
#        plt.figure(facecolor="white")
#        plt.scatter(x, ssa_raw, color="black")
#        plt.scatter(x, ssa_filter1, color="red", alpha=0.5)
#        plt.plot(x, ssa_filter2, color="blue", lw=2, alpha=0.5)
#        plt.plot(x, ssa, color="orange", lw=3)
#        plt.show()

    def _calculate_uncertainty(self, l2):
        """
        Components that add to sea surface anomaly uncertainty
        - mss uncertainty (if known)
        - uncertainty of lead elevations
        - distance to next lead tiepoint
        """

        # short cuts to options
        max_distance = self._options.uncertainty_tiepoints_distance_max
        ssa_unc_min = self._options.uncertainty_minimum
        ssa_unc_max = self._options.uncertainty_maximum

        # prepare parameter arrays
        lead_indices = l2.surface_type.lead.indices
        lead_elevation = np.full((l2.n_records), np.nan, dtype=np.float32)
        lead_elevation[lead_indices] = self._value[lead_indices]

        # Compute distance to next lead tiepoint in meter
        tiepoint_distance = get_tiepoint_distance(np.isfinite(lead_elevation))
        tiepoint_distance = tiepoint_distance.astype(np.float32)
        tiepoint_distance *= self._options.smooth_filter_width_footprint_size

        # Compute the influence of distance to next tiepoints
        # in the range of 0: minimum influence to 1: maximum influence
        # It is assumed that the uncertainty has a quadratic dependance
        # on tiepoint distance
        tiepoint_distance_scalefact = tiepoint_distance / max_distance
        above_distance_limit = np.where(tiepoint_distance_scalefact > 1.)[0]
        tiepoint_distance_scalefact[above_distance_limit] = 1.
        tiepoint_distance_scalefact = tiepoint_distance_scalefact**2.

        # Compute the ssa uncertainty based on a min/max approach
        # scaled by
        ssa_unc_range = ssa_unc_max - ssa_unc_min
        ssa_unc = ssa_unc_min + ssa_unc_range * tiepoint_distance_scalefact

        # Save in class
        self._uncertainty = ssa_unc

#        XXX: Debug plot code
#        import matplotlib.pyplot as plt
#        plt.figure(figsize=(10, 6))
#        plt.plot(lead_elevation, "o", label="leads")
#        plt.plot(l2.elev-l2.mss, "+", label="elevations above mss")
#        plt.plot(ssa_unc, lw=2, label="ssha uncertainty")
#        plt.legend()
#        plt.ylim(-0.5, 1.0)
#        plt.ylabel("meter")
#        plt.tight_layout()
#        plt.show()


    @property
    def filter_width(self):
        filter_width = self._options.smooth_filter_width_m / \
            self._options.smooth_filter_width_footprint_size
        # Make sure filter width is odd integer
        filter_width = np.floor(filter_width) // 2 * 2 + 1
        filter_width = filter_width.astype(int)
        return filter_width


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


def rolling_window(a, window):
    """
    Recipe for rapid computations of rolling windows using numpy
    from: http://www.rigtorp.se/2011/01/01/rolling-statistics-numpy.html
    """
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)


def get_tiepoint_distance(is_tiepoint):
    """
    Calculates the distance to the next tiepoints (boolean array)
    """
    distance_forward = get_tiepoints_oneway_distance(is_tiepoint)
    distance_reverse = get_tiepoints_oneway_distance(is_tiepoint, reverse=True)
    return np.minimum(distance_forward, distance_reverse)


def get_tiepoints_oneway_distance(a, reverse=False):
    """ loops through array and determines distance to latest flag=true """
    n = len(a)
    distance = np.full(a.shape, n+1, dtype=np.int32)
    if reverse:
        a = a[::-1]
    dist = n
    for i in np.arange(n):
        if a[i]:
            dist = 0
        elif dist == n:
            pass
        else:
            dist += 1
        distance[i] = dist
    if reverse:
        distance = distance[::-1]
    return distance


def get_l2_ssh_class(name):
    return globals()[name]()
