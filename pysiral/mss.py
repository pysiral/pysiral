# -*- coding: utf-8 -*-
"""
Created on Sat Aug 01 17:03:19 2015

@author: Stefan
"""

from pysiral.errorhandler import ErrorStatus
from pysiral.auxdata import AuxdataBaseClass
from pysiral.iotools import ReadNC
from pysiral.filter import (fill_nan, idl_smooth)

from treedict import TreeDict
import scipy.ndimage as ndimage
import numpy as np


class BaseMSS(AuxdataBaseClass):

    def __init__(self):
        super(BaseMSS, self).__init__()
        self._roi = None

    def set_roi(self, roi):
        self._roi = roi

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

    def _initialize(self):

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
        self.error = ErrorStatus()

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
        self.error.caller_id = self.__class__.__name__

    def _interpolate(self, l2):
        self._linear_smoothed_interpolation_between_tiepoints(l2)
        self._calculate_uncertainty(l2)
        if "marine_segment_filter" in self._options:
            self._marine_segment_filter(l2)
        if "tiepoint_maxdist_filter" in self._options:
            self._tiepoint_maxdist_filter(l2)

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


##        # TODO: Make example plot of individual filter steps
#        import matplotlib.pyplot as plt
#        x = np.arange(l2.n_records)
#
#        plt.figure("land")
#        plt.plot(l2.surface_type.land.flag)
#
#        plt.figure()
#        plt.plot(x, l2.mss)
#        plt.plot(x, l2.mss+ssa)
#        plt.scatter(x, l2.elev, marker="+", alpha=0.5)
#
#        plt.figure(facecolor="white")
#        # plt.scatter(x, ssa_raw, color="black")
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

        # get tiepoint distance
        tiepoint_distance = self.get_tiepoint_distance(l2)

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

    def _marine_segment_filter(self, l2):
        """ Check all sections divided by land masses for reliable
        information content """

        filter_options = self._options.marine_segment_filtering
        minimum_lead_number = filter_options.minimum_lead_number
        footprint_size = self._options.smooth_filter_width_footprint_size
        section_prop = {"i0": 0.0, "i1": 0.0,
                        "width": 0.0, "n_tiepoints": 0,
                        "land_before": (9999.0, 0),
                        "land_after": (9999.0, 0)}

        # Find sea ice clusters
        land = l2.surface_type.land

        # No land -> nothing to do
        if land.num == 0:
            return

        # Get indices for land sections
        lead_flag = l2.surface_type.lead.flag
        land_flag = land.flag.astype(int)
        land_start = np.where(np.ediff1d(land_flag) > 0)[0]
        land_stop = np.where(np.ediff1d(land_flag) < 0)[0]

        # It is assumed here, that the l1b orbit segment never starts
        # or end with land. Thus the number of land start and land stop
        # events need to be identical.
        if len(land_start) != len(land_stop):
            code = "l2-crop-error"
            msg = "l2 segments either starts or ends with land"
            self.error.add_error(code, msg)
            return

        # Add artificial large land sections on beginning and end of profile
        n_marine_segments = len(land_start) + 1
        n = l2.n_records
        land_start = np.concatenate(([-1000], land_start, [n-1]))
        land_stop = np.concatenate(([-1], land_stop, [n+1000]))

        # Loop over marine segments and collect information
        marine_segments = []
        for i in np.arange(n_marine_segments):

            marine_segment = section_prop.copy()

            # Get the start stop indices for marine section
            i0 = land_stop[i]+1
            i1 = land_start[i+1]
            marine_segment["i0"] = i0
            marine_segment["i1"] = i1
            marine_segment["width"] = (i1-i0) * footprint_size

            # get the number of leads
            marine_section_indices = np.arange(i0, i1+1)
            n_tiepoints = np.where(lead_flag[marine_section_indices])[0].size
            marine_segment["n_tiepoints"] = n_tiepoints

            if marine_segment["n_tiepoints"] < minimum_lead_number:
                self._value[marine_section_indices] = np.nan

            marine_segments.append(marine_segment)

    def _tiepoint_maxdist_filter(self, l2):
        """  A filter that does not removes ssa values which distance to
        the next ssh tiepoint exceeds a defined threshold """

        # Get options
        filter_options = self._options.tiepoint_maxdist_filter
        edges_only = filter_options.edges_only
        distance_threshold = filter_options.maximum_distance_to_tiepoint

        # Compute distance to next tiepoint
        tiepoint_distance = self.get_tiepoint_distance(l2)

        # Get indices
        invalid_indices = np.where(tiepoint_distance > distance_threshold)[0]

        # Only remove ssa values at the edges (if edges_only:True)
        if edges_only:
            lead_indices = l2.surface_type.lead.indices
            before_first_lead = invalid_indices < lead_indices[0]
            after_last_lead = invalid_indices > lead_indices[-1]
            edge_condition = np.logical_or(before_first_lead, after_last_lead)
            invalid_indices = invalid_indices[np.where(edge_condition)[0]]

        # Validity check
        if len(invalid_indices) == 0:
            return

        # Set the values
        self._value[invalid_indices] = np.nan

    def get_tiepoint_distance(self, l2):
        """ Returns the distance in meter to the next ssh tiepoint for
        each record """

        # prepare parameter arrays
        lead_indices = l2.surface_type.lead.indices
        lead_elevation = np.full((l2.n_records), np.nan, dtype=np.float32)
        lead_elevation[lead_indices] = self._value[lead_indices]

        # Compute distance to next lead tiepoint in meter
        tiepoint_distance = get_tiepoint_distance(np.isfinite(lead_elevation))
        tiepoint_distance = tiepoint_distance.astype(np.float32)
        tiepoint_distance *= self._options.smooth_filter_width_footprint_size
        return tiepoint_distance

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
    pyclass = globals().get(name, None)
    if pyclass is not None:
        return pyclass()
    else:
        return pyclass
