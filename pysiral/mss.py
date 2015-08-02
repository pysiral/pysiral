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
        self._linear_smoothed_interpolation_between_tiepoints(l2)
        self._calculate_uncertainty()

    def _linear_smoothed_interpolation_between_tiepoints(self, l2):
        # Use ocean and lead elevations
        self._ssh_tiepoints = l2.surface_type.lead.indices
        if self._options.use_ocean_wfm:
            self._ssh_tiepoints.append(l2.surface_type.ocean.indices)
            self._ssh_tiepoints = np.sort(self._ssh_tiepoints)
        # Get initial elevation at tiepoint locations
        mss_frb = l2.elev - l2.mss
        ssa_raw = np.ndarray(shape=(l2.n_records))*np.nan
        ssa_raw[self._ssh_tiepoints] = mss_frb[self._ssh_tiepoints]
        non_tiepoints = np.where(np.isnan(ssa_raw))
        # Prepare filtering, get filter width
        #   filter width need to in points
        #   -> need typical footprint size to be defined in filter options
        filter_width = self._options.smooth_filter_width_m / \
            self._options.smooth_filter_width_footprint_size
        # Make sure filter width is odd number
        filter_width = np.floor(filter_width) // 2 * 2 + 1
        filter_width = filter_width.astype(int)
        # Filtered raw values
        # Use custom implementation of IDL SMOOTH:
        # idl_smooth(x, w) equivalent to SMOOTH(x, w, /edge_truncate, /nan)
        ssa_filter1 = idl_smooth(ssa_raw, filter_width)
        # Leave only the original ssh tie points
        ssa_filter1[non_tiepoints] = np.nan
        # Fill nans with linear interpolation and contant values at borders
        # python: fill_nan(x) = IDL: FILL_NAN(x, /NEIGHBOUR)
        ssa_filter2 = fill_nan(ssa_filter1)
        # Final smoothing
        ssa = idl_smooth(ssa_filter2, filter_width)
        self._value = ssa

        # TODO: Make example plot of individual filter steps
#        import matplotlib.pyplot as plt
#        x = np.arange(l2.n_records)
#        plt.figure(facecolor="white")
#        plt.scatter(x, ssa_raw, color="black")
#        plt.scatter(x, ssa_filter1, color="red", alpha=0.5)
#        plt.plot(x, ssa_filter2, color="blue", lw=2, alpha=0.5)
#        plt.plot(x, ssa, color="orange", lw=3)
#        plt.show()

    def _calculate_uncertainty(self):
        pass

#      IDL Code for getting ssh uncertainty
#      STD = REPLICATE(!VALUES.F_NAN, CS_L1B.N_RECORDS)
#      FOR I=0,CS_L1B.N_RECORDS-1 DO BEGIN
#        IF I LT SM/2 THEN BEGIN
#          TMP = WHERE(FINITE(SSHA_RAW[0:I+SM/2]) EQ 1, NN)
#          IF NN GT 0 THEN STD[I] = 1/SQRT(NN) * STDDEV([REPLICATE(!VALUES.D_NAN,ABS(I-SM/2)),SSHA_RAW[0:I+SM/2]],/NAN)
#          NN = 0
#        ENDIF
#        IF I GE SM/2 AND I LE CS_L1B.N_RECORDS-1-SM/2 THEN BEGIN
#          TMP = WHERE(FINITE(SSHA_RAW[I-SM/2:I+SM/2]) EQ 1, NN)
#          IF NN GT 0 THEN STD[I] = 1/SQRT(NN) * STDDEV(SSHA_RAW[I-SM/2:I+SM/2],/NAN)
#          NN = 0
#        ENDIF
#        IF I GT CS_L1B.N_RECORDS-1-SM/2 THEN BEGIN
#          TMP = WHERE(FINITE(SSHA_RAW[I-SM/2:CS_L1B.N_RECORDS-1]) EQ 1, NN)
#          IF NN GT 0 THEN STD[I] = 1/SQRT(NN) * STDDEV([SSHA_RAW[I-SM/2:CS_L1B.N_RECORDS-1],REPLICATE(!VALUES.D_NAN,ABS(CS_L1B.N_RECORDS-1+SM/2))],/NAN)
#          NN = 0
#        ENDIF
#      ENDFOR
#
#      TMP = WHERE(FINITE(STD,/NAN),NN)
#      IF NN GT 0 THEN BEGIN
#        STD_MAX = ABS(SSHA_VALID - SMOOTH(ELEVATION,SM,/NAN,/EDGE_TRUNCATE))
#        STD[TMP] = STD_MAX[TMP]
#      ENDIF
#
#      SSH_UNCERTAINTY = SMOOTH(STD,SM,/NAN,/EDGE_TRUNCATE)


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


def numpy_smooth(x, window):
    return np.convolve(x, np.ones(window)/window)


def scipy_smooth(x, window):
    """ Numpy implementation of the IDL SMOOTH function """
    from scipy.ndimage.filters import uniform_filter
    return uniform_filter(x, size=window)


def astropy_smooth(x, window):
    from astropy.convolution import convolve, Box1DKernel
    kernel = Box1DKernel(window)
    smoothed = convolve(x, kernel, boundary="extend", normalize_kernel=True)
    return smoothed


def idl_smooth(x, window):
    """ Implementation of the IDL smooth(x, window, /EDGGE_TRUNCATE, /NAN)"""
    smoothed = np.copy(x)*np.nan
    n = len(x)
    for i in np.arange(n):
        kernel_halfsize = np.floor((window-1)/2).astype(int)
        if (i < kernel_halfsize):
            kernel_halfsize = i
        if (n-1-i < kernel_halfsize):
            kernel_halfsize = n-1-i
        smoothed[i] = np.nanmean(x[i-kernel_halfsize:i+kernel_halfsize+1])
    return smoothed


# TODO: Custom implementation of a gaussian weighted smoother?
def gauss_smooth(x, window):
    pass
#    from astropy.convolution import Gaussian1DKernel
#    smoothed = np.copy(x)*np.nan
#    # get a gaussi
#    kernel = Gaussian1DKernel(10)
#    n = len(x)
#    for i in np.arange(n):
#        kernel_halfsize = np.floor((window-1)/2).astype(int)
#        if (i < kernel_halfsize):
#            kernel_halfsize = i
#        if (n-1-i < kernel_halfsize):
#            kernel_halfsize = n-1-i
#        smoothed[i] = np.nanmean(x[i-kernel_halfsize:i+kernel_halfsize+1])
#    return smoothed


def spline_smooth(y, window):
    from scipy.interpolate import UnivariateSpline
    x = np.arange(len(y))
    no_nan = np.where(np.isfinite(y))
    spl = UnivariateSpline(x[no_nan], y[no_nan], s=len(y)/window, k=4)
    spl.set_smoothing_factor(0.25)
    return spl(x)


def fill_nan(y):
    result = np.copy(y)
    # Get first and last valid index
    no_nan = np.where(np.isfinite(y))[0]
    valid0, valid1 = np.amin(no_nan), np.amax(no_nan)
    # Cut the "inside" section that is bounded by valid measurements
    y_inside = y[valid0:valid1+1]
    x = np.arange(len(y_inside))
    valid_inside = np.where(np.isfinite(y_inside))[0]
    # Interpolate inside range of valid entries
    func = interp1d(x[valid_inside], y_inside[valid_inside],
                    bounds_error=False)
    interpolated_inside = func(x)
    result[valid0:valid1+1] = interpolated_inside
    # fill nan-borders with first/last valid value
    if valid0 != 0:
        result[0:valid0] = y[valid0]
    if valid1 != len(y)-1:
        result[valid1+1:len(y)+1] = y[valid1]
    return result
