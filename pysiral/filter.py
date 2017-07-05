# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 15:30:53 2016

@author: Stefan
"""
from pysiral.config import options_from_dictionary
from pysiral.logging import DefaultLoggingClass
from pysiral.flag import FlagContainer, ORCondition

from scipy.interpolate import interp1d
import numpy as np


class FilterBaseClass(DefaultLoggingClass):

    def __init__(self):
        super(FilterBaseClass, self).__init__(self.__class__.__name__)
        self._flag = None

    def set_options(self, **opt_dict):
        self._options = options_from_dictionary(**opt_dict)

    def apply_filter(self, *args, **kwargs):
        self._apply_filter(*args, **kwargs)

    @property
    def flag(self):
        return self._flag

    @property
    def options(self):
        return self._options


# %% Filter for Level2Processor

class L1bBackscatterDriftCorrection(FilterBaseClass):
    """ Very specific filter to correct backscatter drift. The filter
    applies a monthly linear correction based on the drift factor and
    base period. """

    def __init__(self):
        super(L1bBackscatterDriftCorrection, self).__init__()
        self.log.name = self.__class__.__name__

    def _apply_filter(self, l1b):
        """ Apply the backscatter correction """
        # Get the backscatter value
        datagroup = self.options.l1b_data_group
        name = self.options.l1b_parameter_name
        value = l1b.get_parameter_by_name(datagroup, name)
        # Compute the drift correction
        year_base, month_base = self.options.backscatter_base_period
        year, month = l1b.info.start_time.year, l1b.info.start_time.month
        time_shift_factor = (year_base - year) * 12 + \
            (month_base - month)
        sigma0_drift_factor = self.options.backscatter_drift_factor
        sigma0_drift = time_shift_factor * sigma0_drift_factor
        # Apply and update the backscatter
        value += sigma0_drift
        l1b.set_parameter_by_name(datagroup, name, value)


class L2ParameterValidRange(FilterBaseClass):
    """ Filters freeboard outliers by simple min/max thresholding
    Requires l2 data container and target (either: "afrb", "rfrb")
    """

    def __init__(self):
        super(L2ParameterValidRange, self).__init__()

    def _apply_filter(self, l2, target):
        parameter = getattr(l2, target)
        invalid = ORCondition()
        invalid.add(parameter < self._options.valid_minimum_point_value)
        invalid.add(parameter > self._options.valid_maximum_point_value)
        self._flag = FlagContainer(invalid.flag)


def get_filter(name):
    return globals()[name]()


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
    if len(no_nan) == 0:
        return result
    valid0, valid1 = np.amin(no_nan), np.amax(no_nan)
    # Cut the "inside" section that is bounded by valid measurements
    y_inside = y[valid0:valid1+1]
    x = np.arange(len(y_inside))
    valid_inside = np.where(np.isfinite(y_inside))[0]
    # Interpolate inside range of valid entries
    try:
        func = interp1d(x[valid_inside], y_inside[valid_inside],
                        bounds_error=False)
        interpolated_inside = func(x)
        result[valid0:valid1+1] = interpolated_inside
    except ValueError:
        # May not be applicable (no inner nan ranges)
        pass
    # fill nan-borders with first/last valid value
    if valid0 != 0:
        result[0:valid0] = y[valid0]
    if valid1 != len(y)-1:
        result[valid1+1:len(y)+1] = y[valid1]
    return result
