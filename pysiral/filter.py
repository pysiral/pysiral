# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 15:30:53 2016

@author: Stefan
"""

import numpy as np
from scipy.interpolate import interp1d
from astropy.convolution import convolve

from pysiral.flag import FlagContainer, ORCondition
from pysiral.l2proc.procsteps import Level2ProcessorStep


class L1bEnvisatBackscatterDriftCorrection(Level2ProcessorStep):
    """
    Very specific filter to correct backscatter drift. The filter
    applies a monthly linear correction based on the drift factor and
    base period.

    Notes:

        1. This class is very specifically designed for correcting
           the sigma0 value of the degrading Envisat RA-2 antenna
    """

    def __init__(self, *args, **kwargs):
        super(L1bEnvisatBackscatterDriftCorrection, self).__init__(*args, **kwargs)

    def execute_procstep(self, l1b, l2):
        """
        Apply a backscatter correction
        :param l1b:
        :param l2:
        :return:
        """

        # Get the default error status
        error_status = self.get_clean_error_status(l2.n_records)

        # Get the backcatter value
        datagroup = self.cfg.options.l1b_data_group
        name = self.cfg.options.l1b_parameter_name
        value = l1b.get_parameter_by_name(datagroup, name)

        # Compute the drift correction
        year_base, month_base = self.cfg.options.backscatter_base_period
        year, month = l1b.info.start_time.year, l1b.info.start_time.month
        time_shift_factor = (year_base - year) * 12 + (month_base - month)
        sigma0_drift_factor = self.cfg.options.backscatter_drift_factor
        sigma0_drift = time_shift_factor * sigma0_drift_factor

        # Apply and update the backscatter by changing the l1b data object in-place
        value += sigma0_drift
        l1b.set_parameter_by_name(datagroup, name, value)

        # Done
        return error_status

    @property
    def l2_input_vars(self):
        return []

    @property
    def l2_output_vars(self):
        return []


class L2ParameterValidRange(Level2ProcessorStep):
    """
    Filters L2 variables that are outside a valid range.

    Usage in Level-2 processor definition files:

    -   module: filter
        pyclass: L2ParameterValidRange
        options:
            source_variable: <source_variable>
            target_variables: [<target_variables>]
            valid_minimum_point_value: <min_val>
            valid_maximum_point_value: <max_val>

    will lead to that all target_variables will be set to np.nan if either

        source_variable < min_val  of source_variable > max_val

    """

    def __init__(self, *args, **kwargs):
        super(L2ParameterValidRange, self).__init__(*args, **kwargs)

    def execute_procstep(self, l1b, l2):
        """
        Mandatory method for Level2ProcessorStep will change l2 data object in place
        :param l1b:
        :param l2:
        :return:
        """

        # Step 1: Get the filter flag
        parameter = getattr(l2, self.cfg.options.source_variable)
        invalid = ORCondition()
        invalid.add(parameter < self.cfg.options.valid_minimum_point_value)
        invalid.add(parameter > self.cfg.options.valid_maximum_point_value)
        filter_flag = FlagContainer(invalid.flag)

        # Check if any flag was raised
        if filter_flag.num == 0:
            return filter_flag.flag

        # Modify target variables
        for parameter_id in self.cfg.options.target_variables:
            var = getattr(l2, parameter_id)
            var.set_nan_indices(filter_flag.flag)
            setattr(l2, parameter_id, var)

        # The filter flag can be used as error_status
        return filter_flag.flag

    @property
    def l2_input_vars(self):
        input_vars = [self.cfg.options.source_variable]
        input_vars.extend(self.cfg.options.target_variables)
        return input_vars

    @property
    def l2_output_vars(self):
        return []


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
        if i < kernel_halfsize:
            kernel_halfsize = i
        if (n-1-i) < kernel_halfsize:
            kernel_halfsize = n-1-i
        smoothed[i] = np.nanmean(x[i-kernel_halfsize:i+kernel_halfsize+1])
    return smoothed


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


def smooth_2darray(array, filter_width=5, preserve_gaps=True):
    """ A simple smoothing filter """

    # presever nan's and mask
    nan_list = np.where(np.isnan(array))

    try:
        mask = array.mask
    except AttributeError:
        mask = np.zeros(array.shape)

    # astropy concolve seems to be the only one treating nan's
    # as missing values
    kernel = np.ones((filter_width, filter_width))
    convolved = convolve(array, kernel, normalize_kernel=True)
    output = np.ma.array(convolved)

    if preserve_gaps:
        output[nan_list] = np.nan
        output.mask = mask

    return output
