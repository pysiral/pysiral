# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 15:30:53 2016

@author: Stefan
"""

import numpy as np
import bottleneck as bn
from typing import List, Union
from scipy.interpolate import interp1d, UnivariateSpline
from astropy.convolution import convolve, Box1DKernel

import statsmodels.api as sm

from pysiral.l2data import Level2Data
from pysiral.l1bdata import Level1bData
from pysiral.core.flags import FlagContainer, ORCondition
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

    @property
    def error_bit(self):
        return self.error_flag_bit_dict["other"]


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

    @property
    def error_bit(self):
        return self.error_flag_bit_dict["filter"]


class ParameterSmoother(Level2ProcessorStep):
    """
    Creates a filtered/smoothed copy of a given parameter.

    Usage in Level-2 processor definition files:

    -   module: filter
        pyclass: ParameterSmoother
        options:
            source_variable: <source_variable>
            target_variable_name: <output variable auxiliary data name>
            target_variable_id: <output variable auxiliary data id>
            smoothing_method: <the method for the smoothing: box_filter|gaussian_process>
            smoother_args: dict
            retain_nan_mask: bool

    will register a new auxiliary variable output_variable_name|output_variable_id to the
    Level-2 data object
    """
    def __init__(self, *args, **kwargs):
        super(ParameterSmoother, self).__init__(*args, **kwargs)

    def execute_procstep(self, l1b, l2):
        """
        Compute the smoother
        :param l1b:
        :param l2:
        :return:
        """

        # Get the error flag
        error_status = self.get_clean_error_status(l2.n_records)

        # Get the output name and id of the resulting auxiliary data variable
        auxid = self.cfg.options.target_variable_id
        auxname = self.cfg.options.target_variable_name

        # Get a dictionary of the filter function
        filter_func = {"box_filter": self.box_filter_smoother,
                       "lowess": self.lowess_smoother}

        # check if requested smoothing method is implemented
        if self.cfg.options.smoothing_method not in filter_func.keys():
            error_status[:] = True
            msg = "Not-Implemented: Smoothing method: {}".format(self.cfg.options.smoothing_method)
            self.error.add_error("filter-not-implemented", msg)
            l2.register_auxvar(auxid, auxname, np.full(error_status.shape, np.nan), None)
            return error_status

        # Get the source parameter
        var = getattr(l2, self.cfg.options.source_variable)
        result = filter_func[self.cfg.options.smoothing_method](var[:],  **self.cfg.options.smoother_args)
        var_filtered, var_filtered_unc = result

        # Add the result as auxiliary data set without an uncertainty
        l2.set_auxiliary_parameter(auxid, auxname, var_filtered)
        return error_status

    @staticmethod
    def box_filter_smoother(x, window=5):
        return astropy_smooth(np.array(x), window)

    @staticmethod
    def lowess_smoother(y, **smoother_args):
        lowess = sm.nonparametric.lowess
        x = np.arange(y.shape[0])

        # Compute the data fraction to use at each y-value
        # TODO: Value of 1501 is hard coded -> move to the settings
        data_fraction = min(1501./float(y.shape[0]), 1.)

        # Compute the lowess filter with potential keywords from the config file
        filter_props = smoother_args.get("filter_props", {})
        z = lowess(y, x, frac=data_fraction, return_sorted=False, **filter_props)

        return z, None

    @property
    def l2_input_vars(self):
        return [self.cfg.options.source_variable]

    @property
    def l2_output_vars(self):
        return [self.cfg.options.target_variable_id]

    @property
    def error_bit(self):
        return self.error_flag_bit_dict["filter"]


class MarginalIceZoneFilterFlag(Level2ProcessorStep):
    """
    Create a flag value that can be used to filter freeboard/thickness values
    that are affected by surface waves penetrating the marginal ice zone.
    The filter flag is not applied to any geophysical variables, but added
    to the l2 data container.

    The flag can take the following values:

        0: not in marginal ice zone
        1: in marginal ice zone: light to none wave influence detected
        2: in marginal ice zone: medium wave influence detected
        3: in marginal ice zone: severe wave influence detected

    The flag values depend on:

        - leading edge with of ocean waveforms at the ice edge
        - sea ice freeboard gradient as a function of distance to the ice edge
        - sea ice freeboard value

    Thresholds to determine the flag values need to be specified in the
    options of the Level2 processor definition file:

    -   module: filter
        pyclass: MarginalIceZoneFilterFlag
        options:
            distance_to_ocean_maximum: <maximum size of the marginal ice zone>
            leading_edge_width_miz_gradient: [<medium_threshold>, <severe_threshold>]
            leading_edge_width_ocean_value [<medium_threshold>, <severe_threshold>]

    """

    def __init__(self, *args, **kwargs):
        super(MarginalIceZoneFilterFlag, self).__init__(*args, **kwargs)

    def execute_procstep(self,
                         l1b: "Level1bData",
                         l2: "Level2Data"
                         ) -> np.ndarray:
        """
        API method for Level2ProcessorStep subclasses. Computes and add the filter flag to the l2 data object
        :param l1b:
        :param l2:
        :return:
        """

        filter_flag_miz = np.full(l2.n_records, 0, dtype=int)
        filter_flag_miz_error = self.get_clean_error_status(l2.n_records)

        # Step 1: Detect ice edge(s) in the trajectory data.
        ice_edge_idx = self._get_ice_edge_idxs(l2)
        if ice_edge_idx is None:
            l2.set_auxiliary_parameter("fmiz", "filter_flag_miz", filter_flag_miz, None)
            return filter_flag_miz_error

        breakpoint()

    def _get_ice_edge_idxs(self, l2: "Level2Data") -> Union[List, None]:
        """
        Find any transitions between sea ice / open water areas
        :param l2:
        :return: List of indices marking the ice edge (e.g. first ice waveform).
        """

        # Initial sanity check
        if l2.surface_type.ocean.num == 0:
            return None

        import pandas as pd
        # some sample data
        ts = pd.Series(l2.pp)

        # add the 20 day rolling standard deviation:
        # pd.Series(l2.pp).rolling(window=51).skew().plot(style='r')
        ts.plot()
        ts.rolling(window=51, center=True, min_periods=1).std().plot(style='r')
        # ts.rolling(window=51, center=True).min().plot(style='b')
        # ts.rolling(window=51, center=True).max().plot(style='b')

        import matplotlib.pyplot as plt
        # plt.figure()
        # plt.plot(l2.sic)
        # plt.plot(l2.surface_type.ocean.flag*50)
        #
        # plt.figure("Pulse Peakiness")
        # plt.plot(l2.pp)
        #
        # plt.figure("ssd")
        # plt.plot(l2.ssd)
        #
        # plt.figure("lew")
        # plt.plot(l2.lew)

        plt.show()
        breakpoint()

    @property
    def l2_input_vars(self):
        return ["surface_type",
                "leading_edge_width",
                "sea_ice_freeboard",
                "distance_to_ocean"]

    @property
    def l2_output_vars(self):
        return ["filter_flag_miz"]

    @property
    def error_bit(self):
        return self.error_flag_bit_dict["filter"]


def numpy_smooth(x, window):
    return np.convolve(x, np.ones(window)/window)


def scipy_smooth(x, window):
    """ Numpy implementation of the IDL SMOOTH function """
    from scipy.ndimage.filters import uniform_filter
    return uniform_filter(x, size=window)


def astropy_smooth(x, window, boundary="extend", normalize_kernel=True):
    kernel = Box1DKernel(window)
    smoothed = convolve(x, kernel, boundary=boundary, normalize_kernel=normalize_kernel)
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
        smoothed[i] = bn.nanmean(x[i-kernel_halfsize:i+kernel_halfsize+1])
    return smoothed


def spline_smooth(y, window):
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
