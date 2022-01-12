# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 15:30:53 2016

@author: Stefan
"""
import matplotlib.pyplot as plt
import numpy as np
from dataclasses import dataclass
import bottleneck as bn
from typing import List, Union, Tuple
from loguru import logger
from scipy.interpolate import interp1d, UnivariateSpline
from astropy.convolution import convolve, Box1DKernel

import statsmodels.api as sm

from pysiral.l2data import Level2Data
from pysiral.l1bdata import Level1bData
from pysiral.core.flags import ANDCondition
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


@dataclass
class MarginalIceZoneProperties:
    """ A lean data class for MarginalIceZoneFilterFlag """
    ice_edge_idx: int
    sea_ice_is_left: bool
    ocean_idxs: np.ndarray
    seaice_idxs: np.ndarray
    is_valid: bool = False
    ocean_leading_edge_width: float = None
    sea_ice_freeboard_values: np.ndarray = None
    sea_ice_freeboard_gradient: np.ndarray = None
    flag_value: np.ndarray = None


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

        filter_flag_miz = np.full(l2.n_records, -1, dtype=int)
        filter_flag_miz_error = self.get_clean_error_status(l2.n_records)

        # Check if there is any ocean data at all.
        if l2.surface_type.ocean.num == 0:
            l2.set_auxiliary_parameter("fmiz", "flag_miz", filter_flag_miz, None)
            return filter_flag_miz_error

        if not np.isfinite(l2.frb[:]).any():
            l2.set_auxiliary_parameter("fmiz", "flag_miz", filter_flag_miz, None)
            return filter_flag_miz_error

        # Step 1: Detect ice edge(s) in the trajectory data.
        miz_prop_list = self.get_ice_edge_idxs(l2)
        if miz_prop_list is None:
            l2.set_auxiliary_parameter("fmiz", "flag_miz", filter_flag_miz, None)
            return filter_flag_miz_error

        # Step 2: Determine the waveform properties on both sides of the ice edge.
        # NOTE: There might be multiple transitions per trajectory
        for miz_prop in miz_prop_list:
            properties = self.get_miz_data_properties(l2, miz_prop)
            if properties is None:
                l2.set_auxiliary_parameter("fmiz", "flag_miz", filter_flag_miz, None)
                return filter_flag_miz_error
            l2.set_auxiliary_parameter("fmiz", "flag_miz", properties.flag_value, None)
        return filter_flag_miz_error

    def get_ice_edge_idxs(self, l2: "Level2Data") -> Union[List["MarginalIceZoneProperties"], None]:
        """
        Find any transitions between sea ice / open water areas determined
        by sea ice concentration crossing the 15% threshold.
        :param l2: Level2Data
        :return: List of indices marking the ice edge (e.g. first ice waveform).
        """

        # Initial sanity check

        # Get properties
        sea_ice_filter_size_m = self.cfg.options.get("sea_ice_filter_size_m", None)
        if sea_ice_filter_size_m is None:
            logger.error(f"{self.__class__.__name__}: Missing option `sea_ice_filter_size_m")
            return None

        ocean_filter_size_m = self.cfg.options.get("ocean_filter_size_m", None)
        if ocean_filter_size_m is None:
            logger.error(f"{self.__class__.__name__}: Missing option `ocean_leading_edge_window")
            return None

        # Detect change ice edge by sea ice concentration jumps over the 15% threshold
        sea_ice_concentration = l2.sic[:]
        distance_to_ocean = l2.get_parameter_by_name("distance_to_ocean")
        sic_flag = np.array(sea_ice_concentration > 15.0).astype(float)
        sic_flag_change = np.abs(np.ediff1d(sic_flag))
        ice_edge_candidate_idxs = np.where(sic_flag_change > 0.)[0]

        # Filter indices that are resulting of NaN values and return
        ice_edge_idxs = [
            idx
            for idx in ice_edge_candidate_idxs
            if np.all(np.isfinite([sea_ice_concentration[idx-1: idx+2]]))
        ]

        # Find out which side is which and get ocean/ice data range
        ocean_filter_size_float = ocean_filter_size_m / l2.footprint_spacing
        ocean_filter_size = int(int(ocean_filter_size_float) // 2 * 2 + 1)

        output = []
        for ice_edge_idx in ice_edge_idxs:

            sea_ice_is_left = sea_ice_concentration[ice_edge_idx+1] < 15.

            # sea ice indices are determined by the distance to the ice edge
            # TODO: Add a track distance cut-off
            sea_ice_idxs = np.where(
                np.logical_and(
                    distance_to_ocean <= sea_ice_filter_size_m,
                    np.isfinite(sea_ice_concentration)
                )
            )[0]

            # Ocean indices are determined by window next to the ice edge
            if sea_ice_is_left:
                ocean_idxs = np.arange(
                    ice_edge_idx+1,
                    min(ice_edge_idx + ocean_filter_size, l2.n_records)
                )
            else:
                ocean_idxs = np.arange(
                    max(ice_edge_idx - ocean_filter_size, 0),
                    ice_edge_idx+1
                )

            output.append(MarginalIceZoneProperties(ice_edge_idx,
                                                    sea_ice_is_left,
                                                    ocean_idxs,
                                                    sea_ice_idxs)
                          )

            # plt.figure()
            # plt.plot(sea_ice_concentration)
            # plt.xlim(0, l2.n_records-1)
            #
            # plt.figure()
            # plt.plot(distance_to_ocean)
            # plt.xlim(0, l2.n_records-1)
            #
            # plt.figure()
            # plt.scatter(sea_ice_idxs, sea_ice_concentration[sea_ice_idxs], label="sea ice")
            # plt.scatter(ocean_idxs, sea_ice_concentration[ocean_idxs], label="ocean")
            # plt.legend()
            # plt.xlim(0, l2.n_records - 1)
            #
            # plt.figure()
            # plt.scatter(sea_ice_idxs, distance_to_ocean[sea_ice_idxs], label="sea ice")
            # plt.scatter(ocean_idxs, distance_to_ocean[ocean_idxs], label="ocean")
            # plt.legend()
            # plt.xlim(0, l2.n_records - 1)
            #
            # plt.show()

        return output

    def get_miz_data_properties(self,
                                 l2: "Level2Data",
                                 miz_prop: "MarginalIceZoneProperties"
                                 ) -> Union["MarginalIceZoneProperties", None]:
        """
        Compute the data properties in the marginal ice zone for the filter classification
        :param l2:
        :param miz_prop:
        :return:
        """

        # distance_to_ocean = l2.get_parameter_by_name("distance_to_ocean")
        # # sea_ice_concentration = l2.get_parameter_by_name("sea_ice_concentration")
        sea_ice_freeboard = l2.get_parameter_by_name("sea_ice_freeboard")
        # sea_ice_freeboard_rolling_sdev = l2.get_parameter_by_name("sea_ice_freeboard_rolling_sdev")
        leading_edge_width = l2.get_parameter_by_name("leading_edge_width")

        # First test: There must be a valid ocean leading edge width and freeboard estimate
        miz_frbs = sea_ice_freeboard[miz_prop.seaice_idxs]
        ocean_lews = leading_edge_width[miz_prop.ocean_idxs]
        if not miz_frbs.any() or not ocean_lews.any():
            return None

        # The leading edge value at the ocean side of the ice edge
        ocean_lew_at_ice_edge = np.nanmean(ocean_lews)

        # Compute the gradient of filtered sea ice freeboard as change / km
        # which is expected to be strong in the presence of wave influence
        x_all = np.arange(0, l2.n_records)
        freeboard_smoother_filter_size_m = self.cfg.options.get("freeboard_smoother_filter_size_m", None)
        if freeboard_smoother_filter_size_m is None:
            logger.error(f"{self.__class__.__name__}: Missing option `freeboard_smoother_filter_size_m")
            return None

        # Compute a smoothed & gap-less freeboard time series and its gradient
        freeboard_smoother_filter_size_float = freeboard_smoother_filter_size_m / l2.footprint_spacing
        freeboard_smoother_filter_size = int(int(freeboard_smoother_filter_size_float) // 2 * 2 + 1)
        frb_flt, frb_flt_gradient = self.get_filtered_value_and_gradient(
            x_all,
            sea_ice_freeboard,
            freeboard_smoother_filter_size,
            gradient_scale_factor=1000./l2.footprint_spacing
        )

        # Find the gradient value next to the ice edge
        valid_frb_idx = np.where(np.isfinite(frb_flt_gradient))[0]
        next_to_ice_edge_idx = valid_frb_idx[-1] if miz_prop.sea_ice_is_left else valid_frb_idx[1]
        frb_gradient_at_ice_edge = frb_flt_gradient[next_to_ice_edge_idx]

        # Compute the filter flag
        filter_flag = np.full(l2.n_records, 0).astype(int)

        # Filter flag 1: Is in marginal ice zone
        # All waveforms on the sea ice side within a specified distance to the ocean
        # NOTE: From SIC based distance_to_ocean, not the along-track distance
        filter_flag[miz_prop.seaice_idxs] = 1

        frb_gradient_threshold = self.cfg.options.get("sea_ice_freeboard_miz_gradient")
        lew_ocean_threshold = self.cfg.options.get("leading_edge_width_ocean_value")

        # Find first zero passing of freeboard gradient next to ice edge
        filter_flag_idx = np.array([])
        condition = (
                np.abs(frb_gradient_at_ice_edge) >= frb_gradient_threshold
                or
                ocean_lew_at_ice_edge >= lew_ocean_threshold
        )

        if condition:
            filter_flag_idx = self.get_impacted_range(frb_flt_gradient,
                                                      next_to_ice_edge_idx,
                                                      miz_prop.sea_ice_is_left)
            filter_flag[filter_flag_idx] = 2
        miz_prop.flag_value = filter_flag

        # fig = plt.figure()
        # figManager = plt.get_current_fig_manager()
        # figManager.window.showMaximized()
        #
        # ax1 = fig.add_subplot(2, 2, 1)
        # ax1.set_title(f"Sea Ice Freeboard [is_left: {miz_prop.sea_ice_is_left}]")
        # ax1.scatter(x_all, sea_ice_freeboard, s=1, linewidths=0)
        # if filter_flag_idx.size != 0:
        #     ax1.scatter(x_all[filter_flag_idx], sea_ice_freeboard[filter_flag_idx], s=8, linewidths=0)
        # ax1.plot(x_all, frb_flt)
        #
        # ax1.set_xlim(0, l2.n_records)
        # ax1.axvline(miz_prop.ice_edge_idx, lw=0.5, c="0.0", ls="--")
        # ax1.set_ylim(-0.2, 1.0)
        #
        # ax2 = fig.add_subplot(2, 2, 2)
        # ax2.set_title(f"ocean_lew_at_ice_edge={ocean_lew_at_ice_edge}")
        # ax2.scatter(x_all, leading_edge_width, s=1, linewidths=0)
        # ax2.scatter(x_all[miz_prop.ocean_idxs], leading_edge_width[miz_prop.ocean_idxs], s=8, linewidths=0)
        # ax2.axvline(miz_prop.ice_edge_idx, lw=0.5, c="0.0", ls="--")
        # ax2.axhline(ocean_lew_at_ice_edge)
        # ax2.axhline(lew_ocean_threshold, lw=0.5, c="0.0", ls="--")
        # ax2.set_xlim(0, l2.n_records)
        # ax2.set_ylim(0, 10)
        #
        # ax5 = fig.add_subplot(2, 2, 3)
        # ax5.set_title(f"Sea Ice Freeboard Gradient [{frb_gradient_at_ice_edge}]")
        # ax5.scatter(x_all, frb_flt_gradient, s=1, linewidths=0)
        # # ax5.scatter(x_all[miz_prop.seaice_idxs], frb_flt_gradient[miz_prop.seaice_idxs],
        # #             s=8, linewidths=0)
        # ax5.axhline(frb_gradient_at_ice_edge, lw=1)
        # ax5.axhline(frb_gradient_threshold, lw=0.5, c="0.0", ls="--")
        # ax5.axhline(-1*frb_gradient_threshold, lw=0.5, c="0.0", ls="--")
        # ax5.axvline(miz_prop.ice_edge_idx, lw=0.5, c="0.0", ls="--")
        # ax5.axvline(valid_frb_idx[0])
        # ax5.axvline(valid_frb_idx[-1])
        # ax5.set_xlim(0, l2.n_records)
        #
        # ax6 = fig.add_subplot(2, 2, 4)
        # ax6.set_title("MIZ Filter Flag")
        # ax6.scatter(x_all, filter_flag, s=1, linewidths=0)
        # ax6.axvline(miz_prop.ice_edge_idx, lw=0.5, c="0.0", ls="--")
        # ax6.set_xlim(0, l2.n_records)
        # ax6.set_ylim(-1, 4)
        #
        # plt.show()
        return miz_prop

    @staticmethod
    def get_filtered_value_and_gradient(x: np.ndarray,
                                        y: np.ndarray,
                                        filter_size: int,
                                        gradient_scale_factor: float = 1.0,
                                        ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute a filtered version and its gradient of a noisy and gappy variable
        (In this case likely freeboard)
        :param x:
        :param y:
        :param filter_size:
        :param gradient_scale_factor:
        :return:
        """

        # Step 1: Use a lowess filter to create
        y_filtered = lowess_smooth(x, y, filter_size)

        # Step 2: Fill NaN Gaps in gradient with linear interpolation
        y_filtered_gapless = interp1d_gap_filling(y_filtered)

        # Step 3: Remove minor wiggles which may happen in the lowess filtering
        #         Also ensure that the edges of the valid data are not flattened
        post_filter_size = 11
        y_filtered2 = astropy_smooth(y_filtered_gapless, post_filter_size, boundary="fill", fill_value=np.nan)
        valid_idx = np.where(np.isfinite(y_filtered2))[0]
        i0, i1 = valid_idx[0], valid_idx[-1]
        y_filtered2[i0:i0+post_filter_size] = y_filtered_gapless[i0:i0+post_filter_size]
        y_filtered2[i1-post_filter_size:i1+1] = y_filtered_gapless[i1-post_filter_size:i1+1]

        # Use scaling for gradient (meant to convert to m/m; m/km, etc).
        y_gradient = np.ediff1d(y_filtered2) * gradient_scale_factor
        y_gradient = np.insert(y_gradient, 0, np.nan)

        # plt.figure(dpi=150)
        # # plt.scatter(x, y, linewidths=0)
        # plt.plot(x, y_filtered, lw=3, label="y_filtered")
        # plt.plot(x, y_filtered_gapless, lw=2, label="y_filtered_gapless")
        # plt.plot(x, y_filtered2, lw=1, label="y_filtered2")
        # plt.scatter(x, y_filtered2, s=1, label="y_filtered2 points")
        # plt.legend()
        #
        # plt.figure(dpi=150)
        # plt.plot(x, y_gradient)
        # plt.show()

        return y_filtered2, y_gradient

    @staticmethod
    def get_impacted_range(frb_gradient: np.ndarray,
                           next_to_ice_edge_idx: int,
                           sea_ice_is_left: bool) -> np.ndarray:
        # sourcery skip: inline-immediately-returned-variable
        """
        Compute the impacted area which is defined as the range between the nearest point
        to the ice edge to the point where the sea ice freeboard gradient show the first
        zero crossing (-> reversal of freeboard trend). This method should only be called,
        when it is quite certain that the freeboard is impacted by wave.
        :param frb_gradient:
        :param next_to_ice_edge_idx:
        :param sea_ice_is_left:
        :return: Indices Array (potentially empty)
        """

        # Compute zero crossings (and discard the last incorrect one)
        # NOTE: without the `np.isfinite` check there will be a lot false zero crossings.
        zero_crossing_idx = np.where(
            np.logical_and(
                np.diff(np.sign(frb_gradient)).astype(bool),
                np.isfinite(frb_gradient[:-1]),
            )
        )[0]
        if zero_crossing_idx.size <= 1:
            return np.array([])
        zero_crossing_idx = zero_crossing_idx[:-1]

        return (
            np.arange(zero_crossing_idx[-1], next_to_ice_edge_idx+1)
            if sea_ice_is_left else
            np.arange(next_to_ice_edge_idx, zero_crossing_idx[0])
        )

        # x = np.arange(frb_gradient.shape[0])
        # plt.figure(dpi=150)
        # plt.scatter(x, frb_gradient)
        # plt.scatter(x[a], frb_gradient[a])
        # plt.axhline(0)
        # plt.show()
        #
        # return a

    @property
    def l2_input_vars(self):
        return ["surface_type",
                "leading_edge_width",
                "sea_ice_freeboard",
                "distance_to_ocean"]

    @property
    def l2_output_vars(self):
        return ["flag_miz"]

    @property
    def error_bit(self):
        return self.error_flag_bit_dict["filter"]


def numpy_smooth(x, window):
    return np.convolve(x, np.ones(window)/window)


def scipy_smooth(x, window):
    """ Numpy implementation of the IDL SMOOTH function """
    from scipy.ndimage.filters import uniform_filter
    return uniform_filter(x, size=window)


def astropy_smooth(x, window ,**kwargs):
    kernel = Box1DKernel(window)
    return convolve(x, kernel, **kwargs)


def interp1d_gap_filling(y: np.ndarray, **interp_kwargs) -> np.ndarray:
    """
    Fill not-finite gaps by linear interpolation using `scipy.interp1d`with bounds_error
    turned off by default. Keywords to this function are passed to interp1d.
    :param y:
    :return:
    """
    x = np.arange(y.shape[0])
    valid_idx = np.where(np.isfinite(y))[0]
    f = interp1d(x[valid_idx], y[valid_idx],
                 bounds_error=False,
                 **interp_kwargs
                 )
    return f(x)


def lowess_smooth(x: np.ndarray,
                  y: np.ndarray,
                  window: int) -> np.ndarray:
    """
    Simple lowess smoother
    :param x:
    :param y:
    :param window:
    :return:
    """
    data_fraction = min(float(window) / float(x.shape[0]), 1.0)
    return sm.nonparametric.lowess(y, x, frac=data_fraction, return_sorted=False)


def idl_smooth(x: np.ndarray,
               window: int
               ) -> np.ndarray:
    """ Implementation of the IDL smooth(x, window, /EDGE_TRUNCATE, /NAN)"""
    smoothed = np.full(x.shape, np.nan)
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
    # DEPR
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
