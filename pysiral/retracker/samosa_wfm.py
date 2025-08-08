# -*- coding: utf-8 -*-

"""
This is the SAMOSA+ retracker variant that uses `samosa-waveform-model`, a
re-implementation of the SAMPy package. The objective of the re-implementation
is better performance (code optimization and multiprocessing) and a greater flexibility
for retracking settings (sub-waveform retracking, limiting parameters of the fit)

NOTES:
======

List of fit modes:

- single fit - all parameters
- double fit - all parameters
- single fit - fixed mean square slope

All variants with options for single core and multi-processing

TODO: Implement sigma0/windspeed computation (per switch in config file)
TODO: Computation of residuals should be configurable method

Workflow

1. Colocate all input variables in dedicated list of input data classes
2. Initialize specified fit method (one class per fit method?)
3. Process waveforms with/without multiprocessing
4. Organize output

"""

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"

import multiprocessing
import time

import numpy as np
from functools import partial

from loguru import logger
from typing import Tuple, List, Dict, Optional, Any, Literal, Callable, get_args
from dataclasses import dataclass, field

from scipy.optimize import least_squares, OptimizeResult
from scipy.signal import argrelmin

from samosa_waveform_model import (
    ScenarioData, SensorParameters, SAMOSAWaveformModel, PlatformLocation, SARParameters,
    WaveformModelParameters
)
from samosa_waveform_model.dataclasses import WaveformModelOutput

from pysiral import psrlcfg
from pysiral.core.config import RadarModes
from pysiral.retracker import BaseRetracker
from pysiral.l1data import Level1bData
from pysiral.l2data import Level2Data

# NOTE:
# There are three fit methods that may be chosen for different surface types.
# E.g., for open ocean it has been standard practice to have a first fit only
# for mean square slope to estimate mean square slope and then a second fit
# with invariable mean square slope ("two_fits_mss_first_swh_second").
# This implementation also supports fitting all parameters in single fit
# ("single_fit_mss_swh") of only of significant wave height with mean square
# slope sourced from waveform parameter ("single_fit_mss_preset").
# (Epoch is always fitted).
VALID_METHOD_LITERAL = Literal["samosap_standard", "samosap_specular", "samosap_single"]
VALID_METHODS = get_args(VALID_METHOD_LITERAL)

# Default fit tolerances for the least squares optimization from SAMPy
DEFAULT_FIT_KWARGS = dict(
    loss="linear",
    method="trf",
    ftol=1e-2,
    xtol=2*1e-3,
    gtol=1e-2,
    max_nfev=100,
)

# Flag to enable debug mode for SAMOSA+ retracker
# if set to True, the SAMOSA+ waveform model will store
# all intermediate results in the fit process.
SAMOSA_WFM_COLLECT_FIT_PARAMS = False

# Coefficients for empirical nu - ocog width relation
# These coefficients have been derived from a limited
# number of orbits and are not yet validated.
# TODO: To be refined (currently based on single orbit)
NU_OCOG_COEFS = (1.11110807e+07, 2.50396017e+00)


@dataclass
class WaveformModelParametersFit(WaveformModelParameters):
    samosa_step: Literal["step1", "step2"] = None
    num_ddm_evaluations: int = -1


@dataclass
class NormedWaveform:
    """
    Data container for the input waveform
    """
    power: np.ndarray
    range_bins: np.ndarray
    tau: np.ndarray
    ocog_width: float
    scaling_factor: float
    radar_mode_flag: int
    window_delay: float
    transmit_power: float
    look_angles: np.ndarray
    surface_type: str
    surface_class: str
    thermal_noise: float = 0.0
    first_maximum_index: int = None
    idx: int = None
    absolute_maximum: float = field(init=False)

    def __post_init__(self):
        self.absolute_maximum = np.nanmax(self.power)


@dataclass
class WaveformFitData:
    """
    Dataclass that contains all data from l1b and l2 data objects
    required for the SAMOSA waveform model fit
    """
    idx: int
    waveform_data: NormedWaveform
    scenario_data: ScenarioData


@dataclass
class SAMOSAWaveformFitResult:
    """
    Dataclass for all output parameters of the waveform fitting process
    """
    epoch: float = np.nan
    retracker_range: float = np.nan
    retracker_range_standard_error: float = np.nan
    significant_wave_height: float = np.nan
    significant_wave_height_standard_error: float = np.nan
    mean_square_slope: float = np.nan
    mean_square_slope_standard_error: float = np.nan
    thermal_noise: float = np.nan
    misfit: float = np.nan
    fit_mode: str = "n/a"
    waveform: np.ndarray = None
    waveform_model: np.ndarray = None
    sub_waveform_mask: np.ndarray = None
    misfit_sub_waveform: float = None
    number_of_model_evaluations_step1: int = -1
    number_of_model_evaluations_step2: int = -1
    fit_return_status_step1: int = -2
    fit_return_status_step2: int = -2
    sigma0: float = 0.0

    @property
    def is_sub_waveform_fit(self) -> bool:
        return self.misfit_sub_waveform is not None

    @property
    def samosa_leading_edge_error(self) -> float:
        return get_samosa_leading_edge_error(self.waveform_model, self.waveform)


# TODO: Investigate Sub-Classing
class SAMOSAWaveformFit(object):
    """
    Class for fitting the SAMOSA waveform model to a single waveform.
    The purpose of this class is to minimize the CPU load for the
    fitting procedure.

    This class will therefore initiliaze the waveform model with
    all static parameters and then only re-compute the variables
    that directly depend on the fit parameters.

    Usage:

        fit_cls = SAMOSAWaveformFit(scenario_data, normed_waveform, ...)
        fit_result = scipy.optimize.least_squares(fit_cls.fit_func, ...)
    """

    def __init__(
            self,
            scenario_data: ScenarioData,
            normed_waveform: NormedWaveform,
            waveform_model: Optional[SAMOSAWaveformModel] = None,
            sub_waveform_mask: Optional[np.ndarray] = None,
            waveform_model_kwargs: Optional[Dict] = None,
            step1_fixed_nu_value: float = 0.0,
            step2_fixed_swh_value: float = 0.0,
            amplitude_is_free_param: bool = True,
            method: str = None
    ) -> None:

        # Default waveform model kwargs to empty dict
        waveform_model_kwargs = waveform_model_kwargs if isinstance(waveform_model_kwargs, dict) else {}

        # initialize waveform mode if no instance has been provided.
        self.samosa_waveform_model = (
            waveform_model if isinstance(waveform_model, SAMOSAWaveformModel) else
            SAMOSAWaveformModel(scenario_data, **waveform_model_kwargs)
        )
        self.normed_waveform = normed_waveform

        # The first fit step in the SAMOSA+ retracker uses a fixed nu value
        self.step1_fixed_nu_value = step1_fixed_nu_value
        self.step2_fixed_swh_value = step2_fixed_swh_value
        self.amplitude_is_free_param = amplitude_is_free_param
        self.nu_ocog_coefs = NU_OCOG_COEFS

        # Mask is empty by default (mask=True -> point does not generate a residual value)
        self.sub_waveform_mask = sub_waveform_mask
        self.method = method

    def fit_func(self, fit_args: List[float], *_) -> np.ndarray:
        """
        Fit function for the least square algorithm. Computes and returns
        the residual between waveform and SAMOSA+ waveform model.

        Also stores the last waveform model in class instance, which then
        can be accessed later after the optimization is complete.

        :param fit_args: Waveform model parameter

        :return: Difference between waveform and waveform model.
            be minimized by least square process.
        """

        if self.amplitude_is_free_param:
            epoch, significant_wave_height, nu, amplitude_scale = fit_args
        else:
            epoch, significant_wave_height, nu = fit_args
            amplitude_scale = 1.0

        waveform_model = get_model_from_args(
            self.samosa_waveform_model,
            [epoch * 1e-9, significant_wave_height, nu],
            thermal_noise=self.normed_waveform.thermal_noise,
            amplitude_scale=amplitude_scale
        )
        return self.compute_residuals(waveform_model)

    def fit_func_samosap_standard_step1(self, fit_args: List[float], *_) -> np.ndarray:
        """
        Fit of the first step in the SAMOSA+ fitting process.

        :param fit_args:
        :param _:

        :return: (Masked) residual vector
        """
        if self.amplitude_is_free_param:
            epoch, significant_wave_height, amplitude_scale = fit_args
        else:
            epoch, significant_wave_height = fit_args
            amplitude_scale = 1.0

        nu = self.step1_fixed_nu_value
        waveform_model = get_model_from_args(
            self.samosa_waveform_model,
            [epoch * 1e-9, significant_wave_height, nu],
            thermal_noise=self.normed_waveform.thermal_noise,
            amplitude_scale=amplitude_scale
        )
        return self.compute_residuals(waveform_model)

    def fit_func_samosap_standard_step2(self, fit_args: List[float], *_) -> np.ndarray:
        """
        Fit of the second step in the SAMOSA+ fitting process.

        :param fit_args:
        :param _:

        :return: Masked residual vector
        """
        if self.amplitude_is_free_param:
            epoch_ns, nu, amplitude_scale = fit_args
        else:
            epoch_ns, nu = fit_args
            amplitude_scale = 1.0

        significant_wave_height = self.step2_fixed_swh_value
        waveform_model = get_model_from_args(
            self.samosa_waveform_model,
            [epoch_ns * 1e-9, significant_wave_height, nu],
            thermal_noise=self.normed_waveform.thermal_noise,
            amplitude_scale=amplitude_scale
        )
        return self.compute_residuals(waveform_model)

    def compute_residuals(self, waveform_model: WaveformModelOutput) -> np.ndarray:
        """
        Compute the residuals between the waveform model and the waveform data
        depending on the sub-waveform mask.

        Note: The peak power is scaled to 1 in both the waveform and the model
              Therefore amplitude scaling is not applied and only thermal noise
              is added.

        :param waveform_model: SAMOSA+ waveform model output

        :return: residual vector
        """
        # Compute the full residual
        full_residual = waveform_model.power - self.normed_waveform.power

        # No sub-waveform mask is provided: Return full residual
        if self.sub_waveform_mask is None:
            return full_residual

        # sub-waveform mask is provided: Linear interpolation of valid range gates
        x = np.arange(full_residual.size)
        valid_entries = np.logical_not(self.sub_waveform_mask)
        return np.interp(x, x[valid_entries], full_residual[valid_entries])


class SAMOSAWaveformCollectionFit(object):
    """
    Handles the SAMOSA waveform model fitting of a waveform
    collection. Can be configured to use several fitting
    strategies with or without multiprocessing.

    The notation convention for class methods is `_fit_{fit_method}{_mp: use_multiprocessing: True}`.
    Therefore, each fit method needs to class methods: One for single process and one for parallel
    processes.

    :param fit_method: Name of the fitting method (choices: "mss_swh", "mss_preset", "mss_seperate")
    :param predictor_method:
    :param use_multiprocessing: Flag if to use multiprocessing (default: False)
    :param num_processes: Number of multiprocessing workers
    :param predictor_kwargs: Keyword arguments for the parameter predictor class
    :param least_squares_kwargs: Keyword arguments for scipy.optimize.least_squares
    """

    def __init__(
        self,
        fit_method: VALID_METHOD_LITERAL,
        predictor_method: str,
        use_multiprocessing: bool = False,
        num_processes: Optional[int] = None,
        samosap_fit_kwargs: Optional[Dict] = None,
        predictor_kwargs: Optional[Dict] = None,
        least_squares_kwargs: Optional[Dict] = None,
    ) -> None:

        self.fit_method = fit_method
        self.predictor_method = predictor_method
        self.samosap_fit_kwargs = {} if samosap_fit_kwargs is None else samosap_fit_kwargs
        self.least_squares_kwargs = DEFAULT_FIT_KWARGS if least_squares_kwargs is None else least_squares_kwargs
        self.predictor_kwargs = {} if predictor_kwargs is None else predictor_kwargs
        self.use_multiprocessing = use_multiprocessing
        self.num_processes = num_processes
        self.fit_method = self._select_fit_method()

    def fit(self, waveform_collection: List[WaveformFitData]) -> List[SAMOSAWaveformFitResult]:
        """
        Main entry method for the waveform collection fitting.

        :param waveform_collection: Waveform collection data

        :return: List of waveform fit results
        """
        return self.fit_method(waveform_collection)

    def _select_fit_method(self) -> Callable:
        """
        Selects to appropiate fit method based on class configuration.

        :raises AttributeError": Invalid configuration

        :return: callable fit method
        """
        mp_string = "_mp" if self.use_multiprocessing else ""
        fit_method_name = f"_fit_{self.fit_method}{mp_string}"
        return getattr(self, fit_method_name)

    def _fit_samosap_single(self, waveform_collection: List[WaveformFitData]) -> List[SAMOSAWaveformFitResult]:
        """
        Computes Samosa waveform model fit in the main process (no multiprocessing)

        :param waveform_collection: List of fit input

        :return: List of fit outputs
        """
        return [
            samosa_fit_samosap_single(
                fit_data,
                samosap_fit_kwargs=self.samosap_fit_kwargs,
                least_squares_kwargs=self.least_squares_kwargs,
                predictor_kwargs=self.predictor_kwargs
            )
            for fit_data in waveform_collection
        ]

    def _fit_samosap_single_mp(self, waveform_collection: List[WaveformFitData]) -> List[SAMOSAWaveformFitResult]:
        """
        Computes Samosa waveform model fit with multiprocessing

        :param waveform_collection: List of fit input

        :return: List of fit outputs
        """
        return self._mp_fit(samosa_fit_samosap_single, waveform_collection)

    def _fit_samosap_standard(self, waveform_collection: List[WaveformFitData]) -> List[SAMOSAWaveformFitResult]:
        """
        Computes Samosa waveform model fit in the main process (no multiprocessing)

        :param waveform_collection: List of fit input

        :return: List of fit outputs
        """
        return [
            samosa_fit_samosap_standard(
                fit_data,
                samosap_fit_kwargs=self.samosap_fit_kwargs,
                least_squares_kwargs=self.least_squares_kwargs,
                predictor_kwargs=self.predictor_kwargs
            )
            for fit_data in waveform_collection
        ]

    def _fit_samosap_standard_mp(self, waveform_collection: List[WaveformFitData]) -> List[SAMOSAWaveformFitResult]:
        """
        Computes Samosa waveform model fit with multiprocessing

        :param waveform_collection: List of fit input

        :return: List of fit outputs
        """
        return self._mp_fit(samosa_fit_samosap_standard, waveform_collection)

    def _fit_samosap_specular(self, waveform_collection: List[WaveformFitData]) -> List[SAMOSAWaveformFitResult]:
        """
        Computes Samosa waveform model fit in the main process (no multiprocessing)

        :param waveform_collection: List of fit input

        :return: List of fit outputs
        """
        return [
            samosa_fit_samosap_specular(
                fit_data,
                samosap_fit_kwargs=self.samosap_fit_kwargs,
                least_squares_kwargs=self.least_squares_kwargs,
                predictor_kwargs=self.predictor_kwargs
            )
            for fit_data in waveform_collection
        ]

    def _fit_samosap_specular_mp(self, waveform_collection: List[WaveformFitData]) -> List[SAMOSAWaveformFitResult]:
        """
        Computes Samosa waveform model fit with multiprocessing

        :param waveform_collection: List of fit input

        :return: List of fit outputs
        """
        return self._mp_fit(samosa_fit_samosap_specular, waveform_collection)

    def _mp_fit(self, func, waveform_collection):
        pool = multiprocessing.Pool(self.num_processes)
        fit_func = partial(
            func,
            samosap_fit_kwargs=self.samosap_fit_kwargs,
            least_squares_kwargs=self.least_squares_kwargs,
            predictor_kwargs=self.predictor_kwargs,
        )
        fit_results = pool.map(fit_func, waveform_collection)

        pool.close()
        pool.join()
        return fit_results


class SAMOSAModelParameterPrediction(object):
    """
    Class to get initial guess and bounds of SAMOSA waveform fit parameter

    TODO: Class Docstr
    """

    def __init__(
        self,
        method: VALID_METHOD_LITERAL,
        surface_type: str,
        initial_guess: Dict[str, Tuple[float, float]],
        bounds: Dict[str, Tuple[float, float]],
        amplitude_is_free_parameter: Tuple[bool, bool] = (True, True),
    ):

        # Input Arguments
        self.surface_type = surface_type
        self.amplitude_is_free_parameter = amplitude_is_free_parameter
        self.method = method
        self.bounds = bounds
        self.initial_guess = initial_guess

        # Method specific functions (as properties)
        self.first_guess_method = getattr(self, f"_get_first_guess_{method}")
        self.bounds_method = getattr(self, f"_get_bounds_{method}")

    def get(self, waveform: NormedWaveform, **kwargs) -> Tuple[Tuple, Tuple, Tuple]:
        """
        Get first guess and fit bounds for fit parameters from input waveform. The specific
        parameters and their order in the return tuples is defined by the property `method`
        of this class.

        :param waveform: Waveform data
        :param kwargs: Potential keyword arguments, e.g. fit mode

        :return: Tuples of first guess, lower bounds and upper bounds for the fit parametes.
        """
        first_guess = self.first_guess_method(waveform, **kwargs)
        lower_bounds, upper_bounds = self.bounds_method(waveform.tau, first_guess, **kwargs)
        return first_guess, lower_bounds, upper_bounds

    def _get_first_guess_samosap_specular(self, waveform: NormedWaveform, **_) -> Tuple:
        """
        Estimate the first guess for the SAMOSA+ specular waveform fitting mode (fixed swh, only nu).
        Equivalent to SAMOSA+ standard fit method with mode 2.

        :param waveform: Waveform data

        :return: Parameter first guess (epoch, nu, amplitude)
        """
        return self._get_first_guess_samosap_standard(waveform, mode=2)

    def _get_first_guess_samosap_single(self, waveform: NormedWaveform, **_) -> Tuple:
        """
        Estimate the first guess for the SAMOSA+ specular waveform fitting mode (fixed swh, only nu).
        Equivalent to SAMOSA+ standard fit method with mode 2.

        :param waveform: Waveform data

        :return: Parameter first guess (epoch, nu, amplitude)
        """
        epoch_first_guess = waveform.tau[waveform.first_maximum_index]
        nu_first_guess = get_nu_from_ocog_width(waveform.ocog_width, NU_OCOG_COEFS)
        swh_first_guess = self.initial_guess["swh"]
        return_tuple = epoch_first_guess * 1e9, swh_first_guess, nu_first_guess, self.initial_guess["amplitude"]
        return return_tuple if self.amplitude_is_free_parameter else return_tuple[:-1]

    def _get_bounds_samosap_single(self, tau, first_guess) -> Tuple[Tuple, Tuple]:
        """

        :return: Parameter fit bounds: (lower bounds, upperbounds) [mode 1: epoch, swh, amplitude; mode 2: epoch, nu]
        """
        range_gates_after_fmi = self.bounds["epoch"]["range_gates_after_fmi"]
        epoch_bounds = get_epoch_bounds(tau, first_guess[0], range_gates_after_fmi=range_gates_after_fmi)
        lb = epoch_bounds[0] * 1e9, self.bounds["swh"][0], self.bounds["nu"][0], self.bounds["amplitude"][0]
        ub = epoch_bounds[1] * 1e9, self.bounds["swh"][1], self.bounds["nu"][1], self.bounds["amplitude"][0]
        return (lb, ub) if self.amplitude_is_free_parameter else (lb[:-1], ub[:-1])

    def _get_bounds_samosap_specular(self, tau, first_guess, **_) -> Tuple[Tuple, Tuple]:
        """
        Estimate the fit bounds for the SAMOSA+ specular waveform fitting mode (fixed swh, only nu).
        Equivalent to SAMOSA+ standard fit method with mode 2.

        :param tau: range gates in seconds
        :param first_guess: first guess triplet

        :return: Parameter first guess (epoch, nu, amplitude)
        """
        return self._get_bounds_samosap_standard(tau, first_guess, mode=2)

    def _get_first_guess_samosap_standard(self, waveform: NormedWaveform, mode: int = -1) -> Tuple:
        """
        Estimate the first guess for the two-step SAMOSA+ waveform fitting approch.
        The fit parameter depend on the specific mode/step.

        First fit step (mode=1): [epoch, significant_waveheight, amplitude scale]
        Seconds fit step (mode=2): [epoch, mean square slope]

        :param waveform: Waveform data
        :param mode: The SAMOSA+ waveform model mode (conf.STEP in SAMPy)

        :raises ValueError: mode not 1 or 2

        :return: Parameter first guess [mode 1: epoch, swh, amplitude; mode 2: epoch, nu]
        """
        epoch_first_guess = waveform.tau[waveform.first_maximum_index]

        # Fitting epoch, swh, [amplitude]
        if mode == 1:
            fit_params = epoch_first_guess * 1e9, self.initial_guess["swh"], self.initial_guess["amplitude"]
            return fit_params if self.amplitude_is_free_parameter[0] else fit_params[:-1]

        # Fitting epoch, nu, [amplitude]
        elif mode == 2:
            fit_params = epoch_first_guess * 1e9, self.initial_guess["nu"], self.initial_guess["amplitude"]
            return fit_params if self.amplitude_is_free_parameter[1] else fit_params[:-1]
        else:
            raise ValueError(f"mode={mode} not in [1, 2]")

    def _get_bounds_samosap_standard(self, tau, first_guess, mode: int = -1) -> Tuple[Any, Any]:
        """
        Estimate the parameter bounds for the two-step SAMOSA+ waveform fitting approch.
        The bounds depend on the specific mode/step.

        First fit step (mode=1): [epoch, significant_waveheight, amplitude scale]
        Seconds fit step (mode=2): [epoch, mean square slope, amplitude scale]

        The bounds for the epoch are set to end behind the first maximum. For the other
        parameter bounds need to be provided during the initialization of this class.

        :param tau: range gates in seconds
        :param first_guess: first guess triplet
        :param mode: The SAMOSA+ waveform model mode (conf.STEP in SAMPy)

        :raises ValueError: mode not 1 or 2

        :return: Parameter fit bounds: (lower bounds, upperbounds) [mode 1: epoch, swh, amplitude; mode 2: epoch, nu]
        """
        range_gates_after_fmi = self.bounds["epoch"]["range_gates_after_fmi"]
        epoch_bounds = get_epoch_bounds(tau, first_guess[0], range_gates_after_fmi=range_gates_after_fmi)
        amp_bounds = self.bounds["amplitude"]

        # Fitting epoch, swh, [amplitude]
        if mode == 1:
            swh_bounds = self.bounds["swh"]
            return self._compile_bounds(epoch_bounds, swh_bounds, amp_bounds, mode)

        # Fitting epoch, nu, [amplitude]
        elif mode == 2:
            nu_bounds = self.bounds["nu"]
            return self._compile_bounds(epoch_bounds, nu_bounds, amp_bounds, mode)
        else:
            raise ValueError(f"mode={mode} not in [1, 2]")

    def _compile_bounds(self, epoch_bnds, param_bnds, amp_bnds, mode) -> Tuple[Any, Any]:
        lb = epoch_bnds[0] * 1e9, param_bnds[0], amp_bnds[0]
        ub = epoch_bnds[1] * 1e9, param_bnds[1], amp_bnds[1]
        return (lb, ub) if self.amplitude_is_free_parameter[mode-1] else (lb[:-1], ub[:-1])


class SAMOSAPlusRetracker(BaseRetracker):

    def __init__(self) -> None:
        super(SAMOSAPlusRetracker, self).__init__()
        self._retracker_params = {}

    def l2_retrack(
            self,
            rng: np.ndarray,
            wfm: np.ndarray,
            indices: np.ndarray,
            radar_mode: np.ndarray,
            is_valid: np.ndarray
    ) -> None:
        """
         API method for the retracker interface in the Level-2 processor.

        :param rng: The range per waveform samples for all waveforms in the Level-1 data object
        :param wfm: All waveforms in the Level-1 data object
        :param indices: List of waveforms for the retracker
        :param radar_mode: Flag indicating the radar mode for all waveforms
        :param is_valid: Error flag for all waveforms

        :return: None (Output is added to the instance)
        """

        # Run the retracker
        # NOTE: Output is a SAMOSAFitResult dataclass for each index in indices.
        _, kwargs = self._get_waveform_model_fit_configuration()
        t0 = time.time()
        fit_results = self._samosa_plus_retracker(rng, wfm, indices, radar_mode)
        t1 = time.time()
        duration = t1 - t0
        if kwargs.get("use_multiprocessing", False):
            n_processes = kwargs.get("num_processes", 1)
        else:
            n_processes = 1
        secs_per_waveform = n_processes * duration / len(indices)
        logger.debug(f"Waveform fitting took {duration:.2f} sec with {n_processes=} -> ({secs_per_waveform=:.2f})")

        # Store retracker properties (including range)
        self._store_retracker_properties(fit_results, indices)

        # Set/compute uncertainty
        self._set_range_uncertainty(fit_results, indices)

        # Add range biases (when set in config file)
        # self._set_range_bias(radar_mode)

        # Add filter
        self._apply_filter()

        # Add retracker parameters to the l2 data object
        self._l2_register_retracker_parameters()

    def create_retracker_properties(self, n_records: int) -> None:
        """
        Initialize retracker properties with correct arrays shapes (shape = (n_records, )).
        The list of parameters depends on whether the SAMOSA_DEBUG_MODE flag is set.

        NOTE: The properties are set to an array, but can be accessed as `self.{property_name}`
        via the __getattr__ method.

        :param n_records: Number of Level-2 data points
        """

        parameter = [
            ("misfit", np.nan, np.float32),
            ("misfit_sub_waveform", np.nan, np.float32),
            ("leading_edge_error", np.nan, np.float32),
            ("swh", np.nan, np.float32),
            ("mean_square_slope", np.nan, np.float32),
            ("wind_speed", np.nan, np.float32),
            ("epoch", np.nan, np.float32),
            ("guess", np.nan, np.float32),
            ("fit_num_func_eval_step1", -1, np.int32),
            ("fit_num_func_eval_step2", -1, np.int32),
            ("fit_return_status_step1", -2, np.int32),  # least squares valid range -1 to 4
            ("fit_return_status_step2", -2, np.int32),  # least squares valid range -1 to 4
            ("Pu", np.nan, np.float32),
            ("rval", np.nan, np.float32),
            ("kval", np.nan, np.float32),
            ("pval", np.nan, np.float32),
            ("cval", np.nan, np.float32),
        ]
        for parameter_name, fill_value, dtype in parameter:
            self._retracker_params[parameter_name] = np.full(n_records, fill_value, dtype=dtype)

        parameter = [
            "waveform",
            "waveform_model",
            "sub_waveform_mask"
        ]
        n_range_gates = self._l1b.waveform.power.shape[1]
        for parameter_name in parameter:
            self._retracker_params[parameter_name] = np.full((n_records, n_range_gates), np.nan, dtype=np.float32)

    def _samosa_plus_retracker(
            self,
            rng: np.ndarray,
            wfm: np.ndarray,
            indices: np.ndarray,
            radar_mode: np.ndarray
    ) -> List[SAMOSAWaveformFitResult]:
        """
        Performs the retracking via SAMOSA waveform model fits for
        all target waveforms.

        :param rng: Range window arrays
        :param wfm: Waveform arrays
        :param indices: Target waveform indices
        :param radar_mode: Radar modes

        :return: List of SAMOSA plus waveform fit results
        """

        # Aggregate all necessary input data for the SAMOSA waveform model fit
        waveform_collection = self._get_waveform_fit_data(rng, wfm, indices, radar_mode)

        # Configure and execute waveform fitting for all target waveforms
        args, kwargs = self._get_waveform_model_fit_configuration()
        waveform_fits = SAMOSAWaveformCollectionFit(*args, **kwargs)
        return waveform_fits.fit(waveform_collection)

    def _get_waveform_fit_data(
            self,
            rng: np.ndarray,
            wfm: np.ndarray,
            indices: np.ndarray,
            radar_mode: np.ndarray
    ) -> List[WaveformFitData]:
        """
        Collects all the data required for the waveform model fits.
        Whenever possible,

        :param rng:
        :param wfm:
        :param indices:
        :param radar_mode:
        :return:
        """

        waveform_collection = []
        for idx in indices:
            # Create the SAMOSA Waveform scenario data
            look_angles = get_look_angles(self._l1b, idx)
            scenario_data = self._get_scenario_data(self._l1b, self._l2, idx, look_angles)
            waveform_data = self._get_normed_waveform(
                rng[idx, :],
                wfm[idx, :],
                scenario_data.rp.tau,
                radar_mode[idx],
                self._l1b.classifier.window_delay[idx],
                self._l1b.classifier.transmit_power[idx],
                look_angles,
                self._l1b.classifier.first_maximum_index[idx],
                self._l1b.classifier.ocog_width[idx],
                idx=idx
            )
            waveform_collection.append(WaveformFitData(idx, waveform_data, scenario_data))

        return waveform_collection

    def _get_scenario_data(
            self,
            l1: Level1bData,
            l2: Level2Data,
            idx: int,
            look_angles
    ) -> ScenarioData:
        """
        Creates the scenario data for the SAMOSA waveform model. This dataclass
        is needed for the initialization of the SAMOSA waveform model class and
        contains information on the radar mode and radar processing settings as
        well as the geometry (location, altitude & altitude rate, attitude and
        velocity of the altimeter platform).

        The specific radar parameters are expected to be included in the
        `samosa_waveform_model` and are retrieved by specifying platform
        and radar mode.

        Geometry parameters are available in the Level-1 and Level-2 data
        containers.

        :param l1: Level-1 data container
        :param l2: Level-2 data container
        :param idx: Target Waveform index
        :param look_angles:

        :return: SAMOSA waveform model scenario data for specific waveform
        """

        radar_mode_name = RadarModes.get_name(l1.waveform.radar_mode[idx])
        platform = l2.info.mission

        sp = SensorParameters.get(platform, radar_mode_name)

        # pysiral specific: All waveforms windowed to 256 range gates
        # Here to be changed without zero-padding factor (which is same for SAR and SARin).
        if radar_mode_name == "sin":
            sp.range_gates_per_pulse = 128

        location_data = dict(
            latitude=l2.latitude[idx],
            longitude=l2.longitude[idx],
            altitude=l2.altitude[idx],
            height_rate=l1.time_orbit.altitude_rate[idx],
            pitch=np.radians(l1.time_orbit.antenna_pitch[idx]),
            roll=np.radians(l1.time_orbit.antenna_roll[idx]),
            velocity=total_velocity_from_vector(
                self._l1b.classifier.satellite_velocity_x[idx],
                self._l1b.classifier.satellite_velocity_y[idx],
                self._l1b.classifier.satellite_velocity_z[idx]
            )
        )
        geo = PlatformLocation(**location_data)
        sar = SARParameters(look_angles=look_angles)
        sar.compute_multi_look_parameters(geo=geo, sp=sp)

        return ScenarioData(sp, geo, sar)

    def _get_normed_waveform(
            self,
            rng: np.ndarray,
            wfm: np.ndarray,
            tau: np.ndarray,
            radar_mode: int,
            window_delay: float,
            transmit_power: float,
            look_angles: np.ndarray,
            first_maximum_index: int,
            ocog_width: float,
            idx: int = None
    ) -> NormedWaveform:
        """
        Return a normed representation of the input waveform

        :param rng:
        :param wfm:
        :param radar_mode:

        :return:
        """
        scaling_factor = float(1.0 / wfm[first_maximum_index])
        normed_waveform = wfm * scaling_factor

        return NormedWaveform(
            normed_waveform,
            rng,
            tau,
            ocog_width,
            scaling_factor,
            radar_mode,
            window_delay,
            transmit_power,
            look_angles,
            self._options.get("surface_type"),
            self._options.get("surface_class"),
            first_maximum_index=first_maximum_index,
            idx=idx
        )

    def _get_waveform_model_fit_configuration(self) -> Tuple[List, Dict]:
        """
        Constructs the arguments and keywords for the SAMOSAWaveformCollectionFit class,
        defining the fit method and its configration

        :return: arguments and keyword arguments for the waveform model
        """

        fit_method = self._options.get("fit_method")
        if fit_method is None:
            raise ValueError("Mandatory configuration parameter `fit_method not specified")

        samosap_fit_kwargs = self._options.get("samosap_fit_kwargs")
        if samosap_fit_kwargs is None:
            raise ValueError("Mandatory configuration parameter `samosap_fit_kwargs` not specified")

        predictor_method = self._options.get("predictor_method")
        if fit_method is None:
            raise ValueError("Mandatory configuration parameter `predictor_method` not specified")
        args = [fit_method, predictor_method]

        num_processes = self._options.get("num_processes")
        num_processes = psrlcfg.CPU_COUNT if num_processes == "pysiral-cfg" else int(num_processes)
        kwargs = {
            "use_multiprocessing": self._options.get("use_multiprocessing", False),
            "num_processes": num_processes,
            "samosap_fit_kwargs": samosap_fit_kwargs,
            "least_squares_kwargs": self._options.get("least_squares_kwargs", DEFAULT_FIT_KWARGS),
            "predictor_kwargs": self._options.get("predictor_kwargs", {}),
        }
        return args, kwargs

    def _set_range_uncertainty(
            self,
            fit_results: List[SAMOSAWaveformFitResult],
            indices: np.ndarray
    ) -> None:
        """
        The retracker range uncertainty is computed directly from the
        epoch standard error (from the least squares cost function)
        in the waveform model fitting process.

        :param fit_results: A list of waveform model fit results
        :param indices: List of waveform indices
        """
        for index, fit_result in zip(indices, fit_results):
            self._uncertainty[index] = fit_result.retracker_range_standard_error

    def _apply_filter(self) -> None:
        """
        Apply a filter to the retracker results. E.g., range values exceeding
        a certain misfit value or leading edge error are masked.

        :return: None (Variables are changed in place
        """

        filter_args = self._options.get("misfit_filter")
        if not filter_args:
            return

        # Compute filter flag indicated waveforms with misfit/leading edge error values exceeding threshold
        invalid_misfit = self.misfit_sub_waveform > filter_args["max_misfit"]
        invalid_leading_edge_fit = self.leading_edge_error > filter_args["max_leading_edge_error"]
        radar_freeboard_filter_flag = np.logical_or(invalid_misfit, invalid_leading_edge_fit)

        # Apply filter to range and uncertainty
        logger.info(f"Filtering {np.sum(radar_freeboard_filter_flag)} waveforms due to misfit & leading edge error")
        self._range[radar_freeboard_filter_flag] = np.nan
        self._uncertainty[radar_freeboard_filter_flag] = np.nan

    def _store_retracker_properties(
        self,
        fit_results: List[SAMOSAWaveformFitResult],
        indices: np.ndarray
    ) -> None:
        """
        Store the output of the SAMOSA+ retracker in the class and maps
        the fit result to the corresponding indices in the Level-2 data.

        :param fit_results:
        :param indices:

        :return: None
        """

        for index, fit_result in zip(indices, fit_results):

            # Main retracker properties
            self._range[index] = fit_result.retracker_range
            self._power[index] = fit_result.sigma0

            # Derived physical parameters
            # self.wind_speed[index] = func_wind_speed([fit_result.sigma0])

            # Waveform model fit parameters
            self.swh[index] = fit_result.significant_wave_height
            self.misfit[index] = fit_result.misfit
            self.misfit_sub_waveform[index] = fit_result.misfit_sub_waveform
            self.mean_square_slope[index] = fit_result.mean_square_slope
            self.epoch[index] = fit_result.epoch
            self.leading_edge_error[index] = fit_result.samosa_leading_edge_error

            # Store waveform and waveform model for debugging
            self.waveform[index, :] = fit_result.waveform
            self.waveform_model[index, :] = fit_result.waveform_model
            if fit_result.sub_waveform_mask is not None:
                self.sub_waveform_mask[index, :] = fit_result.sub_waveform_mask

            # Fit method statistics
            self.fit_num_func_eval_step1[index] = fit_result.number_of_model_evaluations_step1
            self.fit_num_func_eval_step2[index] = fit_result.number_of_model_evaluations_step2
            self.fit_return_status_step1[index] = fit_result.fit_return_status_step1
            self.fit_return_status_step2[index] = fit_result.fit_return_status_step2

    def _l2_register_retracker_parameters(self) -> None:
        """
        Add retracker variables to the L2 data object. The specific variables
        depend on surface type.
        """

        # General auxiliary variables
        self.register_auxdata_output("sammf", "samosa_misfit", self.misfit)
        self.register_auxdata_output("samswh", "samosa_swh", self.swh)
        self.register_auxdata_output("sammfsw", "samosa_misfit_sub_waveform", self.misfit_sub_waveform)
        self.register_auxdata_output("samlee", "samosa_leading_edge_error", self.leading_edge_error)
        self.register_auxdata_output("sammss", "samosa_mean_square_slope", self.mean_square_slope)
        self.register_auxdata_output("samfnfe1", "samosa_fit_num_func_eval_step1", self.fit_num_func_eval_step1)
        self.register_auxdata_output("samfnfe2", "samosa_fit_num_func_eval_step2", self.fit_num_func_eval_step2)
        self.register_auxdata_output("samfrs", "samosa_fit_return_status_step1", self.fit_return_status_step1)
        self.register_auxdata_output("samfrs", "samosa_fit_return_status_step2", self.fit_return_status_step2)

        # Waveform and waveform model
        # Register results as auxiliary data variable
        num_range_gates = self.waveform.shape[1]
        dims = {"new_dims": (("range_gates", num_range_gates),),
                "dimensions": ("time", "range_gates"),
                "add_dims": (("range_gates", np.arange(num_range_gates)),)}
        self._l2.set_multidim_auxiliary_parameter("samiwfm", "waveform", self.waveform, dims, update=True)
        self._l2.set_multidim_auxiliary_parameter("samwfm", "waveform_model", self.waveform_model, dims, update=True)
        self._l2.set_multidim_auxiliary_parameter("samswfm", "sub_waveform_mask", self.sub_waveform_mask, dims, update=True)

        # Lead and open ocean surfaces
        surface_class = self._options.get("surface_class", "undefined")
        if surface_class == "polar_ocean":
            self.register_auxdata_output("samwsp", "samosa_wind_speed", self.wind_speed)

        # Sea ice surface types
        elif surface_class == "sea_ice":
            self.register_auxdata_output("samshsd", "samosa_surface_height_standard_deviation", self.swh / 4.)

        else:
            logger.warning("No specific surface type set for SAMOSA+: No auxiliary variables added to L2")

    def __getattr__(self, item: str) -> Any:
        """
        Direct attribute access to the retracker properties dictionary

        :param item: parameter name

        :raise AttributeError: item not in self._retracker_params (see self.create_retracker_properties)

        :return:
        """
        if item in self._retracker_params:
            return self._retracker_params[item]
        else:
            raise AttributeError(f"{self.__class__.__name__} has no attribute {item}")


def samosa_fit_samosap_single(
        fit_data: WaveformFitData,
        samosap_fit_kwargs: Dict = None,
        predictor_kwargs: Dict = None,
        least_squares_kwargs: Dict = None,
        filter_trailing_edge_kwargs: Dict = None
) -> SAMOSAWaveformFitResult:
    """
    Fits the SAMOSA waveform model with all free parameters (epoch, swh, mss, amplitude) using
    the two-step standard fitting approach used for open ocean waveforms with SAMOSA+.

    This fit is intended for waveforms classified as sea ice.

    :param fit_data: Input parameters for waveform fitting process. Mainly
        contains waveform model scenario data and waveform data.
    :param samosap_fit_kwargs:
    :param predictor_kwargs: Input parameter for parameter first guess and fit bounds
    :param least_squares_kwargs: Keyword arguments to `scipy.optimize.least_squares`
    :param filter_trailing_edge_kwargs: Keyword arguments to

    :return: SAMOSA+ waveform model fit result
    """

    # Input validation
    predictor_kwargs = {} if predictor_kwargs is None else predictor_kwargs
    least_squares_kwargs = {} if least_squares_kwargs is None else least_squares_kwargs
    filter_trailing_edge_kwargs = {} if filter_trailing_edge_kwargs is None else filter_trailing_edge_kwargs

    # Unpack for readability
    scenario_data, waveform_data = fit_data.scenario_data, fit_data.waveform_data
    waveform_data.thermal_noise = compute_thermal_noise(waveform_data.power)

    # Get the sub-waveform mask
    # (unless explicitly disabled by `trailing_edge_sub_waveform_filter=False` config file)
    trailing_edge_sub_waveform_filter = samosap_fit_kwargs.get("trailing_edge_sub_waveform_filter", True)
    if trailing_edge_sub_waveform_filter:
        sub_waveform_mask = get_sub_waveform_mask(waveform_data, filter_trailing_edge_kwargs)
    else:
        sub_waveform_mask = None

    # Get first guess of fit parameters and fit bounds
    predictor = SAMOSAModelParameterPrediction("samosap_single", waveform_data.surface_type, **predictor_kwargs)

    # --- Single SAMOSA+ Fit Step 2 ---
    # This step fits a waveform with fixed swh and variable nu
    # This fit will be used
    model_parameters, fitted_model, optimize_result = samosa_fit_samosap_single_step2(
        waveform_data, scenario_data, predictor,
        least_squares_kwargs=least_squares_kwargs,
        sub_waveform_mask=sub_waveform_mask,
        **samosap_fit_kwargs
    )

    # Compute the misfit from residuals in SAMPy fashion
    misfit_subwaveform = sampy_misfit(optimize_result.fun, waveform_scale=waveform_data.absolute_maximum)
    misfit = sampy_misfit(fitted_model.power - waveform_data.power, waveform_scale=waveform_data.absolute_maximum)

    # Convert epoch to range (excluding range corrections)
    retracker_range = epoch2range(model_parameters.epoch, fit_data.waveform_data.range_bins)
    retracker_range_standard_error = 0.5 * 299792458. * model_parameters.epoch_sdev

    return SAMOSAWaveformFitResult(
         epoch=model_parameters.epoch,
         retracker_range=retracker_range,
         retracker_range_standard_error=retracker_range_standard_error,
         significant_wave_height=model_parameters.significant_wave_height,
         significant_wave_height_standard_error=model_parameters.significant_wave_height_sdev,
         mean_square_slope=model_parameters.mean_square_slope,
         mean_square_slope_standard_error=1. / model_parameters.nu_sdev,
         thermal_noise=model_parameters.thermal_noise,
         misfit=misfit,
         misfit_sub_waveform=misfit_subwaveform,
         sub_waveform_mask=sub_waveform_mask,
         fit_mode="samosap_single",
         waveform=waveform_data.power,
         waveform_model=fitted_model.power,
         number_of_model_evaluations_step2=model_parameters.num_ddm_evaluations,
         fit_return_status_step2=optimize_result.status
    )


def samosa_fit_samosap_standard(
        fit_data: WaveformFitData,
        samosap_fit_kwargs: Dict = None,
        predictor_kwargs: Dict = None,
        least_squares_kwargs: Dict = None,
        filter_trailing_edge_kwargs: Dict = None
) -> SAMOSAWaveformFitResult:
    """
    Fits the SAMOSA waveform model with all free parameters (epoch, swh, mss, amplitude) using
    the two-step standard fitting approach used for open ocean waveforms with SAMOSA+.

    This fit is intended for waveforms classified as sea ice.

    :param fit_data: Input parameters for waveform fitting process. Mainly
        contains waveform model scenario data and waveform data.
    :param samosap_fit_kwargs: Input parameters for SAMOSA+ fitting process
    :param predictor_kwargs: Input parameter for parameter first guess and fit bounds
    :param least_squares_kwargs: Keyword arguments to `scipy.optimize.least_squares`
    :param filter_trailing_edge_kwargs: Keyword arguments to

    :return: SAMOSA+ waveform model fit result
    """

    # Input validation
    samosap_fit_kwargs = {} if samosap_fit_kwargs is None else samosap_fit_kwargs
    predictor_kwargs = {} if predictor_kwargs is None else predictor_kwargs
    least_squares_kwargs = {} if least_squares_kwargs is None else least_squares_kwargs
    filter_trailing_edge_kwargs = {} if filter_trailing_edge_kwargs is None else filter_trailing_edge_kwargs

    # Unpack for readability
    scenario_data, waveform_data = fit_data.scenario_data, fit_data.waveform_data
    waveform_data.thermal_noise = compute_thermal_noise(waveform_data.power)

    # Get the sub-waveform mask
    # (unless explicitly disabled by `trailing_edge_sub_waveform_filter=False` config file)
    trailing_edge_sub_waveform_filter = samosap_fit_kwargs.get("trailing_edge_sub_waveform_filter", True)
    if trailing_edge_sub_waveform_filter:
        sub_waveform_mask = get_sub_waveform_mask(waveform_data, filter_trailing_edge_kwargs)
    else:
        sub_waveform_mask = None

    # Get first guess of fit parameters and fit bounds
    predictor = SAMOSAModelParameterPrediction("samosap_standard", waveform_data.surface_type, **predictor_kwargs)

    # --- SAMOSA+ Fit Step 1 ---
    # This step fits a waveform with fixed nu and variable swh
    model_parameters_step1, fitted_model_step1, optimize_result_step1 = samosa_fit_samosap_standard_step1(
        waveform_data, scenario_data, predictor,
        least_squares_kwargs=least_squares_kwargs,
        sub_waveform_mask=sub_waveform_mask,
        step1_fixed_nu_value=samosap_fit_kwargs["step1_fixed_nu_value"],
        amplitude_is_free_param=samosap_fit_kwargs["amplitude_is_free_param"][0]
    )

    # --- SAMOSA+ Fit Step 2 ---
    # This step fits a waveform with fixed swh and variable nu
    # This fit will be used
    model_parameters_step2, fitted_model_step2, optimize_result_step2 = samosa_fit_samosap_standard_step2(
        waveform_data, scenario_data, predictor,
        least_squares_kwargs=least_squares_kwargs,
        sub_waveform_mask=sub_waveform_mask,
        step2_fixed_swh_value=samosap_fit_kwargs["step2_fixed_swh_value"],
        amplitude_is_free_param=samosap_fit_kwargs["amplitude_is_free_param"][1]
    )

    # --- Summarize the result from two fits ---
    # Note that the waveform model from the combination of swh and nu will provide the best fit

    # Compute the misfit from residuals in SAMPy fashion
    misfit_subwaveform = sampy_misfit(optimize_result_step2.fun, waveform_scale=waveform_data.absolute_maximum)
    misfit = sampy_misfit(fitted_model_step2.power - waveform_data.power, waveform_scale=waveform_data.absolute_maximum)

    # Convert epoch to range (excluding range corrections)
    retracker_range = epoch2range(model_parameters_step2.epoch, fit_data.waveform_data.range_bins)
    retracker_range_standard_error = 0.5 * 299792458. * model_parameters_step2.epoch_sdev

    return SAMOSAWaveformFitResult(
         epoch=model_parameters_step2.epoch,
         retracker_range=retracker_range,
         retracker_range_standard_error=retracker_range_standard_error,
         significant_wave_height=model_parameters_step1.significant_wave_height,
         significant_wave_height_standard_error=model_parameters_step1.significant_wave_height_sdev,
         mean_square_slope=model_parameters_step2.mean_square_slope,
         mean_square_slope_standard_error=1. / model_parameters_step2.nu_sdev,
         thermal_noise=model_parameters_step2.thermal_noise,
         misfit=misfit,
         misfit_sub_waveform=misfit_subwaveform,
         sub_waveform_mask=sub_waveform_mask,
         fit_mode="samosap_standard",
         waveform=waveform_data.power,
         waveform_model=fitted_model_step2.power,
         number_of_model_evaluations_step1=model_parameters_step1.num_ddm_evaluations,
         number_of_model_evaluations_step2=model_parameters_step2.num_ddm_evaluations,
         fit_return_status_step1=optimize_result_step1.status,
         fit_return_status_step2=optimize_result_step2.status
    )


def samosa_fit_samosap_specular(
        fit_data: WaveformFitData,
        samosap_fit_kwargs: Dict = None,
        predictor_kwargs: Dict = None,
        least_squares_kwargs: Dict = None
) -> SAMOSAWaveformFitResult:
    """
    Fits the SAMOSA waveform model with (epoch, mss, [amplitude]) and a fixed swh
    in a single step that is equivalent of the fit step 2 in the default fit
    procedure

    This fit is intended for waveforms classified as lead waveforms.

    :param fit_data: Input parameters for waveform fitting process. Mainly
        contains waveform model scenario data and waveform data.
    :param samosap_fit_kwargs: Input parameters for SAMOSA+ fitting process
    :param predictor_kwargs: Input parameter for parameter first guess and fit bounds
    :param least_squares_kwargs: Keyword arguments to `scipy.optimize.least_squares`

    :return: SAMOSA+ waveform model fit result
    """

    # Input validation

    samosap_fit_kwargs = {} if samosap_fit_kwargs is None else samosap_fit_kwargs
    predictor_kwargs = {} if predictor_kwargs is None else predictor_kwargs
    least_squares_kwargs = {} if least_squares_kwargs is None else least_squares_kwargs

    # Unpack for readability
    scenario_data, waveform_data = fit_data.scenario_data, fit_data.waveform_data
    waveform_data.thermal_noise = compute_thermal_noise(waveform_data.power)

    # Get first guess of fit parameters and fit bounds
    predictor = SAMOSAModelParameterPrediction("samosap_specular", waveform_data.surface_type, **predictor_kwargs)

    # --- SAMOSA+ Fit Step 2 ---
    # This step fits a waveform with fixed swh and variable nu
    # This fit will be used
    model_parameters_step2, fitted_model_step2, optimize_result_step2 = samosa_fit_samosap_standard_step2(
        waveform_data, scenario_data, predictor,
        least_squares_kwargs=least_squares_kwargs,
        step2_fixed_swh_value=samosap_fit_kwargs["step2_fixed_swh_value"],
        amplitude_is_free_param=samosap_fit_kwargs["amplitude_is_free_param"][1]
    )

    # --- Summarize the result from two fits ---
    # Compute the misfit from residuals in SAMPy fashion
    misfit = sampy_misfit(fitted_model_step2.power - waveform_data.power, waveform_scale=waveform_data.absolute_maximum)

    # Convert epoch to range (excluding range corrections)
    retracker_range = epoch2range(model_parameters_step2.epoch, fit_data.waveform_data.range_bins)
    retracker_range_standard_error = 0.5 * 299792458. * model_parameters_step2.epoch_sdev

    return SAMOSAWaveformFitResult(
         epoch=model_parameters_step2.epoch,
         retracker_range=retracker_range,
         retracker_range_standard_error=retracker_range_standard_error,
         significant_wave_height=0.0,
         significant_wave_height_standard_error=0.0,
         mean_square_slope=model_parameters_step2.mean_square_slope,
         mean_square_slope_standard_error=1. / model_parameters_step2.nu_sdev,
         thermal_noise=model_parameters_step2.thermal_noise,
         misfit=misfit,
         fit_mode="samosap_specular",
         waveform=waveform_data.power,
         waveform_model=fitted_model_step2.power,
         number_of_model_evaluations_step2=model_parameters_step2.num_ddm_evaluations,
         fit_return_status_step2=optimize_result_step2.status
    )


def samosa_fit_samosap_single_step2(
        waveform_data: NormedWaveform,
        scenario_data: ScenarioData,
        predictor: SAMOSAModelParameterPrediction,
        least_squares_kwargs: Optional[Dict] = None,
        sub_waveform_mask: Optional[np.ndarray] = None,
        amplitude_is_free_param: bool = True
) -> Tuple[WaveformModelParametersFit, WaveformModelOutput, OptimizeResult]:
    """
    Performs fit step 2 of the SAMOSAPlus retracker (as implemented in SAMPY)

    :param waveform_data:
    :param scenario_data:
    :param predictor:
    :param least_squares_kwargs:
    :param sub_waveform_mask:
    :param amplitude_is_free_param:

    :return: Fit result waveform mode
    """

    # Get the fit parameters
    first_guess, lower_bounds, upper_bounds = predictor.get(waveform_data)

    # Update least square kwargs with fit bounds
    # (Fit bounds are dynamic per waveform).
    fit_kwargs = dict(bounds=(lower_bounds, upper_bounds))
    fit_kwargs.update(least_squares_kwargs)

    # Second fit step in SAMOSA+ two-stage fits, which fits
    # three parameters: 1. epoch, 2. significant wave height, 3. Amplitude
    fit_cls = SAMOSAWaveformFit(
        scenario_data,
        waveform_data,
        waveform_model_kwargs=dict(mode=2, collect_fit_params=SAMOSA_WFM_COLLECT_FIT_PARAMS),
        sub_waveform_mask=sub_waveform_mask,
        amplitude_is_free_param=amplitude_is_free_param
    )
    fit_result = least_squares(fit_cls.fit_func, first_guess, **fit_kwargs)
    parameter_sdev = get_least_squares_parameter_sdev(fit_result)

    # --- Collect output ---
    # Summarize the fit result parameters
    model_parameters = WaveformModelParametersFit(
        epoch=fit_result.x[0] * 1e-9,
        epoch_sdev=float(parameter_sdev[0] * 1e-9),
        significant_wave_height=fit_result.x[1],
        significant_wave_height_sdev=float(parameter_sdev[1]),
        nu=fit_result.x[2],
        nu_sdev=float(parameter_sdev[2]),
        amplitude_scale=1.0,
        thermal_noise=waveform_data.thermal_noise,
        samosa_step="step2",
        num_ddm_evaluations=fit_cls.samosa_waveform_model.generate_ddm_counter
    )

    # Compute the waveform model with the fit parameters
    fitted_model = get_model_from_args(
        fit_cls.samosa_waveform_model,
        model_parameters.args_list,
        amplitude_scale=model_parameters.amplitude_scale,
        thermal_noise=model_parameters.thermal_noise
    )

    return model_parameters, fitted_model, fit_result


def samosa_fit_samosap_standard_step1(
        waveform_data: NormedWaveform,
        scenario_data: ScenarioData,
        predictor: SAMOSAModelParameterPrediction,
        step1_fixed_nu_value: float = 0.0,
        amplitude_is_free_param: bool = True,
        least_squares_kwargs: Optional[Dict] = None,
        sub_waveform_mask: Optional[np.ndarray] = None
) -> Tuple[WaveformModelParametersFit, WaveformModelOutput, OptimizeResult]:
    """
    Performs fit step 1 of the SAMOSAPlus retracker (as implemented in SAMPY)

    :param waveform_data:
    :param scenario_data:
    :param predictor:
    :param amplitude_is_free_param:
    :param step1_fixed_nu_value:
    :param least_squares_kwargs:
    :param sub_waveform_mask:

    :return: Waveform model parameters,
    """

    # Get the fit parameters
    first_guess, lower_bounds, upper_bounds = predictor.get(waveform_data, mode=1)

    # Update least square kwargs with fit bounds
    # (Fit bounds are dynamic per waveform).
    fit_kwargs = dict(bounds=(lower_bounds, upper_bounds))
    fit_kwargs.update(least_squares_kwargs)

    # First fit step in SAMOSA+ two-stage fits, which fits
    # three parameters: 1. epoch, 2. significant wave height, 3. Amplitude
    fit_cls = SAMOSAWaveformFit(
        scenario_data,
        waveform_data,
        waveform_model_kwargs=dict(mode=1, collect_fit_params=SAMOSA_WFM_COLLECT_FIT_PARAMS),
        step1_fixed_nu_value=step1_fixed_nu_value,
        sub_waveform_mask=sub_waveform_mask,
        amplitude_is_free_param=amplitude_is_free_param
    )
    fit_result_step1 = least_squares(fit_cls.fit_func_samosap_standard_step1, first_guess, **fit_kwargs)
    parameter_vars = get_least_squares_parameter_sdev(fit_result_step1)

    # --- Collect output ---
    # Summarize the fit result parameters
    amplitude_scale = fit_result_step1.x[2] if amplitude_is_free_param else 1.0
    model_parameters_step1 = WaveformModelParametersFit(
        epoch=fit_result_step1.x[0] * 1e-9,
        epoch_sdev=float(parameter_vars[0] * 1e-9),
        significant_wave_height=fit_result_step1.x[1],
        significant_wave_height_sdev=float(parameter_vars[1]),
        nu=step1_fixed_nu_value,
        nu_sdev=0.0,
        amplitude_scale=amplitude_scale,
        thermal_noise=waveform_data.thermal_noise,
        samosa_step="step1",
        num_ddm_evaluations=fit_cls.samosa_waveform_model.generate_ddm_counter
    )

    # Compute the waveform model with the fit parameters
    fitted_model_step1 = get_model_from_args(
        fit_cls.samosa_waveform_model,
        model_parameters_step1.args_list,
        amplitude_scale=model_parameters_step1.amplitude_scale,
        thermal_noise=model_parameters_step1.thermal_noise
    )

    return model_parameters_step1, fitted_model_step1, fit_result_step1


def samosa_fit_samosap_standard_step2(
        waveform_data: NormedWaveform,
        scenario_data: ScenarioData,
        predictor: SAMOSAModelParameterPrediction,
        amplitude_is_free_param: bool = True,
        step2_fixed_swh_value: float = 0.0,
        least_squares_kwargs: Optional[Dict] = None,
        sub_waveform_mask: Optional[np.ndarray] = None
) -> Tuple[WaveformModelParametersFit, WaveformModelOutput, OptimizeResult]:
    """
    Performs fit step 2 of the SAMOSAPlus retracker (as implemented in SAMPY)

    :param waveform_data:
    :param scenario_data:
    :param predictor:
    :param amplitude_is_free_param:
    :param step2_fixed_swh_value:
    :param least_squares_kwargs:
    :param sub_waveform_mask:

    :return: Fit result waveform mode
    """

    # Get the fit parameters
    first_guess, lower_bounds, upper_bounds = predictor.get(waveform_data, mode=2)

    # Update least square kwargs with fit bounds
    # (Fit bounds are dynamic per waveform).
    fit_kwargs = dict(bounds=(lower_bounds, upper_bounds))
    fit_kwargs.update(least_squares_kwargs)

    # Second fit step in SAMOSA+ two-stage fits, which fits
    # three parameters: 1. epoch, 2. significant wave height, 3. Amplitude
    fit_cls = SAMOSAWaveformFit(
        scenario_data,
        waveform_data,
        waveform_model_kwargs=dict(mode=2, collect_fit_params=SAMOSA_WFM_COLLECT_FIT_PARAMS),
        step2_fixed_swh_value=step2_fixed_swh_value,
        sub_waveform_mask=sub_waveform_mask,
        amplitude_is_free_param=amplitude_is_free_param
    )
    fit_result_step2 = least_squares(fit_cls.fit_func_samosap_standard_step2, first_guess, **fit_kwargs)
    parameter_sdev = get_least_squares_parameter_sdev(fit_result_step2)

    # --- Collect output ---
    # Summarize the fit result parameters
    amplitude_scale = fit_result_step2.x[2] if amplitude_is_free_param else 1.0
    model_parameters_step2 = WaveformModelParametersFit(
        epoch=fit_result_step2.x[0] * 1e-9,
        epoch_sdev=float(parameter_sdev[0] * 1e-9),
        significant_wave_height=step2_fixed_swh_value,
        significant_wave_height_sdev=0.0,
        nu=fit_result_step2.x[1],
        nu_sdev=float(parameter_sdev[1]),
        amplitude_scale=amplitude_scale,
        thermal_noise=waveform_data.thermal_noise,
        samosa_step="step2",
        num_ddm_evaluations=fit_cls.samosa_waveform_model.generate_ddm_counter
    )

    # Compute the waveform model with the fit parameters
    fitted_model_step2 = get_model_from_args(
        fit_cls.samosa_waveform_model,
        model_parameters_step2.args_list,
        amplitude_scale=model_parameters_step2.amplitude_scale,
        thermal_noise=model_parameters_step2.thermal_noise
    )

    # Compute uncertainties (Standard devi
    return model_parameters_step2, fitted_model_step2, fit_result_step2


def get_model_from_args(
        samosa_waveform_model,
        args,
        amplitude_scale: float = 1.0,
        thermal_noise: float = 0.0,
) -> "WaveformModelOutput":
    """
    Compute waveform model with final fit parameters. This function
    allows to compute the waveform model from fit parameter triplets.

    :param samosa_waveform_model: Initialized and configured SAMOSA+
        Waveform model instance.

    :param args: list of [epoch, swh, mss]
    :param amplitude_scale: Optional scaling factor
    :param thermal_noise: Optional thermal noise

    :return: Waveform model
    """
    model_parameter = WaveformModelParameters(
        epoch=args[0],
        significant_wave_height=args[1],
        nu=args[2],
        amplitude_scale=amplitude_scale,
        thermal_noise=thermal_noise
    )
    return samosa_waveform_model.generate_delay_doppler_waveform(model_parameter)


def total_velocity_from_vector(
        velocity_x: float,
        velocity_y: float,
        velocity_z: float
) -> float:
    """
    Compute total velocity from vector

    :param velocity_x: velocity x-component
    :param velocity_y: velocity y-component
    :param velocity_z: velocity z-component

    :return: Total velocity
    """
    return np.sqrt(velocity_x**2. + velocity_y**2. + velocity_z**2.)


def get_valid_epoch_range(
        epoch: np.ndarray,
        earliest_fraction: float = 0.0,
        latest_fraction: float = 1.0
) -> Tuple[float, float]:
    """
    Get bounds for epoch based on the valid fraction range in the range window.

    :param epoch: epoch array
    :param earliest_fraction: first valid epoch as fraction of the range window
    :param latest_fraction: last valid epoch as fraction of the range window

    :raises ValueError: fraction not in interval [0-1] or earliest > latest

    :return: (epoch lower bound, epoch upper bound)
    """

    if earliest_fraction > latest_fraction:
        raise ValueError(f"{earliest_fraction} not smaller {latest_fraction}")

    for value in [earliest_fraction, latest_fraction]:
        if not 0.0 <= value <= 1.0:
            raise ValueError(f"{value} out of bounds [0-1]")

    window_size = epoch[-1] - epoch[0]
    return float(epoch[0] + earliest_fraction * window_size), float(epoch[0] + latest_fraction * window_size)


def get_epoch_bounds(
        tau: np.ndarray,
        epoch_first_guess_ns: float,
        tau_earliest_fraction: float = 0.02,
        tau_latest_fraction: float = 0.8,
        range_gates_after_fmi: int = 20
) -> Tuple[float, float]:
    """
    Compute the fit bounds for the epoch. The lower epoch bound is defined by
    a fraction of the range window, while the upper epoch bound is defined as
    a certain number of range gates after the first guess (first maximum index).

    :param tau: range gates in nanoseconds
    :param epoch_first_guess_ns: First guess for the epoch (ususally location of first maximum)
    :param tau_earliest_fraction: Fraction of the range window for the earliest epoch bound
    :param tau_latest_fraction: Fraction of the range window for the latest epoch bound (not used)
    :param range_gates_after_fmi: Number of trailing range gates after the epoch
        first guess (defines the upper epoch bound)
    """
    epoch_bounds_range_gate = get_valid_epoch_range(tau, tau_earliest_fraction, tau_latest_fraction)
    tau_fmi_offset = (epoch_first_guess_ns * 1e-9) + range_gates_after_fmi * (tau[1] - tau[0])
    return (
        epoch_bounds_range_gate[0],
        np.min([np.max(tau), tau_fmi_offset])
    )


def sampy_misfit(
        residuals: np.ndarray,
        waveform_mask: Optional[np.ndarray] = None,
        waveform_scale: float = None
) -> float:
    """
    Computes the SAMOSA waveform model misfit parameter according to SAMPy with optional
    misfit computation on sub-waveform.

    One deviation from the SAMPy approach is the optional misfit reduction due to
    waveform maximum power > 1. This may happen when the first maximum power is
    less than the maximum power. In this case the misfit value is artificially
    increased compared to the SAMPy approach and scenarious where the first maximum
    is also the absolute maximum power. This cannot be derived from the residuals
    alone and therefore need to be provided as apriori information.

    :param residuals: difference between waveform and waveform model
    :param waveform_mask: numpy index array of sub-waveform mask
    :param waveform_scale: Maximum power the waveform. (Only needed if
        maximum power > 1).

    :return: SAMOSA waveform model misfit
    """
    waveform_mask = waveform_mask if waveform_mask is not None else np.arange(residuals.size)
    misfit = np.sqrt(1. / residuals[waveform_mask].size * np.nansum(residuals[waveform_mask] ** 2)) * 100.
    return misfit if waveform_scale is None else misfit / waveform_scale


def compute_thermal_noise(
        waveform: np.ndarray,
        start_index: int = 4,
        end_index: int = 12,
) -> float:
    """
    Compute thermal noise for waveform collection (from SAMPy). The strategy is
    to compute the media in a slice of the lowest values of the first half
    of the waveform (apparently the end of the waveform contains zeros).

    :param waveform: The waveform array
    :param start_index: start of the slice of sorted waveform power values (forced to valid range)
    :param end_index: end of the slice of sorted waveform power values (forced to valid range)

    :return: thermal noise in waveform power units
    """
    waveform_sorted = np.sort(waveform[:waveform.shape[0] // 2])
    start_index = max(start_index, 0)
    end_index = waveform_sorted.size - 1 if start_index >= waveform.size else end_index
    return float(np.nanmedian(waveform_sorted[start_index:end_index + 1]))


def epoch2range(epoch: float, range_array: np.ndarray) -> float:
    """
    Computes the retracker range, defined as range of spacecraft center or mass
    to retracked elevation.

    :param epoch: retracker epoch in seconds

    :param range_array:

    :return: retracker range in meter.
    """
    factor = 0.5 * 299792458.
    center_idx = range_array.shape[0] // 2
    center_range = float(range_array[center_idx])
    return epoch * factor + center_range


def get_look_angles(l1, index: int) -> np.ndarray:
    """
    Compute the lookangles based on l1b stack information
    # TODO: This functions raises a ValueError for NaN values in the classifiers (LRM?)

    :param l1: The level 1 data object
    :param index: Waveform index

    :return:
    """
    return 90. - np.linspace(np.rad2deg(l1.classifier.look_angle_start[index]),
                             np.rad2deg(l1.classifier.look_angle_stop[index]),
                             num=int(l1.classifier.stack_beams[index]))


def get_least_squares_parameter_sdev(optimize_result: "OptimizeResult") -> np.ndarray:
    """
    Compute standard deviation of parameter from Jacobian.

    Solution from Stack Overflow: https://stackoverflow.com/questions/42388139

    :param optimize_result: Return object from `scipy.optimize.least_squares`

    :return: List of fit parameter standard deviations
    """
    jac = optimize_result.jac
    try:
        cov = np.linalg.inv(jac.T.dot(jac))
    except np.linalg.LinAlgError:
        # If the Jacobian is singular, we cannot compute the covariance matrix
        # and therefore return NaN for all parameters.
        return np.full(optimize_result.x.size, np.nan)
    variance = np.sqrt(np.diagonal(cov))
    return np.sqrt(variance)


def get_sub_waveform_mask(waveform_data: NormedWaveform, filter_trailing_edge_kwargs: Dict) -> np.ndarray:
    """
    Estimates a sub-waveform mask by flagging off-nadir backscatter elements
    that manifests as peaks on the trailing edge. This is done by estimating
    the lower envelope of the trailing edge and exclude all waveform power
    values that

    :param waveform_data: Waveform data
    :param filter_trailing_edge_kwargs: Configuration keyword arguments for the
        trailing edge filter.

    :return: sub-waveform mask (True: masked)
    """
    sub_waveform_mask = get_trailing_edge_lower_envelope_mask(
        waveform_data.power.copy(),
        waveform_data.first_maximum_index,
        **filter_trailing_edge_kwargs
    )
    return np.logical_not(sub_waveform_mask)


def get_trailing_edge_lower_envelope_mask(
        waveform_power: np.ndarray,
        first_maximum_index: int,
        max_minimum_cleaning_passes: int = 5,
        noise_level_normed: float = 0.02,
        return_type: Literal["bool", "indices"] = "bool"
) -> np.ndarray:
    """
    Compute the indices of the waveforms trailing edge lower envelope.
    This method relies on the computation of relative minima via scipy.
    with added filtering.

    :param waveform_power: Waveform power values (any unit)
    :param first_maximum_index: Range gate index of the first maximum.
    :param max_minimum_cleaning_passes: Local minima on the trailing edge
        may not be in strictly decreasing in power depending on the
        off-nadir backscatter contanimation of the trailing edge.
        Local minima are iteratively removed that are not decreasing in
        power.
    :param noise_level_normed: Points within the noise level of the linear
        fit between successive local minima are added to the trailing edge
        mask. Noise level unit is the first maximum power.
    :param return_type: Options are a boolean array with the same dimensions
        as the waveform (default), or an index array.

    :return:
    """

    # Perform operation on trailing edge only normalized by
    # first maximum power (trailing edge is defined as
    # everything after the first maximum).
    wfm_trailing_edge = waveform_power[first_maximum_index:]
    wfm_trailing_edge /= np.nanmax(wfm_trailing_edge)

    # Get local minima and add first maximum and last range gate
    # the list (important for fitting waveform power).
    idx_lmin = argrelmin(wfm_trailing_edge)[0]
    idx_lmin = np.insert(idx_lmin, 0, 0)
    idx_lmin = np.insert(idx_lmin, len(idx_lmin), wfm_trailing_edge.size-1)

    # Throw out local minima that have greater power then the previous
    # one in several iterations
    for _ in np.arange(max_minimum_cleaning_passes):
        wfm_power_change = [
            wfm_trailing_edge[idx_lmin[1:]] - wfm_trailing_edge[idx_lmin[:-1]]
        ][0]
        wfm_change = np.insert(wfm_power_change, 0, -1)
        valid_local_minima = wfm_change < 0.0
        if valid_local_minima.all():
            break
        idx_lmin = idx_lmin[wfm_change < 0.0]

    # Create a leading edge mask. This is required for the
    # following operation
    mask_le = np.full(wfm_trailing_edge.size, False)
    mask_le[idx_lmin] = True

    # Take in again waveform points that are within the noise level
    # (ensures that are enough points for the fit, without deviating
    # from the lower envelope too much).
    # This is done by computing a linear fit between successive
    # local minima and checking the difference of actual power
    # values to the linear fit.
    n_range_gates = len(waveform_power)
    for i in np.arange(idx_lmin.size - 1):

        # --- Method 1: Distance to POCA of linear fit ---
        idx_lmin_scaled = idx_lmin.astype(float)/ n_range_gates
        k = wfm_trailing_edge[idx_lmin[i]]
        m = ((wfm_trailing_edge[idx_lmin[i + 1]] - wfm_trailing_edge[idx_lmin[i]]) /
             (idx_lmin_scaled[i + 1] - idx_lmin_scaled[i]))

        for x, idx in enumerate(np.arange(idx_lmin[i] + 1, idx_lmin[i + 1])):

            x0 = float(x + 1) / n_range_gates
            y0 = wfm_trailing_edge[idx]

            x_poca = (x0 + m * y0 - m * k) / (m ** 2. + 1)
            y_poca = m * (x0 + m * y0 - m * k) / (m ** 2. + 1) + k

            dist = np.sqrt((x_poca - x0)**2. + (y_poca - y0)**2)
            mask_le[idx] = dist <= noise_level_normed

    # Final step: Pad trailing edge mask to full waveform mask
    mask = np.full(waveform_power.size, True)
    mask[first_maximum_index:] = mask_le

    return mask if return_type == "bool" else np.where(mask)[0]


def get_samosa_leading_edge_error(
        waveform_model: np.ndarray,
        waveform: np.ndarray,
        leading_edge_start_power_threshold: float = 0.005
) -> float:
    """
    Compute the SAMOSA+ leading edge error value from the RMSE of the
    leading edge of the waveform and the fitted SAMOSA+ waveform model.

    This error is an additional quality control measure for the SAMOSA+
    range e.g., in cases where the overall fit seems to be nominal,
    but the leading edge (and its tracking point) are not well represented.

    The leading edge will be estmated from the SAMOSA+ waveform model
    (because this is a lot easier than from the actual waveform).

    :param waveform_model: The SAMOSA+ waveform model.
    :param waveform: The actual waveform
    :param leading_edge_start_power_threshold: Power threshold in units of
        the maximum power of the SAMOSA+ waveform model to determine the
        start of the leading edge.

    :return: SAMOSA+ leading edge error value (RMSE between waveform and
        waveform model for leading edge only).
    """
    waveform_model_normed = waveform_model / np.nanmax(waveform_model)
    try:
        leading_edge_start = np.where(waveform_model_normed > leading_edge_start_power_threshold)[0][0]
    except IndexError:
        leading_edge_start = 0
    model_maximum = np.argmax(waveform_model)
    return rmse(waveform_model[leading_edge_start:model_maximum+1], waveform[leading_edge_start:model_maximum+1])


def rmse(predictions: np.ndarray, targets: np.ndarray) -> float:
    """
    Compute root-mean-square-error (rmse) of 2 non-Nan arrays

    :param predictions: predicted values (non-NaN)
    :param targets: target values (non-NaN)

    :return: rmse value
    """
    return np.sqrt(((predictions - targets) ** 2).mean())


def get_nu_from_ocog_width(
            ocog_width: float,
            nu_ocog_coefs: Tuple[float, float],
            ocog_width_max: float = 50,
            nu_min: float = 0,
            nu_max: float = 1e6,
            default_nu: float = 0.0
) -> float:
    """
    Compute the inverse mean square slope (nu) from an empirical parametrization
    based on ocog width. For stability reasons, the nu value is clipped to
    a minimum and maximum value.

    :param ocog_width: OCOG width (in range gates)
    :param nu_ocog_coefs: coeficients for the OCOG width to nu relation (c[0] * ocog_width ** (-c[1]))
    :param ocog_width_max: The maximum OCOG value for which this relation is valid
    :param nu_min: Minimum value for nu
    :param nu_max: Maximum value for nu
    :param default_nu: Default value for nu if the OCOG width is above the maximum value

    :return: Inverse mean square slope (nu)
    """
    if ocog_width > ocog_width_max:
        return default_nu
    nu_predicted = nu_ocog_coefs[0] + ocog_width ** (-1.0 * nu_ocog_coefs[1])
    return np.clip(nu_predicted, nu_min, nu_max)



# def calc_sigma0(
#         Latm,
#         Pu,
#         CST,
#         RDB,
#         GEO,
#         epoch_sec,
#         window_del_20_hr_ku_deuso,
#         lat_20_hr_ku,
#         alt_20_hr_ku,
#         sat_vel_vec_20_hr_ku,
#         transmit_pwr_20_ku
# ):
#     atmospheric_attenuation = 10. ** (np.zeros(np.shape(alt_20_hr_ku)) / 10.)  ### IMP!!! ===>  Atmospheric Correction, set to zero dB
#     # Atm_Atten    = 10.**(Latm/10.)
#     Range = epoch_sec * CST.c0 / 2 + CST.c0 / 2 * window_del_20_hr_ku_deuso  #### this should be the retracked range without geo-corrections
#     Pout = Pu
#     earth_radius = np.sqrt(CST.R_e ** 2.0 * (np.cos(np.deg2rad(lat_20_hr_ku))) ** 2 + CST.b_e ** 2.0 * (
#         np.sin(np.deg2rad(lat_20_hr_ku))) ** 2)
#     kappa = (1. + alt_20_hr_ku / earth_radius)
#
#     Lx = CST.c0 * Range / (2. * sat_vel_vec_20_hr_ku * RDB.f_0 * RDB.Np_burst * 1. / RDB.PRF_SAR)
#     Ly = np.sqrt(CST.c0 * Range * (1. / RDB.Bs) / kappa)
#     A_SAR = Lx * (2. * Ly)
#     C = ((4. * np.pi) ** 3. * (GEO.Height ** 4) * atmospheric_attenuation) / (((CST.c0 / RDB.f_0) ** 2) * (RDB.G_0 ** 2) * A_SAR)
#
#     sigma0 = 10. * np.log10(Pout / transmit_pwr_20_ku) + 10. * np.log10(C) + RDB.bias_sigma0
#
#     return sigma0, np.log10(Pout / transmit_pwr_20_ku), np.log10(C), Range, kappa
