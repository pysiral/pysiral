# -*- coding: utf-8 -*-

"""
This is the SAMOSA+ retracker variant that uses `samosa-waveform-model`, a
re-implementation of the SAMPy package. The objective of the re-implementation
is better performance (code optimization and multiprocessing) and a greater flexibility
for retracking settings (sub-waveform retracking, limiting parameters of the fit)
"""

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"

import numpy as np
from typing import Tuple, List, Optional, Any
from dataclasses import dataclass
from scipy.optimize import least_squares

from samosa_waveform_model import (
    ScenarioData, SensorParameters, SAMOSAWaveformModel, PlatformLocation, SARParameters,
    WaveformModelParameters
)
from samosa_waveform_model.dataclasses import WaveformModelOutput

from pysiral.core.config import RadarModes
from pysiral.retracker import BaseRetracker
from pysiral.l1data import Level1bData
from pysiral.l2data import Level2Data


@dataclass
class NormedWaveform:
    """
    Data container for the waveform fit
    """
    power: np.ndarray
    range_bins: np.ndarray
    scaling_factor: float
    radar_mode_flag: int
    window_delay: float
    transmit_power: float

@dataclass
class SAMOSAWaveformFitResult:
    epoch: float
    significant_wave_height: float
    mean_square_slope: float
    thermal_noise: float
    misfit: float
    fit_mode: str
    waveform: np.ndarray
    waveform_fit: np.ndarray
    is_sub_waveform_fit: bool = False
    misfit_sub_waveform: float = None
    number_of_fit_iterations: int = -1


class SAMOSAWaveformFit(object):
    """
    Class for fitting SAMOSA waveform model with scipy.optimize.least_squares

    Usage:

        fit_cls = SAMOSAWaveformFit(scenario_data, normed_waveform, ...)
        fit_result = scipy.optimize.least_squares(fit_cls.fit_func, ...)
    """

    def __init__(
            self,
            scenario_data: ScenarioData,
            normed_waveform: NormedWaveform,
            sub_waveform_mask: np.ndarray = None,
            method: str = None
    ) -> None:
        self.samosa_waveform_model = SAMOSAWaveformModel(scenario_data)
        self.normed_waveform = normed_waveform
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
        epoch, significant_wave_height, mean_square_slope = fit_args
        model_parameter = WaveformModelParameters(
            epoch=epoch,
            significant_wave_height=significant_wave_height,
            amplitude=1.0,
            mean_square_slope=mean_square_slope
        )
        waveform_model = self.samosa_waveform_model.generate_delay_doppler_waveform(model_parameter)

        # TODO: Computation of residuals should be configurable method
        # TODO: Add thermal noise to computation of residuals
        return waveform_model.power - self.normed_waveform.power


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
        fit_results = self._samosa_plus_retracker(rng, wfm, indices, radar_mode)

        # Store retracker properties (including range)
        self._store_retracker_properties(fit_results, indices)

        # Set/compute uncertainty
        # self._set_range_uncertainty()

        # Add range biases (when set in config file)
        # self._set_range_bias(radar_mode)

        # Add auxiliary variables to the l2 data object
        # self._register_auxiliary_variables()

    def create_retracker_properties(self, n_records: int) -> None:
        """
        Initialize retracker properties with correct arrays shapes (shape = (n_records, )).
        The list of parameters depends on whether the SAMOSA_DEBUG_MODE flag is set.

        NOTE: The properties are set to an array, but can be accessed as `self.{property_name}`
        via the __getattr__ method.

        :param n_records:
        """

        parameter = [
            "misfit",
            "swh",
            "mean_square_slope",
            "wind_speed",
            "epoch",
            "guess",
            "Pu",
            "rval",
            "kval",
            "pval",
            "cval",
        ]
        for parameter_name in parameter:
            self._retracker_params[parameter_name] = np.full(n_records, np.nan, dtype=np.float32)

    def _samosa_plus_retracker(self, rng, wfm, indices, radar_mode) -> List[SAMOSAWaveformFitResult]:
        """
        TODO: Placeholder implementation (needs option to do multi-processing)
        :param rng:
        :param wfm:
        :param indices:
        :param radar_mode:
        :return:
        """
        return [self._retrack_waveform(rng[idx], wfm[idx,:], radar_mode[idx], idx) for idx in indices]

    def _retrack_waveform(
            self,
            rng: np.ndarray,
            wfm: np.ndarray,
            radar_mode:int,
            idx: int
    ) -> SAMOSAWaveformFitResult:
        """
        Sandbox version of the waveform fit.
        TODO: To be replaced with configurable fit methods

        :param rng:
        :param wfm:
        :param radar_mode:
        :param idx:

        :return: Fit result container
        """
        # Create the SAMOSA Waveform scenario data
        scenario_data = self._get_scenario_data(self._l1b, self._l2, idx)
        normed_waveform = self._get_normed_waveform(
            rng,
            wfm,
            radar_mode,
            self._l1b.classifier.window_delay[idx],
            self._l1b.classifier.transmit_power[idx]
        )
        return self._model_fit(normed_waveform, scenario_data)

    def _get_scenario_data(
            self,
            l1: Level1bData,
            l2: Level2Data,
            idx: int
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

        :return: SAMOSA waveform model scenario data for specific waveform
        """

        radar_mode_name = RadarModes.get_name(l1.waveform.radar_mode[idx])
        platform = l2.info.mission

        sp = SensorParameters.get(platform, radar_mode_name)

        location_data = dict(
            latitude=l2.latitude[idx],
            longitude=l2.longitude[idx],
            altitude=l2.altitude[idx],
            height_rate=l1.time_orbit.altitude_rate[idx],
            pitch=l1.time_orbit.antenna_pitch[idx],
            roll=l1.time_orbit.antenna_roll[idx],
            velocity=total_velocity_from_vector(
                self._l1b.classifier.satellite_velocity_x[idx],
                self._l1b.classifier.satellite_velocity_y[idx],
                self._l1b.classifier.satellite_velocity_z[idx]
            )
        )

        geo = PlatformLocation(**location_data)
        sar = SARParameters()
        sar.compute_multi_look_parameters(geo=geo, sp=sp)

        return ScenarioData(sp, geo, sar)

    @staticmethod
    def _get_normed_waveform(
            rng: np.ndarray,
            wfm: np.ndarray,
            radar_mode:int,
            window_delay: float,
            transmit_power: float
    ) -> NormedWaveform:
        """
        Return a normed representation of the input waveform

        :param rng:
        :param wfm:
        :param radar_mode:

        :return:
        """
        scaling_factor = 1.0 / np.nanmax(wfm)
        normed_waveform = wfm * scaling_factor
        return NormedWaveform(
            normed_waveform,
            rng,
            radar_mode,
            scaling_factor,
            window_delay,
            transmit_power
        )

    def _model_fit(
            self,
            normed_waveform: NormedWaveform,
            scenario_data: ScenarioData
    ) -> SAMOSAWaveformFitResult:
        """
        Fit the SAMOSA waveform model to the waveform and computes additional
        fit parameters.

        :param normed_waveform:
        :param scenario_data:

        :return: fit result dataclass
        """

        # Get first guess
        epoch_first_guess = scenario_data.rp.tau[np.argmax(normed_waveform.power)]
        epoch_bounds = get_epoch_bounds(scenario_data.rp.tau, 0.1, 0.8)
        first_guess = [epoch_first_guess, -0.2, 5e-7]
        lower_bounds = [epoch_bounds[0], -0.3, 0.0]
        upper_bounds = [epoch_bounds[1], 0.5, 1e-6]

        fit_kwargs = dict(
            bounds = (lower_bounds, upper_bounds),
            loss = "linear",
            method = "trf",
            ftol = 1e-2,
            xtol = 2e-4,
            gtol = 1e-2,
            max_nfev = None,
        )

        fit_cls = SAMOSAWaveformFit(scenario_data, normed_waveform)
        fit_result = least_squares(fit_cls.fit_func, first_guess, **fit_kwargs)
        fitted_model = self._get_fitted_model(fit_cls.samosa_waveform_model, fit_result.x)
        misfit = sampy_misfit(fit_result.fun)

        return SAMOSAWaveformFitResult(
            fit_result.x[0],
            fit_result.x[1],
            fit_result.x[2],
            0.0,   # placeholder for thermal noise
            misfit,
            "test",
            normed_waveform.power,
            fitted_model.power
        )

    @staticmethod
    def _get_fitted_model(samosa_waveform_model, fit_coefs) -> "WaveformModelOutput":
        """
        Compute waveform model with final fit parameters

        :param samosa_waveform_model:
        :param fit_coefs: Final set of fit coefficients returned by

        :return:
        """
        model_parameter = WaveformModelParameters(
            epoch=fit_coefs[0],
            significant_wave_height=fit_coefs[1],
            amplitude=1.0,
            mean_square_slope=fit_coefs[2]
        )
        return samosa_waveform_model.generate_delay_doppler_waveform(model_parameter)

    def _store_retracker_properties(
        self,
        fit_results: List[SAMOSAWaveformFitResult],
        indices: np.ndarray
    ) -> None:
        """
        Store the output of the SAMOSA+ retracker in the class.

        :param fit_results:
        :param indices:

        :return:
        """

        for index, fit_result in zip(indices, fit_results):

            self._range[index] = fit_result.rng
            self._power[index] = fit_result.sigma0

            # Store additional retracker parameters
            self.swh[index] = fit_result.significant_wave_height
            self.misfit[index] = fit_result.misfit
            # self.wind_speed[index] = func_wind_speed([fit_result.sigma0])
            self.mean_square_slope[index] = fit_result.mean_square_slope
            self.epoch[index] = fit_result.epoch


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


def get_epoch_bounds(
        epoch: np.ndarray,
        earliest_fraction: float,
        latest_fraction: float
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


def sampy_misfit(residuals: np.ndarray, waveform_mask: Optional[np.ndarray] = None) -> float:
    """
    Computes the SAMOSA waveform model misfit parameter according to SAMPy

    :param residuals: difference between waveform and waveform model
    :param waveform_mask:

    :return: SAMOSA waveform model misfit
    """
    return np.sqrt(1. / residuals.size * np.nansum(residuals ** 2)) * 100.