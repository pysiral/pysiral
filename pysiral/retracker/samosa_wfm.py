# -*- coding: utf-8 -*-

"""
This is the SAMOSA+ retracker variant that uses `samosa-waveform-model`, a
re-implementation of the SAMPy package. The objective of the re-implementation
is better performance (code optimization and multiprocessing) and a greater flexibility
for retracking settings (sub-waveform retracking, limiting parameters of the fit)
"""

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"

import numpy as np
import numpy.typing as npt
from typing import Tuple, List
from dataclasses import dataclass
from scipy.optimize import least_squares

from samosa_waveform_model import (
    ScenarioData, SensorParameters, SAMOSAWaveformModel, PlatformLocation, SARParameters,
    WaveformModelParameters
)

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

@dataclass
class SAMOSAWaveformFitResult:
    epoch: float
    significant_waveheight: float
    mean_square_slope: float
    thermal_noise: float
    misfit: float
    fit_mode: str
    waveform: np.ndarray
    waveform_fit: np.ndarray
    is_sub_waveform_fit: bool = False
    misfit_sub_waveform: float = None
    number_of_fit_iterations: int = -1


class SAMOSAWaveformModelFit(BaseRetracker):

    def __init__(self) -> None:
        super(SAMOSAWaveformModelFit, self).__init__()

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
        breakpoint()

    def _samosa_plus_retracker(self, rng, wfm, indices, radar_mode) -> List[SAMOSAWaveformFitResult]:
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
        normed_waveform = self._get_normed_waveform(rng, wfm, radar_mode)
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

    def _get_altitude_velocity_from_l1(self) -> Tuple[npt.NDArray, npt.NDArray]:
        """
        Get altitude (height) rate and satellite velocity from l1 data

        :return:
        """
        hrate = np.zeros_like(self._l1b.time_orbit.altitude_rate)
        vel = np.sqrt(self.get_l1b_parameter("classifier", "satellite_velocity_x")**2
                      + self.get_l1b_parameter("classifier", "satellite_velocity_y")**2
                      + self.get_l1b_parameter("classifier", "satellite_velocity_z")**2)
        return hrate, vel

    @staticmethod
    def _get_normed_waveform(rng: np.ndarray, wfm: np.ndarray, radar_mode:int) -> NormedWaveform:
        """
        Return a normed representation of the input waveform

        :param rng:
        :param wfm:
        :param radar_mode:

        :return:
        """
        scaling_factor = 1.0 / np.nanmax(wfm)
        normed_waveform = wfm * scaling_factor
        return NormedWaveform(normed_waveform, rng, radar_mode, scaling_factor)

    def _model_fit(self, normed_waveform: NormedWaveform, scenario_data: ScenarioData) -> SAMOSAWaveformFitResult:
        """
        Fit the waveform

        :param normed_waveform:
        :param scenario_data:
        :return:
        """


        # Get first guess
        epoch_first_guess = scenario_data.rp.tau[np.argmax(normed_waveform.power)]
        epoch_bounds = get_epoch_bounds(scenario_data.rp.tau, 0.1, 0.8)
        first_guess = [epoch_first_guess, -0.2, 5e-7]
        lower_bounds = [epoch_bounds[0], -0.3, 0.0]
        upper_bounds = [epoch_bounds[1], 0.5, 1e-5]

        fit_kwargs = dict(
            bounds = (lower_bounds, upper_bounds),
            loss = "linear",
            method = "trf",
            ftol = 1e-2,
            xtol = 2e-3,
            gtol = 1e-2,
            max_nfev = None,
        )

        fit_cls = SAMOSAWaveformFit(scenario_data, normed_waveform)
        fit_result = least_squares(fit_cls.fit_func, first_guess, **fit_kwargs)

        breakpoint()


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
        self.last_waveform_model = None

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
        self.last_waveform_model = self.samosa_waveform_model.generate_delay_doppler_waveform(model_parameter)
        return self.last_waveform_model.power - self.normed_waveform.power


def total_velocity_from_vector(
        velocity_x: float,
        velocity_y: float,
        velocity_z: float
) -> float:
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
    return epoch[0] + earliest_fraction * window_size, epoch[0] + latest_fraction * window_size
