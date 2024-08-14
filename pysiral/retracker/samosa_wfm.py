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
from typing import Tuple

from sympy.physics.units import velocity

from samosa_waveform_model import ScenarioData, SensorParameters, SAMOSAWaveformModel, PlatformLocation, SARParameters

from pysiral.retracker import BaseRetracker


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


    def _samosa_plus_retracker(self, rng, wfm, indices, radar_mode):

        for idx in indices:
            fit_result = self._samosa_wfm_retrack(rng[idx], wfm[idx,:], radar_mode[idx], idx)
        breakpoint()


    def _samosa_wfm_retrack(self, rng: float, wfm: np.ndarray, radar_mode:int, idx: int):

        # Create the SAMOSA Waveform scenario data

        sp = SensorParameters.cryosat2_sar()

        location_data = dict(
            latitude=self._l2.latitude[idx],
            longitude=self._l2.longitude[idx],
            altitude=self._l2.altitude[idx],
            height_rate=self._l1b.time_orbit.altitude_rate[idx],
            pitch=self._l1b.time_orbit.antenna_pitch[idx],
            roll=self._l1b.time_orbit.antenna_roll[idx],
            velocity=total_velocity_from_vector(
                self._l1b.classifier.satellite_velocity_x[idx],
                self._l1b.classifier.satellite_velocity_y[idx],
                self._l1b.classifier.satellite_velocity_z[idx]
            )
        )

        geo = PlatformLocation(**location_data)
        sar = SARParameters()
        sar.compute_multi_look_parameters(geo=geo, sp=sp)

        wfm_scenario = ScenarioData(sp, geo, sar)
        breakpoint()

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


def total_velocity_from_vector(velocity_x: float, velocity_y: float, velocity_z: float) -> float:
    return np.sqrt(velocity_x**2. + velocity_y**2. + velocity_z**2.)
