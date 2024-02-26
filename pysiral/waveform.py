# -*- coding: utf-8 -*-
"""
Created on Fri Jul 01 13:07:10 2016

@author: shendric
"""

import multiprocessing
import multiprocessing.pool
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Union

import bottleneck as bn
import numpy as np
import numpy.typing as npt
from loguru import logger
from schema import And, Schema
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.signal import argrelmin

from pysiral import psrlcfg
from pysiral.core.clocks import StopWatch
from pysiral.core.config import RadarModes
from pysiral.core.functions import inverse_power
from pysiral.core.helper import get_multiprocessing_1d_array_chunks
from pysiral.l1data import Level1bData
from pysiral.l1preproc.procitems import L1PProcItem
from pysiral.retracker.tfmra import cTFMRA


def get_waveforms_peak_power(wfm: npt.NDArray, use_db: bool = False) -> npt.NDArray:
    """
    Return the peak power (in input coordinates) of an array of waveforms

    Arguments
    ---------
        wfm (float array)
            echo waveforms, order = (n_records, n_range_bins)

    Returns
    -------
        float array with maximum for each echo
    """
    peak_power = np.amax(wfm, axis=1)
    if use_db:
        peak_power = 10 * np.log10(peak_power)
    return peak_power


def get_footprint_pulse_limited(r: float, band_width: float) -> float:
    """
    Compute the CryoSat-2 LRM footprint for variable range to the surface.

    Applicable Documents:

        Michele Scagliola, CryoSat Footprints ESA/Aresys, v1.2
        ESA document ref: XCRY-GSEG-EOPG-TN-13-0013

    :param r: range from satellite center of mass to surface reflection point
    :param band_width: pulse bandwidth in Hz

    :return area_pl: Radar backscatter coefficient
    """

    c_0 = 299792458.0
    footprint_radius = np.sqrt(r * c_0 / band_width) / 1000.
    return np.pi * footprint_radius ** 2.


def get_footprint_sar(r: float,
                      v_s: float,
                      ptr_width: float,
                      tau_b: float,
                      lambda_0: float,
                      wf: float = 1.0,
                      r_mean: float = 6371000.0
                      ) -> float:
    """
    Compute the radar footprint of SAR altimeters for variable range to the surface.

    Applicable Documents:

        Michele Scagliola, CryoSat Footprints ESA/Aresys, v1.2
        ESA document ref: XCRY-GSEG-EOPG-TN-13-0013

    :param r: range from satellite center of mass to surface reflection point
    :param v_s: satellite along track velocity in meter/sec
    :param ptr_width: 3dB range point target response temporal width in seconds
    :param tau_b: burst length in seconds
    :param lambda_0: radar wavelength in meter
    :param wf: footprint widening factor
    :param r_mean: mean earth radius in meter

    :return area_sar: The SAR footprint in square meters
    """
    c_0 = 299792458.0
    alpha_earth = 1. + (r / r_mean)
    lx = (lambda_0 * r) / (2. * v_s * tau_b)
    ly = np.sqrt((c_0 * r * ptr_width) / alpha_earth)
    return (2. * ly) * (wf * lx)


def get_sigma0_sar(rx_pwr: float,
                   tx_pwr: float,
                   r: float,
                   a: float,
                   lambda_0: float,
                   g_0: float,
                   l_atm: float = 1.0,
                   l_rx: float = 1.0,
                   bias_sigma0: float = 0.0,
                   ) -> float:
    """
    Compute the sigma0 backscatter coefficient according the radar equation, e.g.
    equation 20 in

        Guidelines for reverting Waveform Power to Sigma Nought for
        CryoSat-2 in SAR mode (v2.2), Salvatore Dinardo, 23/06/2016
        XCRY-GSEG-EOPS-TN-14-0012

    :param rx_pwr: received power (waveform maximumn)
    :param tx_pwr: transmitted power
    :param r: range to surface in ms
    :param a: area illuminated by the altimeter
    :param lambda_0: radar wavelength in meter
    :param g_0: antenna gain factor
    :param l_atm: atmospheric loss factor (1.0 -> no loss)
    :param l_rx: receiving chain losses
    :param bias_sigma0: sigma0 bias in dB
    :return: sigma0 in dB
    """

    k = ((4.*np.pi)**3. * r**4. * l_atm * l_rx)/(lambda_0**2. * g_0**2. * a)
    return 10. * np.log10(rx_pwr/tx_pwr) + 10. * np.log10(k) + bias_sigma0


def get_sigma0(wf_peak_power_watt: float,
               tx_pwr: float,
               r: float,
               wf_thermal_noise_watt: float = 0.0,
               lambda_0: float = 0.022084,
               band_width: float = 320000000.0,
               g_0: float = 19054.607179632483,
               bias_sigma0: float = 0.0,
               l_atm: float = 1.0,
               l_rx: float = 1.0,
               c_0: float = 299792458.0,
               ) -> float:
    """
    Compute the radar backscatter coefficient sigma nought (sigma0) for lrm waveforms. Uses the radar
    equation and the computation of the pulse limited footprint

    Applicable Documents:

        Guidelines for reverting Waveform Power to Sigma Nought for
        CryoSat-2 in SAR mode (v2.2), Salvatore Dinardo, 23/06/2016
        XCRY-GSEG-EOPS-TN-14-0012

    :param wf_peak_power_watt:  waveform peak power in watt
    :param tx_pwr: transmitted peak power in watt
    :param r: range from satellite center of mass to surface reflection point
        (to be appoximated by satellite altitude if no retracker range available)
    :param wf_thermal_noise_watt: estimate of thermal noise power in watt (default: 0.0)
        will be used to estimate waveform amplitude (Pu)
    :param band_width: CryoSat-2 pulse bandwidth in Hz
    :param lambda_0: radar wavelength in meter (default: 0.022084 m for CryoSat-2 Ku Band altimeter)
    :param g_0: antenna gain at boresight (default: 10^(4.28) from document)
    :param bias_sigma0: sigma nought bias (default: 0.0)
    :param l_atm: two ways atmosphere losses (to be modelled) (default: 1.0 (no loss))
    :param l_rx: receiving chain (RX) waveguide losses (to be characterized) (default: 1.0 (no loss))
    :param c_0: vacuum light speed in meter/sec

    :return sigma_0: Radar backscatter coefficient
    """

    # In the document it is referred to as "waveform power value in output
    # of the re-tracking stage", however generally it is referred to as
    # "waveform amplitude" that is obtained by a waveform function fit
    # It is the scope of this function to provide a sigma0 estimate without
    # proper retracking, therefore Pu is simply defined by the peak power
    # and the thermal noise_power in watt
    pu = wf_peak_power_watt + wf_thermal_noise_watt

    # Intermediate steps & variables
    footprint_radius = np.sqrt(r * c_0 / band_width) / 1000.
    a_lrm = np.pi * footprint_radius ** 2.
    k = ((4.*np.pi)**3. * r**4. * l_atm * l_rx)/(lambda_0**2. * g_0**2. * a_lrm)

    return 10. * np.log10(pu/tx_pwr) + 10. * np.log10(k) + bias_sigma0


class TFMRALeadingEdgeWidth(object):
    """
    Container for computation of leading edge width by taking differences
    between first maximum power thresholds
    """

    def __init__(self, rng, wfm, radar_mode, retrack_flag, tfmra_options=None):
        """
        Compute filtered waveform and index of first maximum once. Calling this class
        will cause the a preliminary retracking of the waveforms indicated by the
        retracker flag. The information is stored in the class and the leading edge
        width can be extraced with the `get_width_from_thresholds` method.

        :param rng: (np.array, dim=(n_records, n_range_bins))
            Waveform range bins
        :param wfm: (np.array, dim=(n_records, n_range_bins))
            Waveform power in arbitrary units
        :param radar_mode: (np.array, dim=(n_records))
            radar mode flag
        :param retrack_flag: (np.array, dim=(n_records))
            flag indicating which waveforms should be retracked
        :param tfmra_options: (doct)
        """

        # Init the cTFRMA retracker
        self.tfmra = cTFMRA()
        self.tfmra.set_default_options(tfmra_options)

        # Pre-process the waveforms for retracking and store the result to self
        filt_rng, filt_wfm, fmi, norm = self.tfmra.get_preprocessed_wfm(rng, wfm, radar_mode, retrack_flag)
        self.wfm, self.rng, self.fmi = filt_wfm, filt_rng, fmi

    def get_width_from_thresholds(self, thres0, thres1, **kwargs):
        """
        Returns the range difference in range bin units between two thresholds,
        by subtracting the range value of thresh0 from thresh1. This is done
        for all waveforms passed to this class during initialization.
        Intended to compute the width of the leading edge.
        :param thres0: (float) The minimum threshold
        :param thres1: (float) The minimum threshold
        :return:
        """
        return self.tfmra.get_thresholds_distance(self.rng, self.wfm, self.fmi, thres0, thres1, **kwargs)


class L1PLeadingEdgeWidth(L1PProcItem):
    """
    A L1P pre-processor item class for computing leading edge width of a waveform
    using the TFMRA retracker as the difference between two thresholds. The
    unit for leading edge width are range bins """

    def __init__(self, **cfg):
        """
        Init the class with the mandatory options

        :param cfg: (dict) Required options (see self.required.options)
        """
        super(L1PLeadingEdgeWidth, self).__init__(**cfg)

        # Init Required Options
        self.tfmra_leading_edge_start = None
        self.tfmra_leading_edge_end = None
        self.tfmra_options = None

        # Get the option settings from the input
        for option_name in self.required_options:
            option_value = cfg.get(option_name)
            if option_value is None:
                msg = f"Missing option `{option_name}` -> No computation of leading edge width!"
                logger.warning(msg)
            setattr(self, option_name, option_value)

    def apply(self, l1: "Level1bData") -> None:
        """
        API class for the Level-1 pre-processor. Functionality is compute leading edge width (full, first half &
        second half) and adding the result to the classifier data group
        :param l1: A Level-1 data instance
        :return: None, Level-1 object is change in place
        """

        # Prepare input
        radar_mode = l1.waveform.radar_mode
        is_valid = l1.surface_type.get_by_name("ocean").flag

        # Compute the leading edge width (requires TFMRA retracking) using
        # -> Use a wrapper for cTFMRA
        width = TFMRALeadingEdgeWidth(l1.waveform.range,
                                      l1.waveform.power,
                                      radar_mode,
                                      is_valid,
                                      tfmra_options=self.tfmra_options)
        lew = width.get_width_from_thresholds(self.tfmra_leading_edge_start,
                                              self.tfmra_leading_edge_end)

        # Add result to classifier group
        l1.classifier.add(lew, "leading_edge_width")

    @property
    def required_options(self):
        return ["tfmra_leading_edge_start",
                "tfmra_leading_edge_end",
                "tfmra_options"]


class L1PSigma0(L1PProcItem):
    """
    A L1P pre-processor item class for computing the backscatter coefficient (sigma0) from
    waveform data

    """

    def __init__(self, **cfg):
        super(L1PSigma0, self).__init__(**cfg)

    def apply(self, l1):
        """
        API class for the Level-1 pre-processor. Functionality is to compute
        leading edge width (full, first half & second half) and adding the result to
        the classifier data group

        :param l1: A Level-1 data instance

        :return: None, Level-1 object is change in place
        """

        # Get input power

        radar_mode = l1.get_parameter_by_name("waveform", "radar_mode")
        rx_power = get_waveforms_peak_power(l1.waveform.power)
        ocog_amplitude = l1.get_parameter_by_name("classifier", "ocog_amplitude")
        is_lrm = radar_mode == 0
        rx_power[is_lrm] = ocog_amplitude[is_lrm]

        # Get Input parameter from l1 object
        tx_power = l1.get_parameter_by_name("classifier", "transmit_power")
        if tx_power is None:
            msg = "classifier `transmit_power` must exist for this pre-processor item -> aborting"
            logger.warning(msg)
            return

        # The computation of sigma0 requires the range to the surface as input
        # and is a step after the retracker. Here we use the satellite altitude
        # as an approximation for sea ice and ocean surfaces.
        altitude = l1.time_orbit.altitude

        # Compute absolute satellite velocity
        sat_vel_x = l1.get_parameter_by_name("classifier", "satellite_velocity_x")
        sat_vel_y = l1.get_parameter_by_name("classifier", "satellite_velocity_y")
        sat_vel_z = l1.get_parameter_by_name("classifier", "satellite_velocity_z")
        if sat_vel_x is None or sat_vel_y is None or sat_vel_z is None:
            # TODO: This is only strictly true for SAR waveforms, observe if necessary
            msg = "classifier `satellite_velocity_[x|y|z]` must exist for this pre-processor item -> aborting"
            logger.warning(msg)
            return
        velocity = np.sqrt(sat_vel_x**2. + sat_vel_y**2. + sat_vel_z**2.)

        # Compute sigma_0
        sigma0 = self.get_sigma0(rx_power, tx_power, altitude, velocity, radar_mode)

        # Add the classifier
        l1.classifier.add(rx_power, "peak_power")
        l1.classifier.add(sigma0, "sigma0")

    def get_sigma0(self,
                   rx_power: np.ndarray,
                   tx_power: np.ndarray,
                   altitude: np.ndarray,
                   velocity: np.ndarray,
                   radar_mode: np.ndarray
                   ) -> np.ndarray:
        """

        :param rx_power:
        :param tx_power:
        :param altitude:
        :param velocity:
        :param radar_mode:
        :return:
        """

        # The function to compute the footprint and the required keywords
        # depend on the radar mode id.
        footprint_func_dict = {0: get_footprint_pulse_limited,
                               1: get_footprint_sar,
                               2: get_footprint_sar}

        footprint_func_kwargs = {0: self.cfg["footprint_pl_kwargs"],
                                 1: self.cfg["footprint_sar_kwargs"],
                                 2: self.cfg["footprint_sar_kwargs"]}

        # The computation of sigma0 depends on properties
        # of the radar altimeter
        sigma0_kwargs = self.cfg.get("sigma0_kwargs", None)

        # Init the output array
        sigma0 = np.full(rx_power.shape, np.nan)

        # Compute sigma0 per waveform
        for i in np.arange(rx_power.shape[0]):

            # Compute the footprint area
            args = (altitude[i], ) if radar_mode[i] == 0 else (altitude[i], velocity[i])
            func = footprint_func_dict[radar_mode[i]]
            footprint_area = func(*args, **footprint_func_kwargs[radar_mode[i]])

            # Compute the backscatter coefficient
            sigma0[i] = get_sigma0_sar(rx_power[i],
                                       tx_power[i],
                                       altitude[i],
                                       footprint_area,
                                       **sigma0_kwargs)

        # Eliminate infinite values
        sigma0[np.isinf(sigma0)] = np.nan

        return sigma0

    @property
    def required_options(self):
        return ["footprint_pl_kwargs", "footprint_sar_kwargs", "sigma0_kwargs", "sigma0_bias"]


class L1PWaveformPeakiness(L1PProcItem):
    """
    A L1P pre-processor item class for computing pulse peakiness """

    def __init__(self,
                 skip_first_range_bins: int = 0,
                 norm_is_range_bin: bool = True
                 ):

        cfg = {"skip_first_range_bins": skip_first_range_bins,
               "norm_is_range_bin": norm_is_range_bin
               }
        super(L1PWaveformPeakiness, self).__init__(**cfg)

    def apply(self, l1: Level1bData) -> None:
        """
        Computes pulse peakiness and adds parameter to classifier data group.

        NOTE: The classifier parameter name depends on the `norm_is_range_bin keyword:

            norm_is_range_bin = True -> parameter name: 'peakiness'
            norm_is_range_bin = False -> parameter name: 'peakiness_normed'

        :param l1: l1bdata.Level1bData instance

        :raises None:

        :return: None
        """
        waveforms = l1.waveform.power
        pulse_peakiness = self.compute_for_waveforms(waveforms)
        parameter_target_name = "peakiness" if self.norm_is_range_bin else "peakiness_normed"
        l1.classifier.add(pulse_peakiness, parameter_target_name)

    def compute_for_waveforms(self, waveforms: npt.NDArray) -> npt.NDArray:
        """
        Compute pulse peakiness for a waveform array

        :param waveforms:

        :return: pulse peakiness array
        """

        # Get the waveform
        n_records, n_range_bins = waveforms.shape
        if waveforms.dtype.kind != "f":
            waveforms = waveforms.astype(np.float64)

        # Get the norm (default is range bins)
        norm = n_range_bins if self.norm_is_range_bin else 1.0

        # Init output parameters
        pulse_peakiness = np.full(n_records, np.nan)

        # Compute peakiness for each waveform
        for i in np.arange(n_records):
            waveform = waveforms[i, self.skip_first_range_bins:]
            pulse_peakiness[i] = self._compute(waveform, norm)

        return pulse_peakiness

    def compute_for_waveform(self, waveform: npt.NDArray) -> float:
        """
        Compute pulse peakiness for a single waveform

        :param waveform:

        :return: pulse peakiness
        """

        # Get the waveform
        n_range_bins = waveform.shape
        if waveform.dtype.kind != "f":
            waveform = waveform.astype(np.float64)
        waveform = waveform[self.skip_first_range_bins:]

        # Get the norm (default is range bins)
        norm = n_range_bins if self.norm_is_range_bins else 1.0

        return self._compute(waveform, norm)

    @staticmethod
    def _compute(waveform: npt.NDArray, norm: Union[int, float]) -> float:
        """
        Compute pulse peakiness for a single waveform

        :param waveform: The waveform
        :param norm

        :return: pulse peakiness
        """
        try:
            pulse_peakiness = bn.nanmax(waveform) / (bn.nansum(waveform)) * norm
        except ZeroDivisionError:
            pulse_peakiness = np.nan
        return pulse_peakiness

    @property
    def required_options(self):
        return ["skip_first_range_bins", "norm_is_range_bin"]


class L1PLeadingEdgeQuality(L1PProcItem):
    """
    Class to compute a leading edge width quality indicator
    Requires `first_maximum_index` classifier parameter
    """

    def __init__(self, **cfg):
        super(L1PLeadingEdgeQuality, self).__init__(**cfg)
        for option_name in self.required_options:
            if option_name not in self.cfg.keys():
                logger.error(f"Missing option: {option_name} -> Leading Edge Quality will not be computed")

    def apply(self, l1):
        """
        Adds a quality indicator for the leading edge

        :param l1: l1bdata.Level1bData instance

        :return: None
        """

        # Get the waveform power
        wfm_power = l1.waveform.power

        # Create output parameters
        leq = np.full(l1.info.n_records, np.nan)              # leading edge quality
        fmi = np.full(l1.info.n_records, -1, dtype=int)       # first maximum index
        fmp = np.full(l1.info.n_records, np.nan)              # first maximum power fraction (to peak power)

        # --- Get the required options ---

        # Waveform window in number of range bins before the first maximum
        leading_edge_lookup_window = self.cfg.get("leading_edge_lookup_window", None)
        if leading_edge_lookup_window is None:
            return

        # The window for the quality computation depends on the radar mode
        window = leading_edge_lookup_window.get(l1.radar_modes, None)
        if window is None:
            logger.error(f"leading_edge_lookup_window not defined for radar mode: {l1.radar_modes}")
            return

        # Normalized power threshold for identifications of the first maximum
        first_maximum_normalized_power_threshold = self.cfg.get("first_maximum_normalized_power_threshold", None)
        if first_maximum_normalized_power_threshold is None:
            return

        # The power threshold depends on the radar mode
        power_threshold = first_maximum_normalized_power_threshold.get(l1.radar_modes, None)
        if window is None:
            logger.error(f"first_maximum_normalized_power_threshold not defined for radar mode: {l1.radar_modes}")
            return

        # Normalized power threshold for identifications of the first maximum
        minimum_valid_first_maximum_index = self.cfg.get("minimum_valid_first_maximum_index", None)
        if minimum_valid_first_maximum_index is None:
            return

        # The power threshold depends on the radar mode
        fmi_min = minimum_valid_first_maximum_index.get(l1.radar_modes, None)
        if window is None:
            logger.error(f"minimum_valid_first_maximum_index not defined for radar mode: {l1.radar_modes}")
            return

        # Loop over all waveforms
        for i in np.arange(l1.info.n_records):

            # Prepare the waveform data
            wfm = wfm_power[i, :].astype(float)    # Get the subset and ensure data type float
            wfm /= np.nanmax(wfm)                  # Normalize (requires float)

            # Get the first maximum index
            fmi_idx = cTFMRA.get_first_maximum_index(wfm, power_threshold)
            if fmi_idx == -1 or fmi_idx < fmi_min:
                continue
            fmi[i] = fmi_idx

            # Save the power values
            fmp[i] = wfm[fmi[i]]

            # Get the search range
            i0, i1 = fmi[i]-window, fmi[i]+1
            i0 = max(i0, 1)
            power_diff = wfm[i0:i1]-wfm[i0-1:i1-1]
            positive_power_diff = power_diff[power_diff > 0]
            total_power_raise = np.sum(positive_power_diff) + wfm[i0]

            # Leading edge quality indicator
            leq[i] = total_power_raise / fmp[i]

        # Add the classifier to the l1 object
        l1.classifier.add(leq, "leading_edge_quality")
        l1.classifier.add(fmi, "first_maximum_index")
        l1.classifier.add(fmp, "first_maximum_power")

    @property
    def required_options(self):
        return ["leading_edge_lookup_window", "first_maximum_normalized_power_threshold",
                "minimum_valid_first_maximum_index"]


class L1PLateTail2PeakPower(L1PProcItem):

    def __init__(self, **cfg) -> None:
        super(L1PLateTail2PeakPower, self).__init__(**cfg)

    def apply(self, l1: Level1bData) -> None:
        """
        Compute late tail to peakiness and add to l1 data object
        :param l1:

        :return:
        """

        # Init the classifier
        ltpp = np.full(l1.info.n_records, np.nan)

        # The power threshold depends on the radar mode
        late_tail_window_idx_dict = self.cfg.get("late_tail_window_idx", None)
        if late_tail_window_idx_dict is None:
            logger.error("settings `late_tail_window_idx` is not defined")
            return

        # Loop over all waveforms and compute ltpp
        wfm = l1.waveform.power
        radar_mode = l1.waveform.radar_mode
        invalid_radar_mode_encountered = 0
        for i in np.arange(wfm.shape[0]):
            late_tail_window_idx = late_tail_window_idx_dict.get(
                RadarModes.get_name(radar_mode[i]), None
            )
            if not late_tail_window_idx:
                ltpp[i] = np.nan
                invalid_radar_mode_encountered += 1
                continue
            ltpp[i] = self.late_tail_to_peak_power(wfm[i, :], late_tail_window_idx)

        if invalid_radar_mode_encountered > 0:
            logger.warning(
                f"Could not compute ltpp for {invalid_radar_mode_encountered} waveforms [invalid radar mode]"
            )
        l1.classifier.add(ltpp, "late_tail_to_peak_power")

    @staticmethod
    def late_tail_to_peak_power(
            wfm: npt.NDArray,
            late_tail_window_idxs: List[int]
    ) -> float:
        """

        :param wfm:
        :param late_tail_window_idxs:
        :return:
        """

        # Get Properties
        wfm_max_power = bn.nanmax(wfm)
        wfm_max_power_index = bn.nanargmax(wfm)
        late_tail_window = np.arange(
            wfm_max_power_index + late_tail_window_idxs[0],
            wfm_max_power_index + late_tail_window_idxs[1] + 1
        )

        # Do not compute ltpp if late tail window is incomplete
        if (late_tail_window >= wfm.shape[0]).any():
            return np.nan

        try:
            return bn.nanmean(wfm[late_tail_window])/wfm_max_power
        except ZeroDivisionError:
            return np.nan


@dataclass
class WFMTrailingEdgeData:
    """ Data Class for L1PTrailingEdgeProperties """

    # Input parameters
    waveform_full: np.ndarray
    trailing_edge_start_idx: int
    waveform_trailing_edge_idx: np.ndarray
    trailing_edge_lower_envelope_idx: np.ndarray

    # Decay Parameters
    trailing_edge_fit: np.ndarray = field(default=None)
    trailing_edge_decay: float = field(default=np.nan)
    trailing_edge_decay_fit_quality: float = field(default=np.nan)
    trailing_edge_decay_mean_absolute_difference: float = field(default=np.nan)

    def decay_debug_plot(self, title="") -> None:

        import matplotlib.pyplot as plt
        x = np.arange(self.waveform_full.size)
        plt.figure(dpi=150)

        # Plot Waveform
        plt.plot(x, self.waveform_full, label="Waveform")

        # Plot Trailing Edge Lower Envelope
        idx_lower_envelope = self.trailing_edge_lower_envelope_idx + self.trailing_edge_start_idx
        plt.scatter(
            idx_lower_envelope,
            self.waveform_full[idx_lower_envelope],
            color="none", edgecolors="red", label="Trailing Edge Lower Envelope"
        )

        # Plot Fit
        if self.trailing_edge_fit is not None:
            plt.plot(
                self.waveform_trailing_edge_idx,
                self.trailing_edge_fit
            )
        plt.legend()
        plt.xlabel("Range Gate")
        plt.ylabel("Waveform Power Units")
        plt.title(title)
        plt.show()

    @property
    def trailing_edge_size(self) -> int:
        return self.waveform_trailing_edge_idx.size

    @property
    def waveform_trailing_edge_subset(self) -> np.ndarray:
        return self.waveform_full[self.waveform_trailing_edge_idx]

    @property
    def waveform_trailing_edge_subset_lower_envelope(self) -> np.ndarray:
        return self.waveform_trailing_edge_subset[self.trailing_edge_lower_envelope_idx]

    @property
    def first_maximum_power(self) -> float:
        return self.waveform_full[self.trailing_edge_start_idx]

    @property
    def decay_parameters(self) -> Tuple[float, float, float]:
        return (
            self.trailing_edge_decay,
            self.trailing_edge_decay_fit_quality,
            self.trailing_edge_decay_mean_absolute_difference
        )


class L1PTrailingEdgeProperties(L1PProcItem):

    def __init__(self, **cfg) -> None:
        super(L1PTrailingEdgeProperties, self).__init__(**cfg)
        # TODO: To be implemented as a Pydantic v2 data model
        self.config_validation_schema.validate(cfg)

    def apply(self, l1: Level1bData) -> None:
        """
        Compute late tail to peakiness and add to l1 data object

        :param l1: The Level-1 data container

        :return:
        """

        # Get configuration
        valid_first_maximum_index_range = self.cfg.get("valid_first_maximum_index_range")
        trailing_edge_width_kwargs = self.cfg.get("trailing_edge_width")
        use_multiprocessing = self.cfg.get("use_multiprocessing", False)

        timer = StopWatch().start()

        # Loop over all waveforms and compute parameters
        wfm = l1.waveform.power
        result = WaveFormTrailingEdgeParameter(
            wfm,
            valid_first_maximum_index_range,
            trailing_edge_width_kwargs
        )

        timer.stop()
        logger.debug(f"- Trailing edge properties computed in {timer.get_seconds():.3f} seconds")
        logger.debug(f"- Trailing edge decay fit as failed for {result.params.decay_fit_has_failed.sum()} waveforms")

        l1.classifier.add(result.params.ted, "trailing_edge_decay")
        l1.classifier.add(result.params.tedfq, "trailing_edge_decay_fit_quality")
        l1.classifier.add(result.params.temad, "trailing_edge_mean_absolute_deviation")
        l1.classifier.add(result.params.tew, "trailing_edge_width")
        l1.classifier.add(result.params.teq, "trailing_edge_quality")

    @staticmethod
    def _compute_single_processing(*args) -> "WaveFormTrailingEdgeParameterData":
        """
        Compute waveform trailing edge parameter on a single thread

        :param args: Arguments to WaveFormTrailingEdgeParameter

        :return: Waveform trailing edge parameter object
        """
        return WaveFormTrailingEdgeParameter(*args).params

    # @staticmethod
    # def _compute_multi_processing(
    #         waveforms: npt.NDArray,
    #         valid_first_maximum_index_range: List,
    #         trailing_edge_width_kwargs: Dict
    # ) -> "WaveFormTrailingEdgeParameterData":
    #     """
    #     Compute waveform trailing edge parameter on a single thread
    #
    #     :param args: Arguments to WaveFormTrailingEdgeParameter
    #
    #     :return: Waveform trailing edge parameter object
    #     """
    #
    #     chunk_index_list, n_processes = get_1d_array_multiprocessing_chunks(waveforms.shape[0])
    #     logger.info(f"- Use multi-processing with {n_processes} workers")
    #     logger.disable("pysiral")
    #     with multiprocessing.pool.ThreadPool(n_processes) as process_pool:
    #         result = process_pool.map(
    #             WaveFormTrailingEdgeParameter,
    #             args=(
    #                 waveforms,
    #                 valid_first_maximum_index_range,
    #                 trailing_edge_width_kwargs
    #             )
    #         )
    #     # process_pool.close()
    #     # process_pool.join()
    #     # results = [job.get() for job in jobs]
    #     logger.enable("pysiral")
    #     #
    #     # params = results[0].params
    #     # for result in results[1:]:
    #     #     params.add(result.params)
    #     return params

    @property
    def config_validation_schema(self) -> Schema:

        return Schema(
            {
                "use_multiprocessing": bool,
                "valid_first_maximum_index_range": And(
                    list,
                    lambda x: all(item > 0 for item in x)
                ),
                "trailing_edge_width": {
                    "trailing_edge_end_power_treshold_normed": And(
                        float,
                        lambda x: 0.0 <= x <= 1.0,
                    ),
                    "oversample_factor": And(
                        int,
                        lambda x: x > 0
                    )
                },
            },
            error=f"Invalid config for {self.__class__.__name__}:\n "+"{}"
        )


class WaveFormTrailingEdgeParameterData(object):

    def __init__(self, n_records: int) -> None:
        """

        :param n_records:
        """
        self.tew = np.full(n_records, np.nan)  # Trailing Edge Width
        self.teq = np.full(n_records, np.nan)  # Trailing Edge Quality
        self.ted = np.full(n_records, np.nan)  # Trailing Edge Decay (Müller et al. 2023)
        self.tedfq = np.full(n_records, np.nan)  # Trailing Edge Decay Fit Quality (Müller et al. 2023)
        self.temad = np.full(n_records, np.nan)
        self.decay_fit_has_failed = np.full(n_records, True)

    def add(self, other_params: "WaveFormTrailingEdgeParameterData") -> None:

        self.tew = np.append(self.tew, other_params.tew)
        self.teq = np.append(self.teq, other_params.teq)
        self.ted = np.append(self.ted, other_params.ted)
        self.tedfq = np.append(self.tedfq, other_params.tedfq)
        self.temad = np.append(self.temad, other_params.temad)
        self.decay_fit_has_failed = np.append(self.decay_fit_has_failed, other_params.decay_fit_has_failed)


class WaveFormTrailingEdgeDecayFit(object):

    def __init__(self, fit_kwargs=None):
        self.default_fit_kwargs = {
            "method": "trf",
            "ftol": 1e-3,
            "gtol": 1e-3,
            "xtol": 1e-4,
            "loss": "linear",
            "bounds": ([0.0, 0.0, -np.inf], [np.inf, np.inf, np.inf])
        }
        self.fit_kwargs = self.default_fit_kwargs
        self.p0 = None
        if isinstance(fit_kwargs, dict):
            self.fit_kwargs.update(fit_kwargs)

    def fit(self, data):
        p0 = [np.float64(data.first_maximum_power), np.float64(1e-01), np.float64(0.0)]
        try:
            popt, pcov, *_ = curve_fit(
                inverse_power,
                data.trailing_edge_lower_envelope_idx.astype(np.float64) + 1.0,
                data.waveform_trailing_edge_subset_lower_envelope,
                p0=p0,
                **self.fit_kwargs
            )
            return popt
        # Catch fit errors
        except (TypeError, RuntimeError, ValueError):
            return None

    def fit_chunks(self, chunks: List) -> List:
        """

        :param chunks:
        :return:
        """
        results = []
        for data in chunks:
            try:
                popt, pcov, *_ = curve_fit(
                    inverse_power,
                    data.trailing_edge_lower_envelope_idx.astype(np.float64) + 1.0,
                    data.waveform_trailing_edge_subset_lower_envelope,
                    p0=self.p0,
                    **self.fit_kwargs
                )
                self.p0 = popt
            # Catch fit errors
            except (TypeError, RuntimeError, ValueError):
                popt = (np.nan, np.nan, np.nan)
            results.append(popt)
        return results


class WaveFormTrailingEdgeParameter(object):

    def __init__(self,
                 waveforms: npt.NDArray,
                 valid_first_maximum_index_range: List,
                 trailing_edge_width_kwargs: Dict
                 ):
        """
        Algorithm class to compute trailing edge parameter

        :param waveforms: An array of waveforms.
        :param valid_first_maximum_index_range: F
        :param trailing_edge_width_kwargs:

        :raises None:
        """

        # Input data and config
        self.waveforms = waveforms
        self.valid_first_maximum_index_range = valid_first_maximum_index_range
        self.trailing_edge_width_kwargs = trailing_edge_width_kwargs

        # Fit properties
        self.last_fit_popt = None

        # Result property
        self.params = WaveFormTrailingEdgeParameterData(waveforms.shape[0])

        # Compute the parameter
        self._compute()

    def _compute(self) -> None:
        """
        Loop over all waveforms and compute the parameters
        # TODO: Need documentation and potentially move to sub-modules
        """

        waveform_index_list = []
        waveform_data_stack = []
        num_waveforms = self.waveforms.shape[0]
        for i in np.arange(num_waveforms):

            # Collect data and sanity check
            trailing_edge_data = self.get_waveform_trailing_edge_data(self.waveforms[i, :])
            if (
                    self.valid_first_maximum_index_range[0] > trailing_edge_data.trailing_edge_start_idx or
                    self.valid_first_maximum_index_range[1] < trailing_edge_data.trailing_edge_start_idx
            ):
                continue
            waveform_data_stack.append(trailing_edge_data)
            waveform_index_list.append(i)

        # multiprocessing fit
        fit_class = WaveFormTrailingEdgeDecayFit()
        chunks_idxs, cpu_count = get_multiprocessing_1d_array_chunks(num_waveforms, psrlcfg.CPU_COUNT)
        trailing_edge_data_chunks = [
            waveform_data_stack[chunks_idx[0]:chunks_idx[1] + 1]
            for chunks_idx in chunks_idxs
        ]
        pool = multiprocessing.Pool(cpu_count)
        result_chunks = pool.imap(fit_class.fit_chunks, trailing_edge_data_chunks)
        pool.close()
        pool.join()

        results = [result for results in result_chunks for result in results]

        # Note: Not all waveforms are fitted, e.g. if the trailing edge starts to
        #       late or too early. The list of results therefore might be shorter
        #       than the list of waveforms. `waveform_index_list` therefore maps
        #       the fit index to the original waveform index.
        for fit_index, i in enumerate(waveform_index_list):

            trailing_edge_data = self.get_trailing_edge_decay(waveform_data_stack[fit_index], results[fit_index])

            # In this case the fitting algorithm has failed
            if trailing_edge_data.trailing_edge_fit is None:
                # trailing_edge_data.decay_debug_plot("Trailing Edge fit failed")
                continue
            self.params.ted[i], self.params.tedfq[i], self.params.temad[i] = trailing_edge_data.decay_parameters
            self.params.decay_fit_has_failed[i] = False

            # Compute other trailing edge properties
            self.params.tew[i] = self.get_trailing_edge_width(trailing_edge_data, **self.trailing_edge_width_kwargs)
            self.params.teq[i] = self.get_trailing_edge_quality(trailing_edge_data)

    def get_waveform_trailing_edge_data(
            self,
            waveform_full: np.ndarray,
    ) -> WFMTrailingEdgeData:
        """
        Aggregate all necesary trailing edge properties in a data class.

        :param waveform_full: The full waveform in original units

        :return: Data class
        """

        trailing_edge_start_idx = np.argmax(waveform_full)
        trailing_edge_idx = np.arange(trailing_edge_start_idx, waveform_full.size)

        # Get indices of trailing edge lower envelope
        # Only the lower envelope will be used to compute offset/scale
        trailing_edge_lower_envelope_idx = self.get_trailing_edge_lower_envelope(
            waveform_full[trailing_edge_idx]
        )

        return WFMTrailingEdgeData(
            waveform_full,
            trailing_edge_start_idx,
            trailing_edge_idx,
            trailing_edge_lower_envelope_idx,
        )

    def get_trailing_edge_decay(self, data: WFMTrailingEdgeData, popt: List) -> WFMTrailingEdgeData:
        """
        Fit inverse power law to the lower envelope of the trailing edge
        and compute decay factor, fit quality and residual to trailing
        edge power. Result will be added to data class.

        :param data: Pre-compiled trailing edge data
        :param popt: Fit coefficient (function inverse power low)

        :raises None:

        :return: Data class with updated results
        """

        # # Execute the fit
        # # NOTE: x-values shall not be 0 (0**anything -> ZeroDivisionError)
        # p0 = \
        #     [np.float64(data.first_maximum_power), np.float64(1e-01), np.float64(0.0)] \
        #     if self.last_fit_popt is None else self.last_fit_popt
        #
        # try:
        #     popt, pcov, *_ = curve_fit(
        #         inverse_power,
        #         data.trailing_edge_lower_envelope_idx.astype(np.float64) + 1.0,
        #         data.waveform_trailing_edge_subset_lower_envelope,
        #         p0=p0,
        #         **self.fit_kwargs
        #     )
        # # Catch fit errors
        # except (TypeError, RuntimeError):
        #     return data
        # self.last_fit_popt = popt

        # Compute fit
        x_values = np.arange(data.trailing_edge_size, dtype=np.float64) + 1
        data.trailing_edge_fit = inverse_power(x_values, *popt)

        # Compute parameters
        data.trailing_edge_decay = popt[1]  # second parameter in inverse power function
        data.trailing_edge_decay_mean_absolute_difference = bn.nanmedian(
            np.abs(data.trailing_edge_fit - data.waveform_trailing_edge_subset)
        )
        data.trailing_edge_decay_fit_quality = coeficient_of_determination(
            data.waveform_trailing_edge_subset,
            data.trailing_edge_fit
        )

        return data

    @staticmethod
    def get_trailing_edge_width(
            data: WFMTrailingEdgeData,
            oversample_factor: int,
            trailing_edge_end_power_treshold_normed: float
    ) -> float:
        """
        Compute the width of the trailing edge based on the decay fit.

        :param data:
        :param oversample_factor:
        :param trailing_edge_end_power_treshold_normed:

        :raises None:

        :return: The trailing edge width in range gate units
        """

        # Get oversampled trailing edge fit
        range_gates_oversampled = np.linspace(
            0,
            data.trailing_edge_size - 1,
            data.trailing_edge_size * oversample_factor + 1,
        ) + data.waveform_trailing_edge_idx[0]
        range_gates_oversampled = np.around(range_gates_oversampled, 1)
        oversample_func = interp1d(data.waveform_trailing_edge_idx, data.trailing_edge_fit)
        trailing_edge_oversampled = oversample_func(range_gates_oversampled)

        # Get first (oversampled) range gate, where trailing edge is below
        # normed power threshold defining end of trailing edge
        trailing_edge_end_idx = np.where(
            trailing_edge_oversampled / data.first_maximum_power <= trailing_edge_end_power_treshold_normed
        )[0]

        if trailing_edge_end_idx.size == 0:
            # data.decay_debug_plot("Trailing Edge does not fall below minimum power threshold")
            return np.nan
        trailing_edge_end_range_gate = range_gates_oversampled[trailing_edge_end_idx][0]

        return trailing_edge_end_range_gate - data.trailing_edge_start_idx

    @staticmethod
    def get_trailing_edge_quality(data: WFMTrailingEdgeData) -> float:
        """

        :param data:
        :return:
        """

        power_diff = data.waveform_trailing_edge_subset[1:] - data.waveform_trailing_edge_subset[:-1]
        absolute_negative_power_diff = np.abs(power_diff[power_diff < 0.0])
        total_power_drop = bn.nansum(absolute_negative_power_diff)

        # Trailing edge quality indicator
        return total_power_drop / data.first_maximum_power

    @staticmethod
    def get_trailing_edge_lower_envelope(
            wfm_trailing_edge: npt.NDArray
    ) -> npt.NDArray:
        """

        :param wfm_trailing_edge:
        :return:
        """

        idx_le = argrelmin(wfm_trailing_edge)[0]
        idx_le = np.insert(idx_le, 0, 0)
        idx_le = np.insert(idx_le, len(idx_le), wfm_trailing_edge.size-1)

        wfm_change = [
            wfm_trailing_edge[idx_le[1:]] - wfm_trailing_edge[idx_le[:-1]]
        ][0]
        wfm_change = np.insert(wfm_change, 0, -1)
        return idx_le[wfm_change < 0.0]


class L1PLeadingEdgePeakiness(L1PProcItem):

    def __init__(self, **cfg):
        super(L1PLeadingEdgePeakiness, self).__init__(**cfg)
        for option_name in self.required_options:
            if option_name not in self.cfg.keys():
                logger.error(f"Missing option: {option_name} -> Leading Edge Quality will not be computed")

    def apply(self, l1: Level1bData):
        """
        Mandatory class of a L1 preproceessor item. Computes the leading edge peakiness

        :param l1:
        :return:
        """

        # Init the classifier (leading edge peakiness
        lep = np.full(l1.info.n_records, np.nan)

        # Get the waveform power
        wfm = l1.waveform.power

        # Get the first maximum index
        fmi = l1.classifier.get_parameter("first_maximum_index", raise_on_error=False)
        if fmi is None:
            logger.error("Classifier `first_maximum_index` not available -> skipping leading edge peakiness")
            l1.classifier.add(lep, "leading_edge_peakiness")
            return

        # Get the window for the pulse peakiness computation
        window_size = self.cfg.get("window_size", None)
        if window_size is None:
            logger.error("Option `window size` not available -> skipping leading edge peakiness")
            l1.classifier.add(lep, "leading_edge_peakiness")
            return

        # Loop over all waveforms
        for i in np.arange(wfm.shape[0]):
            if fmi[i] < 0:
                continue
            lep[i] = self.leading_edge_peakiness(wfm[i, :], fmi[i], window_size)

        # Convert inf to nan
        lep[np.isinf(lep)] = np.nan

        # Add to classifier
        l1.classifier.add(lep, "leading_edge_peakiness")

    @staticmethod
    def leading_edge_peakiness(wfm: np.ndarray, fmi: int, window: int) -> float:
        """
        Compute the leading edge peakiness
        :param wfm: Waveform power
        :param fmi: first maximum index
        :param window: the number or leading range bins to the first maximum for the
            peakiness computation
        :return:
        """

        # Get the waveform subset prior to first maximum
        i0 = fmi - window
        i0 = max(i0, 0)
        return wfm[fmi] / bn.nanmean(wfm[i0:fmi]) * (fmi - i0)

    @property
    def required_options(self):
        return ["window_size"]


class CS2OCOGParameter(object):
    """
    Calculate OCOG Parameters (Amplitude, Width) for CryoSat-2 waveform
    counts.
    Algorithm Source: retrack_ocog.pro from CS2AWI lib
    """

    def __init__(self, wfm_counts):
        self._n = np.shape(wfm_counts)[0]
        self._amplitude = np.ndarray(shape=[self._n], dtype=np.float32)
        self._width = np.ndarray(shape=[self._n], dtype=np.float32)
        self._calc_parameters(wfm_counts)

    def _calc_parameters(self, wfm_counts):
        for i in np.arange(self._n):
            y = wfm_counts[i, :].flatten().astype(np.float64)
            y -= bn.nanmean(y[:11])  # Remove Noise
            y[np.where(y < 0.0)[0]] = 0.0  # Set negative counts to zero
            y2 = y ** 2.0
            self._amplitude[i] = np.sqrt((y2 ** 2.0).sum() / y2.sum())
            self._width[i] = ((y2.sum()) ** 2.0) / (y2 ** 2.0).sum()

    @property
    def amplitude(self):
        return self._amplitude

    @property
    def width(self):
        return self._width


class CS2PulsePeakiness(object):
    """
    Calculates Pulse Peakiness (full, left & right) for CryoSat-2 waveform
    counts
    XXX: This is a 1 to 1 legacy implementation of the IDL CS2AWI method,
         consistent method of L1bData or L2Data is required
    """

    def __init__(self, wfm_counts, pad=2):
        shape = np.shape(wfm_counts)
        self._n = shape[0]
        self._n_range_bins = shape[1]
        self._pad = pad
        self._peakiness = np.full(self._n, np.nan).astype(np.float32)
        self._peakiness_r = np.full(self._n, np.nan).astype(np.float32)
        self._peakiness_l = np.full(self._n, np.nan).astype(np.float32)
        self._noise_floor = np.full(self._n, np.nan).astype(np.float32)
        # self.peakiness_no_noise_removal = np.full(self._n, np.nan).astype(np.float32)
        self._peakiness_normed = np.full(self._n, np.nan).astype(np.float32)
        self._calc_parameters(wfm_counts)

    def _calc_parameters(self, wfm_counts):
        for i in np.arange(self._n):
            try:
                y = wfm_counts[i, :].flatten().astype(np.float32)
                self._noise_floor[i] = bn.nanmean(y[:11])
                y_no_noise_removal = y.copy()
                y -= self._noise_floor[i]  # Remove Noise
                y[np.where(y < 0.0)[0]] = 0.0  # Set negative counts to zero
                yp = np.nanmax(y)  # Waveform peak value
                ypi = np.nanargmax(y)  # Waveform peak index
                if 3 * self._pad < ypi < self._n_range_bins - 4 * self._pad:
                    self._peakiness_l[i] = yp / bn.nanmean(y[ypi - 3 * self._pad:ypi - 1 * self._pad + 1]) * 3.0
                    self._peakiness_r[i] = yp / bn.nanmean(y[ypi + 1 * self._pad:ypi + 3 * self._pad + 1]) * 3.0
                    self._peakiness_normed[i] = yp / y.sum()
                    self._peakiness[i] = self.peakiness_normed[i] * self._n_range_bins

                    # self.peakiness_no_noise_removal[i] = np.nanmax(y_no_noise_removal) / y_no_noise_removal.sum() * self._n_range_bins
            except ValueError:
                self._peakiness_l[i] = np.nan
                self._peakiness_r[i] = np.nan
                self._peakiness[i] = np.nan

    @property
    def peakiness(self):
        return self._peakiness

    @property
    def peakiness_normed(self):
        return self._peakiness_normed

    @property
    def noise_floor(self):
        return self._noise_floor

    @property
    def peakiness_r(self):
        return self._peakiness_r

    @property
    def peakiness_l(self):
        return self._peakiness_l


# Late tail to peak power (LTPP) ratio added for fmi needs
class S3LTPP(object):
    """
    Calculates Late-Tail-to-Peak-Power ratio.
    source: Rinne 2016
    """

    def __init__(self, wfm_counts, pad=1):
        # Warning: if 0padding is introduced in S3 L1 processing baseline, pad must be set to 2
        shape = np.shape(wfm_counts)
        self._n = shape[0]
        self._n_range_bins = shape[1]
        self._pad = pad
        dtype = np.float32
        self._ltpp = np.ndarray(shape=[self._n], dtype=dtype) * np.nan
        self._calc_parameters(wfm_counts)

    def _calc_parameters(self, wfm_counts):
        # loop over the waveforms
        for i in np.arange(self._n):
            try:
                y = wfm_counts[i, :].flatten().astype(np.float32)
                y -= bn.nanmean(y[:11])  # Remove Noise
                y[np.where(y < 0.0)[0]] = 0.0  # Set negative counts to zero
                yp = np.nanmax(y)  # Waveform peak value

                if np.isnan(yp):  # if the current wf is nan
                    # no ltpp can be computed
                    self._ltpp[i] = np.nan
                else:
                    ypi = np.nanargmax(y)  # Waveform peak index

                    # gates to compute the late tail:
                    # [ypi+50:ypi+70] if 0padding=2, [ypi+25:ypi+35] if 0padding=1
                    gate_start = ypi + self._pad * 25
                    gate_stop = ypi + self._pad * 35 + 1

                    if gate_start > self._n_range_bins or gate_stop > self._n_range_bins:
                        # not enough gates to compute the LTPP
                        self._ltpp[i] = np.nan
                    else:
                        self._ltpp[i] = np.mean(y[gate_start:gate_stop]) / yp

            except ValueError:
                self._ltpp[i] = np.nan

    @property
    def ltpp(self):
        return self._ltpp


class CS2LTPP(object):
    """
    Calculates Late-Tail-to-Peak-Power ratio.
    """

    def __init__(self, wfm_counts, pad=2):
        shape = np.shape(wfm_counts)
        self._n = shape[0]
        self._n_range_bins = shape[1]
        self._pad = pad
        dtype = np.float32
        self._ltpp = np.ndarray(shape=[self._n], dtype=dtype) * np.nan
        self._calc_parameters(wfm_counts)

    def _calc_parameters(self, wfm_counts):
        # loop over the waveforms
        for i in np.arange(self._n):
            y = wfm_counts[i, :].flatten().astype(np.float32)
            y -= bn.nanmean(y[:11])  # Remove Noise
            y[np.where(y < 0.0)[0]] = 0.0  # Set negative counts to zero
            yp = np.nanmax(y)  # Waveform peak value
            ypi = bn.nanargmax(y)  # Waveform peak index

            # AMANDINE: implementation for wf of 256 bins (128 zero-padded)?
            onediv = float(1) / float(41)

            # AMANDINE: here i seems to be understood as gate index, but it is wf index!?
            if i == 256:
                break
            if [i > (ypi + 100)] and [i < (ypi + 140)]:  # AMANDINE: syntax to be checked
                try:
                    self._ltpp[i] = (onediv * float(y[i])) / float(yp)  # AMANDINE: where is the sum in this formula?
                except ZeroDivisionError:
                    self._ltpp[i] = np.nan

    @property
    def ltpp(self):
        return self._ltpp


class EnvisatWaveformParameter(object):
    """
    Currently only computes pulse peakiness for Envisat waveforms
    from SICCI processor.

    Parameter for Envisat from SICCI Processor
        skip = 5
        bins_after_nominal_tracking_bin = 83
    """

    def __init__(self, wfm, skip=5, bins_after_nominal_tracking_bin=83):
        self.t_n = bins_after_nominal_tracking_bin
        self.skip = skip
        self._n = wfm.shape[0]
        self._n_range_bins = wfm.shape[1]
        self._init_parameter()
        self._calc_parameter(wfm)

    def _init_parameter(self):
        self.peakiness_old = np.ndarray(shape=self._n, dtype=np.float32)
        self.peakiness = np.ndarray(shape=self._n, dtype=np.float32) * np.nan

    def _calc_parameter(self, wfm):

        for i in np.arange(self._n):
            # Discard first bins, they are FFT artefacts anyway
            wave = wfm[i, self.skip:]

            # old peakiness
            try:
                pp = 0.0 + self.t_n * float(max(wave)) / float(sum(wave))
            except ZeroDivisionError:
                pp = np.nan
            self.peakiness_old[i] = pp

            # new peakiness
            try:
                self.peakiness[i] = float(max(wave)) / float(sum(wave)) * self._n_range_bins
            except ZeroDivisionError:
                self.peakiness[i] = np.nan


def coeficient_of_determination(y: np.ndarray, y_fit: np.ndarray) -> float:
    """
    Compute coeficient of determination

    :param y:
    :param y_fit:

    :return:
    """
    # residual sum of squares
    ss_res = bn.nansum((y - y_fit) ** 2)

    # total sum of squares
    ss_tot = bn.nansum((y - bn.nanmean(y)) ** 2)

    return 1 - (ss_res / ss_tot)



