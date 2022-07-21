# -*- coding: utf-8 -*-
"""
Created on Fri Jul 01 13:07:10 2016

@author: shendric
"""

import numpy as np
import numpy.typing as npt
import bottleneck as bn
from loguru import logger

from pysiral.l1bdata import Level1bData
from pysiral.retracker import cTFMRA
from pysiral.l1preproc.procitems import L1PProcItem


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

    # Final computation of sigma nought
    sigma0 = 10. * np.log10(pu/tx_pwr) + 10. * np.log10(k) + bias_sigma0

    return sigma0


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

    def get_width_from_thresholds(self, thres0, thres1):
        """
        Returns the range difference in range bin units between two thresholds,
        by subtracting the range value of thresh0 from thresh1. This is done
        for all waveforms passed to this class during initialization.
        Intended to compute the width of the leading edge.
        :param thres0: (float) The minimum threshold
        :param thres1: (float) The minimum threshold
        :return:
        """
        width = self.tfmra.get_thresholds_distance(self.rng, self.wfm, self.fmi, thres0, thres1)
        return width


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
            option_value = cfg.get(option_name, None)
            if option_value is None:
                msg = "Missing option `%s` -> No computation of leading edge width!" % option_name
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
    A L1P pre-processor item class for computing leading edge width (full, first half, second half)
    using three TFMRA thresholds """

    def __init__(self, **cfg):
        super(L1PWaveformPeakiness, self).__init__(**cfg)
        for option_name in self.required_options:
            option_value = cfg.get(option_name, None)
            if option_value is None:
                msg = "Missing option `%s` -> No computation of peakiness!" % option_name
                logger.warning(msg)
            setattr(self, option_name, option_value)

        # Init Parameters
        self.peakiness = None

    def apply(self, l1):
        """
        Computes pulse peakiness for lrm waveforms (from SICCI v1 processor).
        :param l1: l1bdata.Level1bData instance
        :return: None
        """
        self._calc(l1)
        l1.classifier.add(self.peakiness, "peakiness")

    def _calc(self, l1):
        """ Compute pulse peakiness (from SICCI v1 processor)."""

        # Get the waveform
        wfm = l1.waveform.power
        n_records, n_range_bins = wfm.shape

        # Init output parameters
        self.peakiness = np.full(n_records, np.nan)
        self.peakiness_old = np.full(n_records, np.nan)

        # Compute peakiness for each waveform
        for i in np.arange(n_records):

            # Discard first bins, they are FFT artefacts anyway
            wave = wfm[i, self.skip_first_range_bins:]

            # new peakiness
            try:
                self.peakiness[i] = float(max(wave))/float(sum(wave))*n_range_bins
            except ZeroDivisionError:
                self.peakiness[i] = np.nan

    @property
    def required_options(self):
        return ["skip_first_range_bins"]


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
            i0 = i0 if i0 > 1 else 1   # ensure the lower index stays in the valid range
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
        i0 = i0 if i0 > 0 else 0
        i1 = fmi

        # Compute the peakiness
        lep = wfm[fmi] / bn.nanmean(wfm[i0:i1]) * (i1 - i0)

        return lep

    @property
    def required_options(self):
        return ["window_size"]
