# -*- coding: utf-8 -*-
"""
Created on Fri Jul 01 13:07:10 2016

@author: shendric
"""
import matplotlib

import numpy as np
import bottleneck as bn
from loguru import logger

from pysiral.l1bdata import Level1bData
from pysiral.retracker import cTFMRA
from pysiral.logging import DefaultLoggingClass


def get_waveforms_peak_power(wfm, dB=False):
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
    if dB:
        peak_power = 10 * np.log10(peak_power)
    return peak_power


def get_sar_sigma0(wf_peak_power_watt, tx_pwr, r, v_s, **sigma0_par_dict):
    """ Wrapper function to compute sigma nought for all waveforms """
    n_records = wf_peak_power_watt.shape[0]
    sigma0 = np.ndarray(shape=(n_records))
    for i in np.arange(n_records):
        sigma0[i] = sar_sigma0(wf_peak_power_watt[i], tx_pwr[i], r[i], v_s[i], **sigma0_par_dict)
    return sigma0


def sar_sigma0(wf_peak_power_watt, tx_pwr, r, v_s,
               wf_thermal_noise_watt=0.0, ptr_width=2.819e-09, tau_b=0.00352,
               lambda_0=0.022084, wf=1., g_0=19054.607179632483,
               bias_sigma0=0.0, l_atm=1.0, l_rx=1.0,
               c_0=299792458.0, r_mean=6371000.0):
    """
    returns sigma nought (sigma0) for sar waveforms.

    Applicable Documents
    --------------------

    Guidelines for reverting Waveform Power to Sigma Nought for
    CryoSat-2 in SAR mode (v2.2), Salvatore Dinardo, 23/06/2016

    XCRY-GSEG-EOPS-TN-14-0012


    Arguments
    ---------
        wf_peak_power_watt (float)
            waveform peak power in watt
        tx_pwr (float)
            transmitted peak power in watt
        r (float)
            range from satellite center of mass to surface reflection point
            (to be appoximated by satellite altitude if no retracker range
            available)
        v_s (float)
            satellite along track velocity in meter/sec

    Keywords
    --------
        wf_thermal_noise_watt (float)
            estimate of thermal noise power in watt (default: 0.0)
            will be used to estimate waveform amplitude (Pu)
        ptr_width (float)
            3dB range point target response temporal width in seconds
            (default: 2.819e-09 sec for CryoSat-2 SAR)
        tau_b (float)
            burst length in seconds
            (default: 0.00352 sec for CryoSat-2 SAR)
        lambda_0 (float)
            radar wavelength in meter
            (default: 0.022084 m for CryoSat-2 Ku Band altimeter)
        wf (float)
            footprint widening factor (1.486 * rv in case of Hamming window
            application on burst data; rv: unspecified empirical factor)
            (default: 1 no weighting window application)
        g_0 (float)
            antenna gain at boresight
            (default: 10^(4.28) from document)
        bias_sigma0 (float)
            sigma nought bias
            (default: 0.0)
        l_atm (float)
            two ways atmosphere losses (to be modelled)
            (default: 1.0 (no loss))
        l_rx (float)
            receiving chain (RX) waveguide losses (to be characterized)
            (default: 1.0 (no loss))
        c_0 (float)
            vacuum light speed in meter/sec
        r_mean (float)
            mean earth radius in meter

    Returns
    -------
        sigma_0 (float)
        :param wf_peak_power_watt:

    """

    # XXX: The definition of the variable Pu is not quite clear to me (Stefan)
    # In the document it is referred to as "waveform power value in output
    # of the re-tracking stage", however generally it is referred to as
    # "waveform amplitude" that is obtained by a waveform function fit
    # It is the scope of this function to provide a sigma0 estimate without
    # proper retracking, therefore Pu is simply defined by the peak power
    # and the thermal noise_power in watt
    pu = wf_peak_power_watt + wf_thermal_noise_watt

    # Intermediate steps & variables
    pi = np.pi
    alpha_earth = 1. + (r/r_mean)
    lx = (lambda_0 * r)/(2. * v_s * tau_b)
    ly = np.sqrt((c_0 * r * ptr_width)/alpha_earth)
    a_sar = (2. * ly) * (wf * lx)
    k = ((4.*pi)**3. * r**4. * l_atm * l_rx)/(lambda_0**2. * g_0**2. * a_sar)

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


class L1PLeadingEdgeWidth(DefaultLoggingClass):
    """
    A L1P pre-processor item class for computing leading edge width of a waveform
    using the TFMRA retracker as the difference between two thresholds. The
    unit for leading edge width are range bins """

    def __init__(self, **cfg):
        """
        Init the class with the mandatory options
        :param cfg: (dict) Required options (see self.required.options)Ä
        """
        super(L1PLeadingEdgeWidth, self).__init__(self.__class__.__name__)

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

    def apply(self, l1):
        """
        API class for the Level-1 pre-processor. Functionality is compute leading edge width (full, first half &
        second half) and adding the result to the classifier data group
        :param l1: A Level-1 data instance
        :return: None, Level-1 object is change in place
        """

        # Prepare input
        radar_mode = l1.waveform.radar_mode
        is_ocean = l1.surface_type.get_by_name("ocean").flag

        # Compute the leading edge width (requires TFMRA retracking)
        width = TFMRALeadingEdgeWidth(l1.waveform.range, l1.waveform.power, radar_mode, is_ocean,
                                      tfmra_options=self.tfmra_options)
        lew = width.get_width_from_thresholds(self.tfmra_leading_edge_start, self.tfmra_leading_edge_end)

        # Add result to classifier group
        l1.classifier.add(lew, "leading_edge_width")

    @property
    def required_options(self):
        return ["tfmra_leading_edge_start",
                "tfmra_leading_edge_end",
                "tfmra_options"]


class L1PSigma0(DefaultLoggingClass):
    """
    A L1P pre-processor item class for computing leading edge width (full, first half, second half)
    using three TFMRA thresholds """

    def __init__(self, **cfg):
        super(L1PSigma0, self).__init__(self.__class__.__name__)

    def apply(self, l1):
        """
        API class for the Level-1 pre-processor. Functionality is compute leading edge width (full, first half &
        second half) and adding the result to the classifier data group
        :param l1: A Level-1 data instance
        :return: None, Level-1 object is change in place
        """

        # Compute sigma nought
        peak_power = get_waveforms_peak_power(l1.waveform.power)

        # Get Input parameter from l1 object
        tx_power = l1.get_parameter_by_name("classifier", "transmit_power")
        if tx_power is None:
            msg = "classifier `transmit_power` must exist for this pre-processor item -> aborting"
            logger.warning(msg)
            return
        altitude = l1.time_orbit.altitude

        # Compute absolute satellite velocity
        sat_vel_x = l1.get_parameter_by_name("classifier", "satellite_velocity_x")
        sat_vel_y = l1.get_parameter_by_name("classifier", "satellite_velocity_y")
        sat_vel_z = l1.get_parameter_by_name("classifier", "satellite_velocity_z")
        if sat_vel_x is None or sat_vel_y is None or sat_vel_z is None:
            msg = "classifier `satellite_velocity_[x|y|z]` must exist for this pre-processor item -> aborting"
            logger.warning(msg)
            return
        velocity = np.sqrt(sat_vel_x**2. + sat_vel_y**2. + sat_vel_z**2.)

        # Compute sigma_0
        sigma0 = get_sar_sigma0(peak_power, tx_power, altitude, velocity)
        sigma0[np.isinf(sigma0)] = np.nan

        # Add the classifier
        l1.classifier.add(peak_power, "peak_power")
        l1.classifier.add(sigma0, "sigma0")

    @property
    def required_options(self):
        return ["tfmra_leading_edge_start", "tfmra_leading_edge_center", "tfmra_leading_edge_end"]


class L1PWaveformPeakiness(DefaultLoggingClass):
    """
    A L1P pre-processor item class for computing leading edge width (full, first half, second half)
    using three TFMRA thresholds """

    def __init__(self, **cfg):
        super(L1PWaveformPeakiness, self).__init__(self.__class__.__name__)
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


class L1PLeadingEdgeQuality(DefaultLoggingClass):
    """
    Class to compute a leading edge width quality indicator
    Requires `first_maximum_index` classifier parameter
    """

    def __init__(self, **cfg):
        super(L1PLeadingEdgeQuality, self).__init__(self.__class__.__name__)
        self.cfg = cfg
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

            # import matplotlib.pyplot as plt
            # if leq[i] < 1:
            #     x = np.arange(wfm.shape[0])
            #     plt.figure(dpi=150)
            #     plt.plot(x, wfm)
            #     plt.plot(x[i0:i1], wfm[i0:i1], color="red")
            #     plt.scatter(x[fmi[i]], fmp[i])
            #     plt.annotate(f"radar_mode = {l1.radar_modes}", (5, 0.9))
            #     plt.annotate(f"leq[{i}] = {leq[i]:.3f}", (5, 0.8))
            #     plt.annotate(f"fmi power = {fmp[i]:.3f}", (5, 0.7))
            #     plt.annotate(f"power raise = {total_power_raise:.3f}", (5, 0.6))
            #     plt.show()

        # Add the classifier to the l1 object
        l1.classifier.add(leq, "leading_edge_quality")
        l1.classifier.add(fmi, "first_maximum_index")
        l1.classifier.add(fmp, "first_maximum_power")

    @property
    def required_options(self):
        return ["leading_edge_lookup_window", "first_maximum_normalized_power_threshold",
                "minimum_valid_first_maximum_index"]


class L1PLeadingEdgePeakiness(DefaultLoggingClass):

    def __init__(self, **cfg):
        super(L1PLeadingEdgePeakiness, self).__init__(self.__class__.__name__)
        self.cfg = cfg
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
