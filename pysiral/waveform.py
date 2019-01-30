# -*- coding: utf-8 -*-
"""
Created on Fri Jul 01 13:07:10 2016

@author: shendric
"""

from retracker import SICCI2TfmraEnvisat
import numpy as np

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

    def __init__(self, rng, wfm, radar_mode, is_ocean):
        # Compute filtered waveform and index of first maximum once
        self.tfmra = SICCI2TfmraEnvisat()
        self.tfmra.set_default_options()
        filt_rng, filt_wfm, fmi, norm = self.tfmra.get_preprocessed_wfm(
            rng, wfm, radar_mode, is_ocean)
        self.wfm, self.rng, self.fmi = filt_wfm, filt_rng, fmi

    def get_width_from_thresholds(self, thres0, thres1):
        """ returns the width between two thresholds in the range [0:1] """
        width = self.tfmra.get_thresholds_distance(
            self.rng, self.wfm, self.fmi, thres0, thres1)
        return width


class L1PLeadingEdgeWidth(DefaultLoggingClass):
    """
    A L1P pre-processor item class for computing leading edge width (full, first half, second half)
    using three TFMRA thresholds """

    def __init__(self, **cfg):
        super(L1PLeadingEdgeWidth, self).__init__(self.__class__.__name__)
        for option_name in self.required_options:
            option_value = cfg.get(option_name, None)
            if option_value is None:
                msg = "Missing option `%s` -> No computation of leading edge width!" % option_name
                self.log.warning(msg)
            setattr(self, option_name, option_value)

    def apply(self, l1):
        """
        API class for the Level-1 pre-processor. Functionality is compute leading edge width (full, first half &
        second half) and adding the result to the classifier data group
        :param l1: A Level-1 data instance
        :return: None, Level-1 object is change in place
        """

        # Prepare input
        wfm = l1.waveform.power
        rng = l1.waveform.range
        radar_mode = l1.waveform.radar_mode
        is_ocean = l1.surface_type.get_by_name("ocean").flag
        thrs_start = self.tfmra_leading_edge_start
        thrs_center = self.tfmra_leading_edge_center
        thrs_end = self.tfmra_leading_edge_end

        # Compute the leading edge width (requires TFMRA retracking)
        width = TFMRALeadingEdgeWidth(rng, wfm, radar_mode, is_ocean)
        lew = width.get_width_from_thresholds(thrs_start, thrs_end)
        lew1 = width.get_width_from_thresholds(thrs_start, thrs_center)
        lew2 = width.get_width_from_thresholds(thrs_center, thrs_end)

        # Add result to classifier group
        l1.classifier.add(lew, "leading_edge_width")
        l1.classifier.add(lew1, "leading_edge_width_first_half")
        l1.classifier.add(lew2, "leading_edge_width_second_half")
        l1.classifier.add(width.fmi, "first_maximum_index")

    @property
    def required_options(self):
        return ["tfmra_leading_edge_start", "tfmra_leading_edge_center", "tfmra_leading_edge_end"]


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
            self.log.warning(msg)
            return
        altitude = l1.time_orbit.altitude

        # Compute absolute satellite velocity
        sat_vel_x = l1.get_parameter_by_name("classifier", "satellite_velocity_x")
        sat_vel_y = l1.get_parameter_by_name("classifier", "satellite_velocity_y")
        sat_vel_z = l1.get_parameter_by_name("classifier", "satellite_velocity_z")
        if sat_vel_x is None or sat_vel_y is None or sat_vel_z is None:
            msg = "classifier `satellite_velocity_[x|y|z]` must exist for this pre-processor item -> aborting"
            self.log.warning(msg)
            return
        velocity = np.sqrt(sat_vel_x**2. + sat_vel_y**2. + sat_vel_z**2.)

        # Compute sigma_0
        sigma0 = get_sar_sigma0(peak_power, tx_power, altitude, velocity)

        # Add the classifier
        l1.classifier.add(peak_power, "peak_power")
        l1.classifier.add(sigma0, "sigma0")