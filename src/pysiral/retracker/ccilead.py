# -*- coding: utf-8 -*-

"""
The SICCILead retracker uses a waveform fitting method for pulse-limited Lead waveforms.
Its implementation comes from the first phase of the ESA Climate Change Initiative on Sea Ice
"""

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"

from typing import Any

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

from core.flags import ANDCondition
from retracker import BaseRetracker


class SICCILead(BaseRetracker):

    def __init__(self):
        super(SICCILead, self).__init__()
        self._retracker_properties = {}

    def create_retracker_properties(self, n_records):
        parameter = [
            "retracked_bin",
            "maximum_power_bin",
            "sigma",
            "k",
            "alpha",
            "power_in_echo_tail",
            "rms_echo_and_model"]
        for parameter_name in parameter:
            self._retracker_properties[parameter_name] = np.ndarray(shape=n_records, dtype=np.float32) * np.nan

    def l2_retrack(self, rng, wfm, indices, radar_mode, is_valid):
        # Run the retracker
        self._sicci_lead_retracker(rng, wfm, indices)
        # Filter the results
        if self._options.filter.use_filter:
            self._filter_results()

    def _sicci_lead_retracker(self, rng, wfm, indices):

        # retracker options (see l2 settings file)
        skip = self._options.skip_first_bins
        initial_guess = list(self._options.initial_guess)
        maxfev = self._options.maxfev
        time = np.arange(wfm.shape[1]-skip).astype(float)
        x = np.arange(wfm.shape[1])

        # Loop over lead indices
        for index in indices:
            wave = wfm[index, skip:]
            initial_guess[0] = np.argmax(wave)
            initial_guess[3] = np.max(wave)
            popt, cov, *_ = curve_fit(
                pl_lead_waveform_model,
                time,
                wave.astype(float),
                p0=initial_guess,
                maxfev=maxfev
            )

            # Store retracker parameter for filtering
            # tracking point in units of range bins
            # tracking point in units of range bins
            self.retracked_bin[index] = skip + popt[0]
            self.k[index] = popt[1]
            self.sigma[index] = popt[2]
            self.alpha[index] = popt[3]
            self.maximum_power_bin[index] = np.argmax(wave)

            # Get derived parameter
            self.power_in_echo_tail[index] = power_in_echo_tail(
                wfm[index, :], self.retracked_bin[index], self.alpha[index])
            self.rms_echo_and_model[index] = rms_echo_and_model(
                wfm[index, :], self.retracked_bin[index],
                self.k[index], self.sigma[index], self.alpha[index])

            # Get the range by interpolation of range bin location
            try:
                self._range[index] = interp1d(
                    x, rng[index, :], kind='linear', copy=False)(
                        self.retracked_bin[index])
            except ValueError:
                self._range[index] = np.nan
                self._range[index] = np.nan

    def _filter_results(self):
        """ Filter the lead results based on threshold defined in SICCI """

        thrs = self._options.filter
        clf = self._classifier

        valid = ANDCondition()
        valid.add(self.sigma < thrs.maximum_std_of_gaussion_rise)
        valid.add(self.maximum_power_bin > thrs.minimum_bin_count_maxpower)

        # sea ice backscatter not available for ERS?
        if thrs.minimum_echo_backscatter is not None:
            valid.add(clf.sigma0 > thrs.minimum_echo_backscatter)

        bin_seperation = np.abs(self.retracked_bin - self.maximum_power_bin)
        valid.add(bin_seperation < thrs.maximum_retracker_maxpower_binsep)

        valid.add(self.power_in_echo_tail < thrs.maximum_power_in_echo_tail)
        valid.add(self.rms_echo_and_model < thrs.maximum_rms_echo_model_diff)
        valid.add(self.retracked_bin > thrs.sensible_lead_retracked_bin[0])
        valid.add(self.retracked_bin < thrs.sensible_lead_retracked_bin[1])

        # Error flag is also computed for other surface types, do not
        # overide those
        error_flag = self._flag
        error_flag[self.indices] = np.logical_not(valid.flag[self.indices])
        self._flag = error_flag

    def __getattr__(self, item: str) -> Any:
        """
        Direct attribute access to the cfg dictionary

        :param item:
        :return:
        """
        if item in self._retracker_properties:
            return self._retracker_properties[item]
        else:
            raise AttributeError(f"attribute {item} not found in retracker properties")


def pl_lead_waveform_model(t, t_0, k, sigma, a):
    """ Lead waveform model (for SICCILead curve fitting) """
    # Time for F to be F_L
    t_b = k*sigma**2
    # Helper coefficient
    sq_ktb = np.sqrt(k*t_b)
    # Polynomial coefficients for F_L
    aa = ((5*k*sigma) - (4*sq_ktb)) / (2*sigma*t_b*sq_ktb)
    aaa = ((2*sq_ktb)-(3*k*sigma)) / (2*sigma*t_b*t_b*sq_ktb)
    # We are mostly interested in time in reference to t_0
    t_diff = (t - t_0)
    # F is piecewise. These are (Truth where applicable) * (Function)
    F_1 = (t <= t_0) * (t_diff/sigma)
    F_2 = (t >= (t_b+t_0)) * (np.sqrt(np.abs(t_diff)*k))
    F_L = np.logical_and(t > t_0, t < (t_b+t_0)) * (aaa*t_diff**3 + aa*t_diff**2+t_diff/sigma)
    # Compile F
    F = F_1 + F_L + F_2
    return a*np.exp(-F*F)  # Return e^-f^2(t)


def rms_echo_and_model(wfm, retracked_bin, k, sigma, alpha):
    """
    The root sum squared difference between the echo and the fitted function
    in the lead retracking is computed. The 5 bins before the tracking point
    are used as the echo rise.
    """
    tracking_point = int(retracked_bin)
    time = np.arange(len(wfm)).astype(float)
    modeled_wave = pl_lead_waveform_model(time, retracked_bin, k, sigma, alpha)
    diff = wfm[tracking_point-4:tracking_point+1] - modeled_wave[tracking_point-4:tracking_point+1]
    return np.sqrt(np.sum(diff*diff)/5)/alpha


def power_in_echo_tail(wfm, retracked_bin, alpha, pad=3):
    """
    The tail power is computed by summing the bin count in all bins from
    2 bins beyond the tracking point to the end of the range gate, then
    normalised by dividing by the value of alpha returned from the lead
    retracking
    source: SICCI
    """
    tracking_point = int(retracked_bin)
    return sum(wfm[tracking_point+pad:])/alpha
