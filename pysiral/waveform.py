# -*- coding: utf-8 -*-
"""
Created on Fri Jul 01 13:07:10 2016

@author: shendric
"""

import numpy as np


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


class TFMRALeadingEdgeWidth(object):
    """
    Container for computation of leading edge width by taking differences
    between first maximum power thresholds
    """

    def __init__(self, rng, wfm, radar_mode, is_ocean):
        from retracker import TFMRA
        # Compute filtered waveform and index of first maximum once
        self.tfmra = TFMRA()
        self.tfmra.set_default_options()
        filt_rng, filt_wfm, fmi, norm = self.tfmra.get_preprocessed_wfm(
            rng, wfm, radar_mode, is_ocean)
        self.wfm, self.rng, self.fmi = filt_wfm, filt_rng, fmi

    def get_width_from_thresholds(self, thres0, thres1):
        """ returns the width between two thresholds in the range [0:1] """
        width = self.tfmra.get_thresholds_distance(
            self.rng, self.wfm, self.fmi, thres0, thres1)
        return width
