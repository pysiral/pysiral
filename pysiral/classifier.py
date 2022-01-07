# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 18:09:30 2015

@author: Stefan

Modified by FMI in August 2018:
LTPP added

:TODO: Merge into a waveform module

"""

import numpy as np
import bottleneck as bn


class BaseClassifier(object):

    def __init__(self):
        pass


class CS2OCOGParameter(BaseClassifier):
    """
    Calculate OCOG Parameters (Amplitude, Width) for CryoSat-2 waveform
    counts.
    Algorithm Source: retrack_ocog.pro from CS2AWI lib
    """

    def __init__(self, wfm_counts):
        super(CS2OCOGParameter, self).__init__()
        self._n = np.shape(wfm_counts)[0]
        self._amplitude = np.ndarray(shape=[self._n], dtype=np.float32)
        self._width = np.ndarray(shape=[self._n], dtype=np.float32)
        self._calc_parameters(wfm_counts)

    def _calc_parameters(self, wfm_counts):
        for i in np.arange(self._n):
            y = wfm_counts[i, :].flatten().astype(np.float64)
            y -= bn.nanmean(y[0:11])       # Remove Noise
            y[np.where(y < 0.0)[0]] = 0.0  # Set negative counts to zero
            y2 = y**2.0
            self._amplitude[i] = np.sqrt((y2**2.0).sum() / y2.sum())
            self._width[i] = ((y2.sum())**2.0) / (y2**2.0).sum()

    @property
    def amplitude(self):
        return self._amplitude

    @property
    def width(self):
        return self._width


class CS2PulsePeakiness(BaseClassifier):
    """
    Calculates Pulse Peakiness (full, left & right) for CryoSat-2 waveform
    counts
    XXX: This is a 1 to 1 legacy implementation of the IDL CS2AWI method,
         consistent method of L1bData or L2Data is required
    """
    def __init__(self, wfm_counts, pad=2):
        super(CS2PulsePeakiness, self).__init__()
        shape = np.shape(wfm_counts)
        self._n = shape[0]
        self._n_range_bins = shape[1]
        self._pad = pad
        dtype = np.float32
        self._peakiness = np.ndarray(shape=[self._n], dtype=dtype)*np.nan
        self._peakiness_r = np.ndarray(shape=[self._n], dtype=dtype)*np.nan
        self._peakiness_l = np.ndarray(shape=[self._n], dtype=dtype)*np.nan
        self._calc_parameters(wfm_counts)

    def _calc_parameters(self, wfm_counts):
        for i in np.arange(self._n):
            try:
                y = wfm_counts[i, :].flatten().astype(np.float32)
                y -= bn.nanmean(y[0:11])  # Remove Noise
                y[np.where(y < 0.0)[0]] = 0.0  # Set negative counts to zero
                yp = np.nanmax(y)  # Waveform peak value
                ypi = np.nanargmax(y)  # Waveform peak index
                if 3*self._pad < ypi < self._n_range_bins-4*self._pad:
                    self._peakiness_l[i] = yp/bn.nanmean(y[ypi-3*self._pad:ypi-1*self._pad+1])*3.0
                    self._peakiness_r[i] = yp/bn.nanmean(y[ypi+1*self._pad:ypi+3*self._pad+1])*3.0
                    self._peakiness[i] = yp/y.sum()*self._n_range_bins
            except ValueError:
                self._peakiness_l[i] = np.nan
                self._peakiness_r[i] = np.nan
                self._peakiness[i] = np.nan

    @property
    def peakiness(self):
        return self._peakiness

    @property
    def peakiness_r(self):
        return self._peakiness_r

    @property
    def peakiness_l(self):
        return self._peakiness_l


# Late tail to peak power (LTPP) ratio added for fmi needs
class S3LTPP(BaseClassifier):
    """
    Calculates Late-Tail-to-Peak-Power ratio.
    source: Rinne 2016
    """
    def __init__(self, wfm_counts, pad=1):
        # Warning: if 0padding is introduced in S3 L1 processing baseline, pad must be set to 2
        super(S3LTPP, self).__init__()
        shape = np.shape(wfm_counts)
        self._n = shape[0]
        self._n_range_bins = shape[1]
        self._pad = pad
        dtype = np.float32
        self._ltpp = np.ndarray(shape=[self._n], dtype=dtype)*np.nan
        self._calc_parameters(wfm_counts)

    def _calc_parameters(self, wfm_counts):
        # loop over the waveforms
        for i in np.arange(self._n): 
            try:
                y = wfm_counts[i, :].flatten().astype(np.float32)
                y -= bn.nanmean(y[0:11])  # Remove Noise
                y[np.where(y < 0.0)[0]] = 0.0  # Set negative counts to zero
                yp = np.nanmax(y)  # Waveform peak value
                
                if np.isnan(yp):  # if the current wf is nan
                    # no ltpp can be computed
                    self._ltpp[i] = np.nan
                else:
                    ypi = np.nanargmax(y)  # Waveform peak index
                    
                    # gates to compute the late tail:
                    # [ypi+50:ypi+70] if 0padding=2, [ypi+25:ypi+35] if 0padding=1
                    gate_start = ypi + self._pad*25
                    gate_stop = ypi + self._pad*35 + 1
                    
                    if gate_start > self._n_range_bins or gate_stop > self._n_range_bins:
                        # not enough gates to compute the LTPP
                        self._ltpp[i] = np.nan
                    else:
                        self._ltpp[i] = np.mean(y[gate_start:gate_stop])/yp
                        
            except ValueError:
                self._ltpp[i] = np.nan
                
    @property
    def ltpp(self):
        return self._ltpp        
            

class CS2LTPP(BaseClassifier):
    """
    Calculates Late-Tail-to-Peak-Power ratio.
    """
    def __init__(self, wfm_counts, pad=2):
        super(CS2LTPP, self).__init__()
        shape = np.shape(wfm_counts)
        self._n = shape[0]
        self._n_range_bins = shape[1]
        self._pad = pad
        dtype = np.float32
        self._ltpp = np.ndarray(shape=[self._n], dtype=dtype)*np.nan
        self._calc_parameters(wfm_counts)

    def _calc_parameters(self, wfm_counts):
        # loop over the waveforms
        for i in np.arange(self._n): 
            y = wfm_counts[i, :].flatten().astype(np.float32)
            y -= bn.nanmean(y[0:11])  # Remove Noise
            y[np.where(y < 0.0)[0]] = 0.0  # Set negative counts to zero
            yp = np.nanmax(y)  # Waveform peak value
            ypi = bn.nanargmax(y)  # Waveform peak index

            # AMANDINE: implementation for wf of 256 bins (128 zero-padded)?
            onediv = float(1)/float(41)

            # AMANDINE: here i seems to be understood as gate index but it is wf index!?
            if i == 256: 
                break
            if [i > (ypi + 100)] and [i < (ypi + 140)]:             # AMANDINE: syntax to be checked
                try:
                    self._ltpp[i] = (onediv*float(y[i]))/float(yp)  # AMANDINE: where is the sum in this formula?
                except ZeroDivisionError:
                    self._ltpp[i] = np.nan

    @property
    def ltpp(self):
        return self._ltpp


class EnvisatWaveformParameter(BaseClassifier):
    """
    Currently only computes pulse peakiness for Envisat waveforms
    from SICCI processor.

    Parameter for Envisat from SICCI Processor
        skip = 5
        bins_after_nominal_tracking_bin = 83
    """

    def __init__(self, wfm, skip=5, bins_after_nominal_tracking_bin=83):
        super(EnvisatWaveformParameter, self).__init__()
        self.t_n = bins_after_nominal_tracking_bin
        self.skip = skip
        self._n = wfm.shape[0]
        self._n_range_bins = wfm.shape[1]
        self._init_parameter()
        self._calc_parameter(wfm)

    def _init_parameter(self):
        self.peakiness_old = np.ndarray(shape=self._n, dtype=np.float32)
        self.peakiness = np.ndarray(shape=self._n, dtype=np.float32)*np.nan

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
                self.peakiness[i] = float(max(wave))/float(sum(wave))*self._n_range_bins
            except ZeroDivisionError:
                self.peakiness[i] = np.nan
