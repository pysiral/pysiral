# -*- coding: utf-8 -*-

"""
This is the SAMOSA+ retracker variant that uses `samosa-waveform-model`, a
re-implementation of the SAMPy package. The objective of the re-implementation
is better performance (code optimization and multiprocessing) and a greater flexibility
for retracking settings (sub-waveform retracking, limiting parameters of the fit)
"""

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"

from pysiral.retracker import BaseRetracker


class SAMOSAWaveformModelFit(BaseRetracker):

    def __init__(self) -> None:
        super().__init__()
