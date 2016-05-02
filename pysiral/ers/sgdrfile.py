# -*- coding: utf-8 -*-

from pysiral.errorhandler import FileIOErrorHandler

import numpy as np
import os


class ERSSGDR(object):

    def __init__(self, raise_on_error=False):

        # Error Handling
        self._init_error_handling(raise_on_error)
        self._baseline = None
        self._radar_mode = "lrm"
        self._filename = None
        self.n_records = 0
        self.n_blocks = 20
        self.range_bin_width = 0.4545
        self.nominal_tracking_bin = 32.5

    def parse(self):
        from pysiral.iotools import ReadNC
        self._validate()
        self.nc = ReadNC(self.filename)
#        for attribute in self.nc.attributes:
#            print "attribute: %s = %s" % (
#                attribute, str(getattr(self.nc, attribute)))
#        for parameter in self.nc.parameters:
#            print parameter, getattr(self.nc, parameter).shape

    def get_status(self):
        # XXX: Not much functionality here
        return False

    def post_processing(self):
        """
        The SGDR data structure needs to be streamlined, so that it
        is easy to grab the relevant parameters as indiviual arrays
        """
        self._prepare_waveform_power_and_range()

    @property
    def filename(self):
        return self._filename

    @filename.setter
    def filename(self, filename):
        """ Save and validate filenames for header and product file """
        # Test if valid file first
        self._error.file_undefined = not os.path.isfile(filename)
        if self._error.file_undefined:
            return
        self._filename = filename

    @property
    def radar_mode(self):
        return self._radar_mode

    def _init_error_handling(self, raise_on_error):
        self._error = FileIOErrorHandler()
        self._error.raise_on_error = raise_on_error
        self._error.file_undefined = True

    def _prepare_waveform_power_and_range(self):
        """
        reforms the waveform to computes the corresponding range for each
        range bin
        """
        records_per_block, n_blocks, n_range_bins = self.nc.ku_wf.shape
        n_records = records_per_block*n_blocks
        self.n_records = n_records
        target_shape = (n_records, n_range_bins)
        self.wfm_power = np.reshape(self.nc.ku_wf, target_shape).astype(
            np.uint16)

        # Get the window delay
        # "The tracker_range_20hz is the range measured by the onboard tracker
        #  as the window delay, corrected for instrumental effects and
        #  CoG offset"
        tracker_range_20hz = self.nc.tracker_range_20hz.flatten()
        range_bin_index = np.arange(n_range_bins)
        self.wfm_range = np.ndarray(shape=target_shape, dtype=np.float32)
        for record in np.arange(n_records):
            self.wfm_range[record, :] = tracker_range_20hz[record] + \
                (range_bin_index*self.range_bin_width) - \
                (self.nominal_tracking_bin*self.range_bin_width)

    def _validate(self):
        pass
