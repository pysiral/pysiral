# -*- coding: utf-8 -*-

from pysiral.errorhandler import FileIOErrorHandler
from pysiral.path import folder_from_filename

import numpy as np
import os


class Sentinel3SRALL1b(object):

    def __init__(self, raise_on_error=False):

        # Error Handling
        self._init_error_handling(raise_on_error)
        self.product_info = Sentinel3SRALProductInfo()
        self._radar_mode = "sar"
        self._filename = None
        self.n_records = 0
        self.range_bin_width = 0.234212857813 * 2.
        self.nominal_tracking_bin = 64
        self._xml_header_file = "xfdumanifest.xml"
        self._xml_metadata_object_index = {
            "processing": 0,
            "acquisitionPeriod": 1,
            "platform": 2,
            "generalProductInformation": 3,
            "measurementOrbitReference": 4,
            "measurementQualityInformation": 5,
            "measurementFrameSet": 6,
            "sralProductInformation": 7}

    def parse(self):
        self._parse_xml_header()
        self._parse_measurement_nc()

    def _parse_xml_header(self):
        """
        Parse the Sentinel-3 XML header file and extract key attributes
        for filtering
        """

        filename_header = os.path.join(
            folder_from_filename(self.filename), self._xml_header_file)
        self._xmlh = parse_sentinel3_l1b_xml_header(filename_header)

        # Extract Metadata
        metadata = self._xmlh["metadataSection"]["metadataObject"]

        # Extract Product Info
        index = self._xml_metadata_object_index["sralProductInformation"]
        sral_product_info = metadata[index]["metadataWrap"]["xmlData"]
        sral_product_info = sral_product_info["sralProductInformation"]

        # Save in
        sar_mode_percentage = sral_product_info["sral:sarModePercentage"]
        self.product_info.sar_mode_percentage = float(sar_mode_percentage)

        open_ocean_percentage = sral_product_info["sral:openOceanPercentage"]
        self.product_info.open_ocean_percentage = float(open_ocean_percentage)

        # print "sar_mode_percentage = %.1f" % int(sar_mode_percentage)

    def _parse_measurement_nc(self):

        from pysiral.iotools import ReadNC
        self._validate()
        self.nc = ReadNC(self.filename)
#        for attribute in self.nc.attributes:
#            print "attribute: %s = %s" % (
#                attribute, str(getattr(self.nc, attribute)))
#        stop
        for parameter in self.nc.parameters:
            print parameter, getattr(self.nc, parameter).shape
        # stop

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

        self.wfm_power = self.nc.waveform_20_ku
        n_records, n_range_bins = np.shape(self.wfm_power)
        self.n_records = n_records
        # Get the window delay
        # "The tracker_range_20hz is the range measured by the onboard tracker
        #  as the window delay, corrected for instrumental effects and
        #  CoG offset"
        tracker_range_20hz = self.nc.tracker_range_20_ku

        self.wfm_range = np.ndarray(shape=self.wfm_power.shape,
                                    dtype=np.float32)
        range_bin_index = np.arange(n_range_bins)
        for record in np.arange(n_records):
            self.wfm_range[record, :] = tracker_range_20hz[record] + \
                (range_bin_index*self.range_bin_width) - \
                (self.nominal_tracking_bin*self.range_bin_width)

    def _validate(self):
        pass


class Sentinel3SRALProductInfo(object):

    def __init__(self):

        self.sar_mode_percentage = None
        self.open_ocean_percentage = None
        self.start_time = None
        self.stop_time = None
        self.lat_min = None
        self.lat_max = None
        self.lon_min = None
        self.lon_max = None


def parse_sentinel3_l1b_xml_header(filename):
    """
    Reads the XML header file of a Sentinel 3 L1b Data set
    and returns the contents as an OrderedDict
    """
    import xmltodict
    with open(filename) as fd:
        content_odereddict = xmltodict.parse(fd.read())
    return content_odereddict[u'xfdu:XFDU']
