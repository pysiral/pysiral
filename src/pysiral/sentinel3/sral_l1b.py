# -*- coding: utf-8 -*-


from pathlib import Path

import numpy as np


class Sentinel3SRALL1b(object):

    def __init__(self):

        # Error Handling
        self.product_info = Sentinel3SRALProductInfo()
        self._filename = None
        self.n_records = 0

    def parse(self):
        self._parse_measurement_nc()

    def parse_xml_header(self, settings):
        """
        Parse the Sentinel-3 XML header file and extract key attributes
        for filtering
        """

        # Retrieve header information from mission settings
        xml_header_file = settings.xml_header_file
        xml_metadata_object_index = settings.xml_metadata_object_index

        filename_header = Path(self.filename).parent / xml_header_file
        self._xmlh = parse_sentinel3_l1b_xml_header(filename_header)

        # Extract Metadata
        metadata = self._xmlh["metadataSection"]["metadataObject"]

        # Extract General Product Info
        index = xml_metadata_object_index["generalProductInformation"]
        product_info = metadata[index]["metadataWrap"]["xmlData"]
        product_info = product_info["sentinel3:generalProductInformation"]

        timeliness_dict = dict(NR="NRT", ST="STC", NT="NTC")
        timeliness = timeliness_dict[product_info["sentinel3:timeliness"]]
        self.product_info.timeliness = timeliness

        # Extract SRAL Product Info
        index = xml_metadata_object_index["sralProductInformation"]
        sral_product_info = metadata[index]["metadataWrap"]["xmlData"]
        sral_product_info = sral_product_info["sralProductInformation"]

        # Save in
        sar_mode_percentage = sral_product_info["sral:sarModePercentage"]
        self.product_info.sar_mode_percentage = float(sar_mode_percentage)

        open_ocean_percentage = sral_product_info["sral:openOceanPercentage"]
        self.product_info.open_ocean_percentage = float(open_ocean_percentage)

        # Get regional coverage (important for NRT data that comes in 10
        # minute granules)
        index = xml_metadata_object_index["measurementFrameSet"]
        measurement_frameset = metadata[index]["metadataWrap"]["xmlData"]
        positions = measurement_frameset["sentinel-safe:frameSet"]
        positions = positions["sentinel-safe:footPrint"]["gml:posList"]
        lons = [float(p) for p in positions.split()[1::2]]
        lats = [float(p) for p in positions.split()[0::2]]
        self.product_info.lat_max = np.nanmax(lats)
        self.product_info.lat_min = np.nanmin(lats)
        self.product_info.lon_max = np.nanmax(lons)
        self.product_info.lon_min = np.nanmin(lons)

    def _parse_measurement_nc(self):

        from pysiral.core.iotools import ReadNC
        self._validate()

        # Read the L2 netCDF file
        self.nc = ReadNC(self.filename, nan_fill_value=True)

    def get_status(self):
        # XXX: Not much functionality here
        return False

    def post_processing(self):
        pass

    @property
    def filename(self):
        return self._filename

    @filename.setter
    def filename(self, filename):
        """ Save and validate filenames for header and product file """
        # Test if valid file first
        if not Path(filename).is_file():
            return
        self._filename = filename

    @property
    def radar_mode(self):
        return self._radar_mode

    def _prepare_waveform_power_and_range(self):
        """
        reforms the waveform to computes the corresponding range for each
        range bin
        """

        # self.wfm_power = self.nc.waveform_20_ku
        self.wfm_power = self.nc.i2q2_meas_ku_l1b_echo_sar_ku

        n_records, n_range_bins = shape = self.wfm_power.shape
        self.n_records = n_records

        # Get the window delay
        # "The tracker_range_20hz is the range measured by the onboard tracker
        #  as the window delay, corrected for instrumental effects and
        #  CoG offset"
        tracker_range_20hz = self.nc.range_ku_l1b_echo_sar_ku

        self.wfm_range = np.ndarray(shape=shape, dtype=np.float32)
        range_bin_index = np.arange(n_range_bins)
        for record in np.arange(n_records):
            self.wfm_range[record, :] = tracker_range_20hz[record] + \
                (range_bin_index*self.range_bin_width) - \
                (self.nominal_tracking_bin * self.range_bin_width)

    def _validate(self):
        pass


class Sentinel3SRALProductInfo(object):

    def __init__(self):

        self.sar_mode_percentage = None
        self.open_ocean_percentage = None
        self.timeliness = None
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
    with open(str(filename)) as fd:
        content_odereddict = xmltodict.parse(fd.read())
    return content_odereddict[u'xfdu:XFDU']
