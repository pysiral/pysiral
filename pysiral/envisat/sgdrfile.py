# -*- coding: utf-8 -*-

from pysiral.errorhandler import FileIOErrorHandler
from pysiral.envisat.sgdr_mds_def import envisat_get_mds_def
from pysiral.esa.header import (ESAProductHeader, ESAScienceDataSetDescriptors)

import numpy as np
import os
import re


class Envisat18HzArrays(object):
    """ Reforms the grouped SGDR data into continuous arrays """

    def __init__(self):
        self.registered_parameters = []
        self.n_records = 0
        self.n_blocks = 20

    def reform_timestamp(self, mds):
        """ Creates an array of datetime objects for each 18Hz record """
        # XXX: Current no microsecond correction
        timestamp = np.ndarray(shape=(self.n_records), dtype=object)
        mdsr_timestamp = get_structarr_attr(mds, "utc_timestamp")
        for i in range(self.n_records):
            timestamp[i] = mdsr_timestamp_to_datetime(mdsr_timestamp[i])
        self.timestamp = np.repeat(timestamp, self.n_blocks)

    def reform_position(self, mds):
        # 0) Get the relevant data blocks
        time_orbit = get_structarr_attr(mds, "time_orbit")
        range_block = get_structarr_attr(mds, "range_information")
        # 1) get the 1Hz data from the time_orbit group
        longitude = get_structarr_attr(time_orbit, "longitude")
        latitude = get_structarr_attr(time_orbit, "latitude")
        altitude = get_structarr_attr(time_orbit, "altitude")
        # 2) Get the increments
        lon_inc = get_structarr_attr(range_block, "18hz_longitude_differences")
        lat_inc = get_structarr_attr(range_block, "18hz_latitude_differences")
        alt_inc = get_structarr_attr(time_orbit, "18hz_altitude_differences")
        # 3) Expand the 1Hz position arrays
        self.longitude = np.repeat(longitude, self.n_blocks)
        self.latitude = np.repeat(latitude, self.n_blocks)
        self.altitude = np.repeat(altitude, self.n_blocks)
        # 4) Apply the increments
        # XXX: Current version of get_struct_arr returns datatype objects
        #      for listcontainers => manually set dtype
        self._apply_18Hz_increment(self.longitude, lon_inc.astype(np.float32))
        self._apply_18Hz_increment(self.latitude, lat_inc.astype(np.float32))
        self._apply_18Hz_increment(self.altitude, alt_inc.astype(np.float32))

    def reform_waveform(self, mds):
        n = self.n_records * self.n_blocks
        self.power = np.ndarray(shape=(n, 128), dtype=np.float32)
        for dsd in range(self.n_records):
            for block in range(self.n_blocks):
                i = dsd*self.n_blocks + block
                self.power[i, :] = mds[dsd].wfm[block].average_wfm_if_corr_ku

    def _apply_18Hz_increment(self, data, inc):
        for i in range(self.n_records):
            i0, i1 = i*self.n_blocks, (i+1)*self.n_blocks
            data[i0:i1] += inc[i, :]


class EnvisatSGDR(object):

    _DS_NAME = {
        "ra2": "RA2_DATA_SET_FOR_LEVEL_2",
        "wfm18hz": "RA2_AVERAGE_WAVEFORMS"}

    def __init__(self, raise_on_error=False):

        # Error Handling
        self._init_error_handling(raise_on_error)
        self._baseline = None
        self._radar_mode = "lrm"
        self._filename = None
        self.mph = None
        self.sph = None
        self.dsd = None
        self.mds_ra2 = None
        self.mds_wfm18hz = None
        self.mds_18hz = Envisat18HzArrays()
        # XXX: nor sure if needed
        # self.mds_mwr = None
        # self.mds_wfmburst = None

    def parse(self):
        self._validate()
        with open(self._filename, "r") as self._fh:
            self._parse_mph()
            self._parse_sph()
            self._parse_dsd()
        with open(self._filename, "rb") as self._fh:
            self._parse_mds("ra2")
            self._parse_mds("wfm18hz")
            # XXX: mwr: not necessary?, wfmburst: not in the files?
            # self._parse_mds("mwr")
            # self._parse_mds("wfmburst")

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
    def baseline(self):
        return self._baseline

    @property
    def radar_mode(self):
        return self._radar_mode

    def _init_error_handling(self, raise_on_error):
        self._error = FileIOErrorHandler()
        self._error.raise_on_error = raise_on_error
        self._error.file_undefined = True

    def _parse_mph(self):
        self.mph = EnvisatMainProductHeader()
        self._read_header_lines(self.mph)

    def _parse_sph(self):
        self.sph = EnvisatSpecificProductHeader()
        self._read_header_lines(self.sph)

    def _parse_dsd(self):
        """ Reads the Data Set Descriptors dsd's in the SGDR header """
        self.dsd = ESAScienceDataSetDescriptors()
        self.n_dsd_lines = self.dsd.get_num_lines(self.mph.num_dsd)
        for i in np.arange(self.n_dsd_lines+1):
            line = self._fh.readline()
            self.dsd.parse_line(line)

    def _parse_mds(self, mds_target):
        """ Read the data blocks """
        # Just reopened the file in binary mode -
        # > get start byte and number of data set records
        l1b_data_set_name = self._DS_NAME[mds_target].lower()
        data_set_descriptor = self.dsd.get_by_fieldname(l1b_data_set_name)
        startbyte = int(data_set_descriptor["ds_offset"])
        self.n_msd_records = int(data_set_descriptor["num_dsr"])
        # Set the file pointer
        self._fh.seek(startbyte)
        # Get the parser (depending on MDS target)
        self.mds_definition = envisat_get_mds_def(
            self.n_msd_records, mds_target)
        mds_parser = self.mds_definition.get_mds_parser()
        # Parser the binary part of the .DBL file
        mds = mds_parser.parse(self._fh.read(mds_parser.sizeof()))
        setattr(self, "mds_"+mds_target, mds)

    def _read_header_lines(self, header):
        """ Method to read the MPH and SPH headers are identical """
        line_index = 0
        while line_index < 1000:
            line = self._fh.readline()
            header.scan_line(line)
            # Break if the first SPH field is reached
            if header.last_field(line):
                break
            line_index = +1

    def _validate(self):
        pass


class EnvisatMainProductHeader(ESAProductHeader):
    """
    Container for the Main Product Header of Envisat SGDR data
    """

    # Datatypes for automatic conversion
    # A field that is in here is converted to an 32bit integer
    _INT_LIST = [
        "PHASE", "CYCLE",  "REL_ORBIT", "ABS_ORBIT", "SAT_BINARY_TIME",
        "CLOCK_STEP", "LEAP_SIGN", "LEAP_ERR", "PRODUCT_ERR", "TOT_SIZE",
        "SPH_SIZE", "NUM_DSD", "DSD_SIZE", "NUM_DATA_SETS"]

    # A field that is in here is converted to an 32bit float
    _FLOAT_LIST = [
        "DELTA_UT1", "X_POSITION", "Y_POSITION",
        "Z_POSITION", "X_VELOCITY", "Y_VELOCITY", "Z_VELOCITY"]

    def __init__(self):
        super(EnvisatMainProductHeader, self).__init__()

    def last_field(self, line):
        return re.search("NUM_DATA_SETS=", line)


class EnvisatSpecificProductHeader(ESAProductHeader):
    """
    Container for the Specific Product Header of Envisat SGDR science data
    """

    # Datatypes for automatic conversion
    # A field that is in here is converted to an 32bit integer
    _INT_LIST = [
        "RA2_L2_PROC_FLAG", "RA2_L1B_PROC_FLAG", "RA2_L1B_HEADER_FLAG",
        "RA2_FLAG_MANOEUVER", "RA2_IF_MASK_SEL", "RA2_IF_MASK_PROC",
        "RA2_USO_SEL", "RA2_USO_PROC", "AVERAGE_GLOBAL_PRESSURE",
        "SOLAR_ACTIVITY_INDEX", "MWR_L2_PROC_FLAG", "MWR_L1B_PROC_FLAG",
        "MWR_L1B_HEADER_FLAG", "MWR_L1B_TELEMETRY_FLAG"]

    # A field that is in here is converted to an 32bit float
    _FLOAT_LIST = [
        "RA2_FIRST_LAT", "RA2_FIRST_LONG", "RA2_LAST_LAT",
        "RA2_LAST_LONG", "RA2_L2_PROCESSING_QUALITY",
        "RA2_L1B_PROCESSING_QUALITY", "RA2_L1B_HEADER_QUALITY",
        "RA2_L2_PROC_THRESH", "RA2_L1B_PROC_THRESH",
        "RA2_L1B_HEADER_THRESH", "RA2_MEASUREMENT_PERCENT",
        "RA2_320_BAND_PERCENT", "RA2_80_BAND_PERCENT", "RA2_20_BAND_PERCENT",
        "RA2_OCEAN_KU_RETRACK_PERCENT", "RA2_OCEAN_S_RETRACK_PERCENT",
        "RA2_ICE1_KU_RETRACK_PERCENT", "RA2_ICE1_S_RETRACK_PERCENT",
        "RA2_ICE2_KU_RETRACK_PERCENT", "RA2_ICE2_S_RETRACK_PERCENT",
        "RA2_SEAICE_KU_RETRACK_PERCENT", "RA2_PEAKINESS_LOW_PERCENT",
        "RA2_PEAKINESS_HIGH_PERCENT", "MWR_BT_OPTIMAL_INTERPOLATION_PERCENT",
        "RA2_TIME_SHIFT_MIDFRAME", "RA2_TIME_INTERVAL", "MWR_FIRST_LAT",
        "MWR_FIRST_LONG", "MWR_LAST_LAT", "MWR_LAST_LONG",
        "MWR_L2_PROC_QUALITY", "MWR_L1B_PROC_QUALITY", "MWR_L1B_HEAD_QUALITY",
        "MWR_L1B_TELEM_QUALITY", "MWR_L2_PROC_THRESH", "MWR_L1B_PROC_THRESH",
        "MWR_L1B_HEAD_THRESH", "MWR_L1B_TELEM_THRESH",
        "RA2_WS_OPTIMAL_INTERPOLATION_PERCENT", "MWR_LANDFLAG_PERCENT",
        "MWR_SEAFLAG_PERCENT"]

    def __init__(self):
        super(EnvisatSpecificProductHeader, self).__init__()

    def last_field(self, line):
        return re.search("MWR_SEAFLAG_PERCENT=", line)
