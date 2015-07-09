# -*- coding: utf-8 -*-
"""
Created on Tue Jul 07 14:10:34 2015

@author: Stefan
"""

import os
import re
import numpy as np
from treedict import TreeDict
import parse
from dateutil import parser as dtparser
import xmltodict


class L1bData(object):
    """
    Parent Class for all L1b data sets
    """
    def __init__(self):

        self._waveform_collection = None
        self._product_header = None


class CryoSatL1B(L1bData):
    """
    Handles CryoSat L1b Data (Header and Product file)
    """

    _VALID_BASELINES = ["B001", "C001"]

    def __init__(self):

        super(CryoSatL1B, self).__init__()

        self._baseline = None
        self._filename_header = None
        self._filename_product = None
        self._mph = None
        self._sph = None
        self._xmlh = None

    @property
    def filename(self):
        return self._filename_product

    @filename.setter
    def filename(self, filename):
        self._filename_product = filename
        self._filename_header = os.path.splitext(filename)[0]+".HDR"

    def parse(self):
        """ Parse the content of the L1B file """
        # First detect the baseline from the file name
        self._detect_baseline()
        # Parse the xml header file
        self._parse_header_file()
        # Parse the product file
        self._parse_product_file()

    def _detect_baseline(self):
        info = parse_cryosat_l1b_filename(self._filename_product)
        self.baseline = info.baseline

    def _parse_header_file(self):
        self._xmlh = parse_cryosat_l1b_xml_header(self._filename_header)

    def _parse_product_file(self):
        """
        Reads the content of *L1B*.DBL files
            - main product header (mph)
            - specific product header (sph)
            - data set descriptors (dsd)
            - data_block
        """
        with open(self._filename_product, "r") as self._fh:
            self._parse_mph()
            self._parse_sph()
            self._parse_dsd()
            self._parse_data_blocks()

    def _parse_mph(self):
        """
        Reads the main product header (mph) of a CryoSat-2 L1b file
        """
        # Save the mph information in its own class
        self.mph = CS2L1bMainProductHeader()
        self._read_header_lines(self.mph)

    def _parse_sph(self):
        """
        Reads the main product header (mph) of a CryoSat-2 L1b file
        """
        # Save the mph information in its own class
        self.sph = CS2L1bSpecificProductHeader()
        self._read_header_lines(self.sph)

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

    def _parse_dsd(self):
        pass

    def _parse_data_blocks(self):
        pass


class CS2L1bBaseHeader(object):
    """
    Parent class for both MPH and SPH header informations
    """

    def __init__(self):
        self.unit_dict = {}
        self.dtype_dict = {}
        self.field_list = []

    def _add_field(self, tag, value):
        setattr(self, tag, value)

    def scan_line(self, line):
        """ Parses a line from the MPH header of the data product file """
        parser = parse.compile("{tag}={value}")
        match = parser.parse(line)
        if match:
            tag = match["tag"].lower()
            # clean string expression and get unit of field value
            value, unit = self._parse_field_line(match["value"])
            # convert the unit to native data type
            value, dtype = self._value_convert_dtype(tag, value)
            # Save to data structure
            setattr(self, tag, value)
            self.field_list.append(tag)
            self.dtype_dict[tag] = dtype
            self.unit_dict[tag] = unit

    def _parse_field_line(self, value):
        """ Processes the raw value string of an MPH field """
        # clear double quotes and extra whitespaces from strings
        value = value.replace("\"", "").strip()
        # check if unit is given
        unit_parser = parse.compile("{value}<{unit}>")
        unit_match = unit_parser.parse(value)
        unit = None
        if unit_match:
            value = unit_match["value"]
            unit = unit_match["unit"]
        # Automatical convert in the unit to the proper data type
        return value, unit

    def _value_convert_dtype(self, tag, value):
        """ Convert the value to either integer or float """
        # TODO: The explicit list for field is not a good solution
        #       Proper way needed to do this
        dtype = "str"
        if tag.upper() in self._INT_LIST:
            dtype = "int32"
            value = np.int32(value)
        if tag.upper() in self._FLOAT_LIST:
            dtype = "float32"
            value = np.float32(value)
        return value, dtype

    def get_by_fieldname(self, fieldname):
        """
        Convinience function that returns the value, data type
        of an MPH records.

        Arguments:
            fieldname: (str):
                Name of the MPH field, case insensitive

        Output:
            value:
                value of the field in its respective data type
                None: if invalid field name

            unit: (str)
                Unit label
                None: if invalid field name
                None: if specific unit not given in the mph
        """
        if fieldname not in self.field_list:
            return None, None
        return getattr(self, fieldname), self.unit_dict[fieldname]

    def __repr__(self):
        output = ""
        for field in self.field_list:
            output += "%s=%s" % (field, getattr(self, field))
            output += ", unit="+str(self.unit_dict[field])
            output += ", dtype="+str(self.dtype_dict[field])
            output += "\n"
        return output


class CS2L1bMainProductHeader(CS2L1bBaseHeader):
    """
    Container for the Main Product Header of CryoSat-2 L1B science data
    """

    # Datatypes for automatic conversion
    # A field that is in here is converted to an 32bit integer
    _INT_LIST = [
        "PHASE", "CYCLE",  "REL_ORBIT", "ABS_ORBIT", "SAT_BINARY_TIME",
        "CLOCK_STEP", "LEAP_SIGN", "LEAP_ERR", "PRODUCT_ERR", "TOT_SIZE",
        "SPH_SIZE", "NUM_DSD", "DSD_SIZE", "NUM_DATA_SETS", "CRC",
        "ABS_ORBIT_START", "ABS_ORBIT_STOP", "L0_PROC_FLAG"]

    # A field that is in here is converted to an 32bit float
    _FLOAT_LIST = [
        "REL_TIME_ASC_NODE_START", "DELTA_UT1", "X_POSITION", "Y_POSITION",
        "Z_POSITION", "X_VELOCITY", "Y_VELOCITY", "Z_VELOCITY",
        "L0_PROCESSING_QUALITY", "L0_PROC_THRESH", "L0_GAPS_FLAG",
        "L0_GAPS_NUM", "OPEN_OCEAN_PERCENT", "CLOSE_SEA_PERCENT",
        "CONTINENT_ICE_PERCENT", "LAND_PERCENT", "L1B_PROD_STATUS",
        "L1B_PROC_FLAG", "L1B_PROCESSING_QUALITY", "L1B_PROC_THRESH",
        "REL_TIME_ASC_NODE_STOP", "EQUATOR_CROSS_LONG", "START_LAT",
        "START_LONG", "STOP_LAT", "STOP_LONG"]

    def __init__(self):
        super(CS2L1bMainProductHeader, self).__init__()

    def last_field(self, line):
        return re.search("CRC=", line)


class CS2L1bSpecificProductHeader(CS2L1bBaseHeader):
    """
    Container for the Specific Product Header of CryoSat-2 L1B science data
    """

    # Datatypes for automatic conversion
    # A field that is in here is converted to an 32bit integer
    _INT_LIST = ["ABS_ORBIT_START", "ABS_ORBIT_STOP", "L0_PROC_FLAG",
                 "L0_GAPS_FLAG", "L0_GAPS_NUM", "L1B_PROD_STATUS",
                 "L1B_PROC_FLAG", ]

    # A field that is in here is converted to an 32bit float
    _FLOAT_LIST = ["REL_TIME_ASC_NODE_START", "REL_TIME_ASC_NODE_STOP",
                   "EQUATOR_CROSS_LONG", "START_LAT", "START_LONG",
                   "STOP_LAT", "STOP_LONG", "L0_PROCESSING_QUALITY",
                   "L0_PROC_THRESH", "OPEN_OCEAN_PERCENT",
                   "CLOSE_SEA_PERCENT", "CONTINENT_ICE_PERCENT",
                   "LAND_PERCENT", "L1B_PROCESSING_QUALITY",
                   "L1B_PROC_THRESH"]

    def __init__(self):
        super(CS2L1bSpecificProductHeader, self).__init__()

    def last_field(self, line):
        return re.search("L1B_PROC_THRESH=", line)


def parse_cryosat_l1b_filename(filename):
    """
    Returns the information in the CryoSat-2 l1b filename
    """
    # Strip path and file extension
    filename = os.path.basename(filename)
    filename = os.path.splitext(filename)[0]
    # Construct the parser
    parser_str = "CS_{proc_stage}_"
    parser_str += "{instrument}_"
    parser_str += "{radar_mode}_"
    parser_str += "{data_level}_"
    parser_str += "{start_dt}_"
    parser_str += "{stop_dt}_"
    parser_str += "{baseline}"
    parser = parse.compile(parser_str)
    # Parse the filename
    result = parser.parse(filename)
    # Do some post-editing
    # - parse is not that smart when it comes to work with date strings
    # - naming conventions as the rest of pysiral
    info = TreeDict()
    info.mission = "cryosat2"
    info.instrument = result["instrument"].lower()
    info.radar_mode = result["radar_mode"].lower()
    info.data_level = "L"+result["data_level"]
    info.start_dt = dtparser.parse(result["start_dt"])
    info.stop_dt = dtparser.parse(result["stop_dt"])
    info.baseline = result["baseline"]
    return info


def parse_cryosat_l1b_xml_header(filename):
    """
    Reads the XML header file of a CryoSat-2 L1b Data set
    and returns the contents as an OrderedDict
    """
    with open(filename) as fd:
        content_odereddict = xmltodict.parse(fd.read())
    return content_odereddict[u'Earth_Explorer_Header']
