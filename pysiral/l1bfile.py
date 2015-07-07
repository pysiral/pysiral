# -*- coding: utf-8 -*-
"""
Created on Tue Jul 07 14:10:34 2015

@author: Stefan
"""

import os
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
        pass


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
