# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 13:51:52 2015

@author: Stefan
"""

from pysiral.errorhandler import FileIOErrorHandler
from pysiral.cryosat2.l1b_mds_def import cryosat2_get_mds_def
from pysiral.cryosat2.functions import (parse_cryosat_l1b_filename,
                                        parse_cryosat_l1b_xml_header,
                                        parse_cryosat2_l1b_header_field)
import os
import re
import numpy as np
import parse


class CryoSatL1B(object):
    """

    Purpose
    -------
    Container for CryoSat L1b Science Data including the xml header
    file (.HDR) and the header/data group file (.DBL). This class
    inherits all attributes and methods from **L1bData**.

    Applicable Documents
    --------------------
    - CRYOSAT Ground Segment, Instrument Processing Facility L1B,
      Products Specification Format, [PROD-FMT], CS-RS-ACS-GS-5106
      (issue used: 6.3, Date: 10/02/2015)

    Usage
    -----
    After initialisation of the class the filename of a .DBL file
    has to be supplied. Invoking the ``parse`` method will start
    the parsing of both the .HDR and .DBL files::

        l1b = CryoSatL1B()
        l1b.filename = filename
        l1b.parse()
        status = l1b.get_status()


    Header Attributes
    -----------------
    Available after ``parse`` method is called. Before that all are of type
    ``None``

    + xmlh
        Content of the xml header file (.HDR)
        The content is currently only saved as a nested OrderedDict
        object. This may change when the xml header information
        will be actually used.

    + mph
        Main Product Header block of the .DBL file. The attribute
        ``mph`` is of type **CS2L1bBaseHeader**.
        The mph fields consists of the parameters: name (tag), value,
        unit and datatype (dtype).

        Accessing fields:
            Each field can be either accessed directly. Example::

                orbit = l1b.mph.abs_orbit

            or with the `get_by_fieldname` method::

                orbit = l1b.mph.get_by_fieldname("abs_orbit")

            A list of all header field types is accessible in the attribute
            ``l1b.mph.field_list``.

        Units:
            The units are stored in a dictionary ``l1b.mph.unit_dict``
            where the keys are the field names::

                unit = l1b.mph.unit_dict["abs_orbit"]

            If no information on the unit is given in the MPH, the unit is
            of type ``None`` .

        Daty Types:
            The data type dictionary is mainly a debug parameter and
            contains the data types that have been used for converting
            the raw strings of the header into variables. Currently
            only types ``int32`` and ``float32`` are used for header
            number values. The default value is ``str``.
            The conversion is based on a specific list of data types for
            each header field. This list is located in the class
            definitions of **CS2L1bMainProductHeader** and
            **CS2L1bSpecificProductHeader** in ``pysiral/l1bfile.py``.

            The datatypes can be accessed like the units::

                dtype = l1b.mph.dtype_dict["abs_orbit"]

        Notes:
            - A quickview of the fields is given by
              ``print l1b.mph``
            - by convention all field names are lower case
            - No unit conversion is done as of yet if the unit contains
              a scaling (e.g. "10^-2%")

    + sph
        Specific Product Header block of the .DBL file. The attribute
        ``sph`` is of type `CS2L1bBaseHeader`. Accessing fields is
        identical to the main product header (see ``mph``).

    + dsd
        Data Set Descriptors in the .DBL file as an object of type
        **CS2L1bScienceDataSetDescriptors**.

        Each DSD in the dbl file is represented as an attribute of ``dsd`` that
        is of type ``dict``. The keys of the dictionary are given by the
        field names of each DSD, namely::

            ds_offset
            ds_tpye
            ds_size
            num_dsr
            filename

       A list of all DSD's is given in the string list ``l1b.dsd.field_list``.
       The information for each DSD can be accessed by directly calling the
       dictionary. E.g.::

           tide_file = l1b.dsd.pole_tide_file["filename"]



    Data Attributes
    ---------------


    Changelog
    ---------

    """

    _VALID_BASELINES = ["C001"]
    _VALID_RADAR_MODES = ["sar"]

    def __init__(self, raise_on_error=False):

        super(CryoSatL1B, self).__init__()
        # Error Handling
        self._init_error_handling(raise_on_error)
        self._baseline = None
        self._filename_header = None
        self._filename_product = None
        self.xmlh = None
        self.mph = None
        self.sph = None
        self.dsd = None

    @property
    def filename(self):
        return self._filename_product

    @filename.setter
    def filename(self, filename):
        """ Save and validate filenames for header and product file """
        # Test if valid file first
        self._error.file_undefined = not os.path.isfile(filename)
        if self._error.file_undefined:
            return
        # Split filenames in product and header file
        self._filename_product = filename
        self._filename_header = os.path.splitext(filename)[0]+".HDR"

    def parse(self):
        """ Parse the content of the L1B file """
        # Validate input and either return or raise when input not ok
        if self._error.test_errors():
            self._error.validate()
            return
        # First detect the baseline from the file name
        self._detect_filetype()
        # Parse the xml header file
        try:
            self._parse_header_file()
        except:
            self._error.io_failed = True
        # Parse the product file
        # XXX: Disabled during development
        self._parse_product_file()
#        try:
#            self._parse_product_file()
#        except:
#            self._error.io_failed = True
        # Validate the parsing
        self._error.validate()

    def get_status(self):
        pass

    def _init_error_handling(self, raise_on_error):
        self._error = FileIOErrorHandler()
        self._error.raise_on_error = raise_on_error
        self._error.file_undefined = True

    def _detect_filetype(self):
        """ Detect and validate the baseline """
        info = parse_cryosat_l1b_filename(self._filename_product)
        self.baseline = info.baseline
        self.radar_mode = info.radar_mode
        # Validate
        if self.baseline not in self._VALID_BASELINES:
            self._error.format_not_supported = True
        if self.radar_mode not in self._VALID_RADAR_MODES:
            self._error.format_not_supported = True

    def _parse_header_file(self):
        self._xmlh = parse_cryosat_l1b_xml_header(self._filename_header)

    def _parse_product_file(self):
        """
        Reads the content of *L1B*.DBL files
            - main product header (mph)
            - specific product header (sph)
            - data set descriptors (dsd)
            - (msd)
        """
        with open(self._filename_product, "r") as self._fh:
            self._parse_mph()
            self._parse_sph()
            self._parse_dsd()
        with open(self._filename_product, "rb") as self._fh:
            self._parse_msd()

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

    def _parse_dsd(self):
        """ Reads the Data Set Descriptors dsd's in the L1b header """
        self.dsd = CS2L1bScienceDataSetDescriptors()
        n_dsd_lines = self.dsd.get_num_lines(self.mph.num_dsd)
        for i in np.arange(n_dsd_lines+1):
            line = self._fh.readline()
            self.dsd.parse_line(line)

    def _parse_msd(self):
        """ Read the data blocks """
        # Just reopened the file in binary mode -
        # > get start byte and number of data set records
        l1b_data_set_name = "sir_l1b_"+self.radar_mode
        data_set_descriptor = self.dsd.get_by_fieldname(l1b_data_set_name)
        startbyte = int(data_set_descriptor["ds_offset"])
        self.n_msd_records = int(data_set_descriptor["num_dsr"])
        # Set the file pointer
        self._fh.seek(startbyte)
        # Get the parser
        mds_parser = cryosat2_get_mds_def(
            self.radar_mode, self.baseline, self.n_msd_records)
        # Parser the binary part of the .DBL file
        self.mds = mds_parser.parse(self._fh.read(mds_parser.sizeof()))

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


class CS2L1bScienceDataSetDescriptors(object):
    """
    Container for the Data Set Descriptors (dsd's) in the L1B file
    For each dsd there is an attribute with its name (dsd_name)
    that stores the information as a dictionary
    """

    _DSD_REFDICT = {
        "ds_type": "",
        "filename": "",
        "ds_offset": 0,
        "ds_size": 0,
        "num_dsr": 0,
        "dsr_size": 0
    }

    _DSD_UNITS = {
        "ds_type": None,
        "filename": None,
        "ds_offset": "bytes",
        "ds_size": "bytes",
        "num_dsr": None,
        "dsr_size": "bytes"
    }

    # There are 8 lines to read for each DSD records
    # 7 parameters + 1 empty line
    _N_RECORDS_PER_DSD = 8

    def __init__(self):
        self._current_dsd = None
        self.unit_dict = self._DSD_UNITS
        self.field_list = []

    def _is_dsd_start(self, tag):
        """ Checks if the current tag name indicated the start of a new dsd """
        return tag == "ds_name"

    def _append_field(self, tag, value):
        """ Adds a field to the current data set descriptor """
        tag = tag.lower()
        dsd_dict = getattr(self, self._current_dsd)
        dsd_dict[tag] = value
        setattr(self, self._current_dsd, dsd_dict)

    def get_num_lines(self, n_dsd):
        return n_dsd*self._N_RECORDS_PER_DSD

    def parse_line(self, line):
        """
        Parses a given line in the dsd block and addes the content
        to the dsd object
        """
        # Parse the line
        tag, value, unit = parse_cryosat2_l1b_header_field(line)
        # Skip empty line
        if tag is None:
            return
        # Check if start of new dsd and create a new one
        if self._is_dsd_start(tag):
            value = value.lower()
            setattr(self, value, self._DSD_REFDICT.copy())
            self._current_dsd = value
            self.field_list.append(value)
            return
        # save the field
        self._append_field(tag, value)

    def get_by_fieldname(self, fieldname):
        """
        Convinience function that returns the dict of an DSD records.

        Arguments:
            fieldname: (str):
                Name of the DSD field, case insensitive

        Output:
            value (dict):
                DSD dict

        """
        if fieldname not in self.field_list:
            return None
        return getattr(self, fieldname)

    def __repr__(self):
        """ String representation of the dsd's """
        output = ""
        for field in self.field_list:
            output += "DSD: %s\n" % (field)
            dsd_dict = getattr(self, field)
            for key in dsd_dict.keys():
                output += "  %s: %s, unit: %s\n" % (
                    key, str(dsd_dict[key]), str(self.unit_dict[key]))
        return output
