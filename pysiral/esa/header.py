# -*- coding: utf-8 -*-

import parse
import numpy as np


class ESAProductHeader(object):
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


class ESAScienceDataSetDescriptors(object):
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

    _VALID_TAGS = ["ds_name"]
    _VALID_TAGS.extend(_DSD_REFDICT.keys())

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
        tag, value, unit = parse_header_field(line)
        # Skip empty line

        if tag not in self._VALID_TAGS:
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


def parse_header_field(line):
    """
    Parses a line of the CryoSat-2 L1B ascii header and
    returns tag, value and unit as strings.

    Note: double quotes in string items are removed

    e.g.

    NUM_DSD=+0000000044
    -> "num_dsd", "+0000000044", None

    DSD_SIZE=+0000000280<bytes>
    -> "dsd_size", "+0000000280", "bytes"
    """
    tag, value, unit = None, None, None
    parser = parse.compile("{tag}={value}")
    match = parser.parse(line)
    if match:
        tag = match["tag"].lower()
        value = match["value"]
        value = value.replace("\"", "").strip()
        # check if unit is given
        unit_parser = parse.compile("{value}<{unit}>")
        unit_match = unit_parser.parse(value)
        if unit_match:
            value = unit_match["value"]
            unit = unit_match["unit"]
    return tag, value, unit
