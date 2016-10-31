# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 15:25:45 2015

@author: Stefan
"""
import sys
from collections import OrderedDict


PYSIRAL_ERROR_CODES = OrderedDict([
    ("auxdata_missing_sic", "Missing sea ice concentration data set"),
    ("auxdata_missing_sitype", "Missing ice type/MYI fraction data set"),
    ("auxdata_missing_snow", "Missing snow depth data set(s)"),
    ("l2proc_invalid_l1b", "Invalid l1bdata input data"),
    ("l2proc_surface_type_discarded", "Discarded by surface type validator")])


class ErrorStatus(object):

    def __init__(self, caller_id=""):
        self.caller_id = caller_id
        self.reset()

    def add_error(self, code, message):
        self.status = True
        self.codes.append(code)
        self.messages.append(message)

    def raise_on_error(self):
        if self.status:
            output = "%s Critical Error(s): (%g)\n" % (self.caller_id,
                                                       len(self.codes))
            for i in range(len(self.codes)):
                output += "  [%s] %s" % (self.codes[i], self.messages[i])
                output += "\n"
            print output
            sys.exit(1)

    def get_all_messages(self):
        output = []
        if self.status:
            for i in range(len(self.codes)):
                error_message = "%s error: [%s] %s" % (
                    self.caller_id, self.codes[i], self.messages[i])
                output.append(error_message)
        return output

    def reset(self):
        """ Remove all error messages and set to clean status """
        self.status = False
        self.codes = []
        self.messages = []

    @property
    def message(self):
        return ",".join(self.messages)


class ErrorHandler(object):
    """
    Parent class for all Errors (very early development phase)
    """
    def __init__(self):
        self._raise_on_error = False

    @property
    def raise_on_error(self):
        return self._error_dict["file_undefined"]

    @raise_on_error.setter
    def raise_on_error(self, value):
        if type(value) is not bool:
            value = True
        self._raise_on_error = value

    def validate(self):
        """
        Check all error states and raise Exception when
        ``raise_on_error=True``
        """
        for error_name in self._error_dict.keys():
            if self._error_dict[error_name] and self._raise_on_error:
                print self._error_dict, self._raise_on_error
                raise self._exception_type

    def test_errors(self):
        """ Returns True if any error is True, else False """
        output = False
        for error_name in self._error_dict.keys():
            if self._error_dict[error_name] and self._raise_on_error:
                output = True
        return output

    def _validate_flag(self, value):
        if type(value) is not bool:
            value = True
        return value


class FileIOErrorHandler(ErrorHandler):
    """
    Error Handler for reading files
    """
    def __init__(self):

        super(FileIOErrorHandler, self).__init__()
        self._exception_type = IOError
        self._error_dict = {
            "file_undefined": False,
            "io_failed": False,
            "format_not_supported": False}

    @property
    def file_undefined(self):
        return self._error_dict["file_undefined"]

    @file_undefined.setter
    def file_undefined(self, value):
        self._error_dict["file_undefined"] = self._validate_flag(value)

    @property
    def io_failed(self):
        return self._error_dict["io_failed"]

    @io_failed.setter
    def io_failed(self, value):
        self._error_dict["io_failed"] = self._validate_flag(value)

    @property
    def format_not_supported(self):
        return self._error_dict["format_not_supported"]

    @format_not_supported.setter
    def format_not_supported(self, value):
        self._error_dict["format_not_supported"] = self._validate_flag(value)
