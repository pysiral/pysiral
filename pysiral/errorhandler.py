# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 15:25:45 2015

@author: Stefan
"""


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
        if type(value) is not 'bool':
            value = True
        self._raise_on_error = value

    def validate(self):
        """
        Check all error states and raise Exception when
        ``raise_on_error=True``
        """
        for error_name in self._error_dict.keys():
            if self._error_dict[error_name] and self._raise_on_error:
                raise self._exception_type


class FileIOErrorHandler(ErrorHandler):
    """
    Error Handler for reading files
    """
    def __init__(self):

        super(FileIOErrorHandler, self).__init__()
        self._exception_type = IOError
        self._error_dict = {
            "file_undefined": False,
            "io_failed": False}

    @property
    def file_undefined(self):
        return self._error_dict["file_undefined"]

    @file_undefined.setter
    def file_undefined(self, value):
        if type(value) is not 'bool':
            value = True
        self._error_dict["file_undefined"] = value

    @property
    def io_failed(self):
        return self._error_dict["io_failed"]

    @io_failed.setter
    def io_failed(self, value):
        if type(value) is not 'bool':
            value = True
        self._error_dict["io_failed"] = value
