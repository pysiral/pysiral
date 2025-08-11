# -*- coding: utf-8 -*-

"""
This module provides access to deprecated legacy classes and functions that are still used in some parts of the codebase.
These classes and functions are maintained for backward compatibility but should not be used in new development.
"""

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"


import sys
from collections import OrderedDict, UserDict
from inspect import getframeinfo, stack
from pathlib import Path

from loguru import logger


class AttrDict(UserDict):
    """
    Short implementation of attrdict.AttrDict using UserDict. The code is based on the solutions shared here
    https://stackoverflow.com/a/76231823 and has been modified to allow nestesd AttrDict instances.
    """

    def __getattr__(self, key):
        item = self.__getitem__(key)
        return AttrDict(**item) if isinstance(item, dict) else item

    def __setattr__(self, key, value):
        if key == "data":
            return super().__setattr__(key, value)
        return self.__setitem__(key, value)


class ErrorStatus(object):

    def __init__(self, caller_id=""):
        self.caller_id = caller_id
        self.status = False
        self.codes = []
        self.messages = []
        self.reset()

    def add_error(self, code, message):
        """ Add an error. Error code and messages are arbitrary """
        self.status = True
        self.codes.append(code)
        self.messages.append(message)

    def raise_on_error(self):
        """ print error messages and exit program on existing error(s) """

        caller = getframeinfo(stack()[1][0])
        filename = Path(caller.filename).name
        if self.status:
            output = "{} Critical Error(s): {:g} [raised in {} L{}]\n"
            output = output.format(self.caller_id, len(self.codes), filename, caller.lineno)
            for i in range(len(self.codes)):
                output += f"  [{self.codes[i]}] {self.messages[i]}"
                output += "\n"
            logger.error(output)
            sys.exit(1)

    def get_all_messages(self):
        output = []
        if self.status:
            for i in range(len(self.codes)):
                error_message = f"{self.caller_id} error: [{self.codes[i]}] {self.messages[i]}"
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


class DefaultLoggingClass(object):
    """
    Template for default pysiral class with logging/error handling capabilities
    """

    def __init__(self, cls_name=None):
        """
        Init the class with a (loguru) logger and an ErrorStatus error handler
        :param cls_name:
        """
        self.error = ErrorStatus(cls_name)
        self.log = logger
