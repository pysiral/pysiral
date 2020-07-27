# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 15:25:45 2015

@author: Stefan
"""
import sys
from loguru import logger
from pathlib import Path
from inspect import getframeinfo, stack
from collections import OrderedDict


# TODO: This is also obsolete
PYSIRAL_ERROR_CODES = OrderedDict([
    ("auxdata_invalid_class", "Invalid auxdata class name [%s]"),
    ("auxdata_invalid_class_name", "Auxdata class does not exist [%s]"),
    ("auxdata_missing_definition", "auxdata_def: Missing definition [%s:%s]"),
    ("auxdata_missing_localrepo_def", "local_machine_def: Missing definition [%s:%s]"),
    ("auxdata_missing_sic", "Missing sea ice concentration data set"),
    ("auxdata_missing_sitype", "Missing ice type/MYI fraction data set"),
    ("auxdata_missing_snow", "Missing snow depth data set(s)"),
    ("l2proc_invalid_l1b", "Invalid l1bdata input data"),
    ("l2proc_surface_type_discarded", "Discarded by surface type validator"),
    ("warren99-invalid-hemisphere", "Warren99 snow climatology only valid for northern hemisphere")])


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
                output += "  [%s] %s" % (self.codes[i], self.messages[i])
                output += "\n"
            logger.error(output)
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
