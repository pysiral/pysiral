# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 13:19:45 2016

@author: Stefan
"""

from pysiral.config import options_from_dictionary


class ValidatorBaseClass(object):

    def __init__(self):
        pass

    def set_options(self, **opt_dict):
        self._options = options_from_dictionary(**opt_dict)

    def validate(self, *args, **kwargs):
        error_flag, error_message = self._validate(*args, **kwargs)
        return error_flag, error_message


# %% Surface Type Validators

class LeadFixedMinimumNumber(ValidatorBaseClass):

    """
    Validate the number of leads of a surface type classification.
    -> Invalid if number of leads falls below fixed number independent of
       length of profile
    """

    def __init__(self):
        super(LeadFixedMinimumNumber, self).__init__()

    def _validate(self, l2):
        error_status = False
        error_message = ""
        if l2.surface_type.lead.num < self._options.minimum_n_leads:
            error_status = True
            error_message = "%s: lead count (%g) below threshold (%g)" % (
                type(self).__name__, l2.surface_type.lead.num,
                self._options.minimum_n_leads)
        return error_status, error_message


# %% Public Functions

def get_validator(name):
    return globals()[name]()
