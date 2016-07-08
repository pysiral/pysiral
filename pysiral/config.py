# -*- coding: utf-8 -*-
#
# Copyright © 2015 Stefan Hendricks
#
# Licensed under the terms of the GNU GENERAL PUBLIC LICENSE
#
# (see LICENSE for details)

"""
Purpose:
    Returns content of configuration and definition files

Created on Mon Jul 06 10:38:41 2015

@author: Stefan
"""

from pysiral.errorhandler import ErrorStatus

from datetime import datetime
from dateutil.relativedelta import relativedelta


import os
import sys
import yaml
from treedict import TreeDict


class ConfigInfo(object):
    """
    Container for the content of the pysiral definition files
    (in pysiral/configration) and the local machine definition file
    (local_machine_definition.yaml)
    """

    PYSIRAL_VERSION = "0.2.pre"

    # Global variables
    _DEFINITION_FILES = {
        "mission": "mission_def.yaml",
        "area": "area_def.yaml",
        "auxdata": "auxdata_def.yaml",
        "product": "product_def.yaml",
        "parameter": "parameter_def.yaml",
    }

    _LOCAL_MACHINE_DEF_FILE = "local_machine_def.yaml"

    def __init__(self):
        """ Read all definition files """
        # Store the main path on this machine
        self.pysiral_local_path = get_pysiral_local_path()
        # read the definition files in the config folder
        self._read_config_files()
        # read the local machine definition file
        self._read_local_machine_file()

    @property
    def doc_path(self):
        """ Returns the local path to the document folder"""
        return self._return_path("doc")

    @property
    def config_path(self):
        """ Returns the local path to the document folder"""
        return self._return_path("config")

    @property
    def mission_ids(self):
        return self.mission.missions

    def get_mission_defaults(self, mission):
        mission_options = self.mission[mission].options
        defaults = {}
        names, options = td_branches(mission_options)
        for name, option in zip(names, options):
            defaults[name] = option.default
        return defaults

    def get_mission_options(self, mission):
        mission_options = self.mission[mission].options
        return mission_options

    def get_mission_info(self, mission):
        mission_info = self.mission[mission]
        if mission_info.data_period.start is None:
            mission_info.data_period.start = datetime.utcnow()
        if mission_info.data_period.stop is None:
            mission_info.data_period.stop = datetime.utcnow()
        return mission_info

    def _read_config_files(self):
        for key in self._DEFINITION_FILES.keys():
            content = get_yaml_config(
                os.path.join(
                    self.config_path,
                    self._DEFINITION_FILES[key]))
            setattr(self, key, content)

    def _read_local_machine_file(self):
        content = get_yaml_config(
            os.path.join(
                self.pysiral_local_path, self._LOCAL_MACHINE_DEF_FILE))
        setattr(self, "local_machine", content)

    def _return_path(self, subfolder):
        return os.path.join(self.pysiral_local_path, subfolder)


class RadarModes(object):

    flag_dict = {"lrm": 0, "sar": 1, "sin": 2}

    def __init__(self):
        pass

    def get_flag(self, mode_name):
        try:
            return self.flag_dict[mode_name]
        except:
            return None

    def name(self, index):
        i = self.flag_dict.values().index(index)
        return self.flag_dict.keys()[i]

    @property
    def num(self):
        return len(self.flag_dict.keys())


class TimeRangeRequest(object):

    _PERIODS = ["monthly"]

    def __init__(self):
        self._start_dt = None
        self._stop_dt = None
        self._period = self._default_period
        self._exclude_month = []
        self.error = ErrorStatus()

    def __repr__(self):
        output = "TimeRangeRequest object:\n"
        for field in ["_start_dt", "_stop_dt", "_period", "_exclude_month"]:
            output += "%12s: %s" % (field, getattr(self, field))
            output += "\n"
        return output

    def set_range(self, start_date, stop_date):

        # Decode start and stop time definition

        # Check if datetime objects
        valid_start, valid_stop = False, False
        if isinstance(start_date, datetime):
            self._start_dt = start_date
            valid_start = True
        if isinstance(stop_date, datetime):
            self._stopt_dt = stop_date
            valid_stop = True

        if valid_start and valid_stop:
            return

        # Check and decode integer lists
        if isinstance(start_date, list):
            if all(isinstance(item, int) for item in start_date):
                self._start_dt = self._decode_int_list(start_date, "start")
            else:
                error_message = "invalid start time (not integer list)"
                self.error.append(self.__class__.__name__, error_message)

        if isinstance(stop_date, list):
            if all(isinstance(item, int) for item in stop_date):
                self._stop_dt = self._decode_int_list(stop_date, "stop")
            else:
                error_message = "invalid stop time (non integer list)"
                self.error.append(self.__class__.__name__, error_message)

    def clip_to_range(self, range_start, range_stop):

        if self._start_dt < range_start and self._stop_dt > range_start:
            self._start_dt = range_start
        elif self._start_dt < range_start and self._stop_dt < range_start:
            self._start_dt = None
            self._stop_dt = None

        if self._stop_dt > range_stop and self._start_dt < range_stop:
            self._stop_dt = range_stop
        elif self._stop_dt > range_stop and self._start_dt > range_stop:
            self._start_dt = None
            self._stop_dt = None

    def set_period(self, period):
        if period in self._PERIODS:
            self._period = period
        else:
            raise ValueError("Invalid TimeRangeRequest period: %s" % period)

    def raise_on_error(self):
        if self.error.status:
            print self.error.message
            sys.exit(1)

    def _decode_int_list(self, int_list, start_or_stop):

        # XXX: Currently only yyyy mm [dd] (day is optional) are allowed
        n_entries = len(int_list)
        if n_entries < 2 or n_entries > 3:
            error_message = "%s date integer list must be yyyy mm [dd]"
            self.error.append(self.__class__.__name__, error_message)
            return None

        # Set the day
        day = 1 if n_entries == 2 else int_list[2]

        # Set the datetime object (as if would be start date)
        # Raise error and return none if unsuccessful
        try:
            dt = datetime(int_list[0], int_list[1], day)
        except:
            error_message = "cannot convert integer list to datetime: %s" % (
                str(int_list))
            self.error.append(self.__class__.__name__, error_message)
            return None

        # if stop time: add one period
        if start_or_stop == "stop":
            if n_entries == 2:
                extra_period = relativedelta(months=1, seconds=-1)
            else:
                extra_period = relativedelta(days=1, seconds=-1)
            dt = dt + extra_period

        return dt

    @property
    def _default_period(self):
        return self._PERIODS[0]

    @property
    def start_dt(self):
        return self._start_dt

    @property
    def stop_dt(self):
        return self._stop_dt


class DefaultCommandLineArguments(object):

    def __init__(self):

        config = ConfigInfo()

        self._args = {

            # Mission id
            "mission": {
                "action": 'store',
                "dest": 'mission_id',
                "choices": config.mission_ids,
                "required": True,
                "help": "pysiral recognized mission id"},

            # Default date parameter
            "date": {
                "action": "store",
                "dest": "stop_date",
                "nargs": "+",
                "type": int,
                "required": False,
                "help": 'list as year and month and day (optional)'},

            # Default date parameter
            "hemisphere": {
                "action": "store",
                "dest": "hemisphere",
                "choices": ["global", "north", "south"],
                "default": "global",
                "required": False,
                "help": 'hemisphere flag for processing)'},

            # List of month to exclude from monthly iterations
            "exclude-month": {
                "action": "store",
                "dest": "exclude_month",
                "nargs": "+",
                "type": int,
                "required": False,
                "default": None,
                "help": 'list of months to be excluded from processing'}}

    def get_argparse_dict(self, name, destination, required):
        options = self._args[name]
        options["dest"] = destination
        options["required"] = required
        return options


def get_yaml_config(filename, output="treedict"):
    """
    Parses the contents of a configuration file in .yaml format
    and returns the content in various formats

    Arguments:
        filename (str)
            path the configuration file

    Keywords:
        output (str)
            "treedict" (default): Returns a treedict object
            "dict": Returns a python dictionary
    """
    with open(filename, 'r') as f:
        content_dict = yaml.load(f)

    if output == "treedict":
        return TreeDict.fromdict(content_dict, expand_nested=True)
    else:
        return content_dict


def get_pysiral_local_path():
    """ Returns pysiral's main directory for the local machine """
    directory = os.path.dirname(__file__)
    directory = os.path.abspath(os.path.join(directory, '..'))
    return directory


def get_parameter_definitions():
    config = ConfigInfo()
    return config.parameter.definition


def td_branches(t):
    """ Convinience function to get only the branches of a treedict object """
    try:
        branch_names = list(t.iterkeys(recursive=False, branch_mode='only'))
        branch_objects = list(t.iterbranches())
    except:
        branch_names = []
        branch_objects = []
    return branch_names, branch_objects


def options_from_dictionary(**opt_dict):
    """ Function for converting option dictionary in Treedict """
    return TreeDict.fromdict(opt_dict, expand_nested=True)
