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

from pysiral.logging import DefaultLoggingClass
from pysiral.errorhandler import ErrorStatus
from pysiral.helper import month_iterator, get_month_time_range

from datetime import datetime
from dateutil.relativedelta import relativedelta


import re
import os
import sys
import yaml
import socket
from treedict import TreeDict

import numpy as np

PYSIRAL_VERSION = "0.4.0-dev"
PYSIRAL_VERSION_FILENAME = "040dev"
HOSTNAME = socket.gethostname()


class ConfigInfo(DefaultLoggingClass):
    """
    Container for the content of the pysiral definition files
    (in pysiral/configration) and the local machine definition file
    (local_machine_definition.yaml)
    """

    # Global variables
    _DEFINITION_FILES = {
        "mission": "mission_def.yaml",
        "area": "area_def.yaml",
        "auxdata": "auxdata_def.yaml",
        "product": "product_def.yaml",
        "parameter": "parameter_def.yaml",
    }

    _LOCAL_MACHINE_DEF_FILE = "local_machine_def.yaml"

    VALID_DATA_LEVEL_IDS = ["l2", "l3"]

    def __init__(self):
        """ Read all definition files """
        super(ConfigInfo, self).__init__(self.__class__.__name__)

        self.error = ErrorStatus(self.__class__.__name__)
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
    def sandbox_path(self):
        """ Returns the local path to the document folder"""
        return self._return_path("sandbox")

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

    def get_mission_settings(self, mission):
        mission_options = self.mission[mission].settings
        return mission_options

    def get_mission_info(self, mission):
        mission_info = self.mission[mission]
        if mission_info.data_period.start is None:
            mission_info.data_period.start = datetime.utcnow()
        if mission_info.data_period.stop is None:
            mission_info.data_period.stop = datetime.utcnow()
        return mission_info

    def get_settings_file(self, data_level, setting_id_or_filename):
        """ Returns a processor settings file for a given data level.
        (data level: l2 or l3). The second argument can either be an
        direct filename (which validity will be checked) or an id, for
        which the corresponding file (id.yaml) will be looked up in
        the default directory """

        if data_level not in self.VALID_DATA_LEVEL_IDS:
            return None

        # Check if filename
        if os.path.isfile(setting_id_or_filename):
            return setting_id_or_filename

        # Get all settings files in settings/{data_level} and its
        # subdirectories
        lookup_directory = self.get_local_setting_path(data_level)
        ids, files = self.get_yaml_setting_filelist(lookup_directory)

        # Test if ids are unique and return error for the moment
        if len(np.unique(ids)) != len(ids):
            msg = "Non-unique %s setting filename" % data_level
            self.error.add_error("ambigous-setting-files", msg)
            self.error.raise_on_error()

        # Find filename to setting_id
        try:
            index = ids.index(setting_id_or_filename)
            return files[index]
        except:
            return None

    def get_yaml_setting_filelist(self, directory, ignore_obsolete=True):
        """ Retrieve all yaml files from a given directory (including
        subdirectories). Directories named "obsolete" are ignored if
        ignore_obsolete=True (default) """
        setting_ids = []
        setting_files = []
        for root, dirs, files in os.walk(directory):
            if os.path.split(root)[-1] == "obsolete" and ignore_obsolete:
                continue
            for filename in files:
                if re.search("yaml$", filename):
                    setting_ids.append(filename.replace(".yaml", ""))
                    setting_files.append(os.path.join(root, filename))
        return setting_ids, setting_files

    def get_local_setting_path(self, data_level):
        if data_level in self.VALID_DATA_LEVEL_IDS:
            return os.path.join(self.pysiral_local_path, "settings",
                                data_level)
        else:
            return None

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

    def get_name(self, flag):
        for mode_name, mode_flag in self.flag_dict.items():
            if flag == mode_flag:
                return mode_name
        return None

    def name(self, index):
        i = self.flag_dict.values().index(index)
        return self.flag_dict.keys()[i]

    @property
    def num(self):
        return len(self.flag_dict.keys())


class TimeRangeRequest(object):

    _PERIODS = ["monthly", "custom"]

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

    def raise_if_empty(self):
        message = ""
        if self._start_dt is None:
            message += "start time is invalid\n"
        if self._stop_dt is None:
            message += "stop time is invalid\n"

        if len(message) > 0:
            message += "Aborting ..."
            print message
            sys.exit(1)

    def set_range(self, start_date, stop_date):

        # Decode start and stop time definition

        # Check if datetime objects
        valid_start, valid_stop = False, False
        if isinstance(start_date, datetime):
            self._start_dt = start_date
            valid_start = True
        if isinstance(stop_date, datetime):
            self._stop_dt = stop_date
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
        """ Clip the current time range to an defined time range """

        is_clipped = False

        if self._start_dt < range_start and self._stop_dt > range_start:
            is_clipped = True
            self._start_dt = range_start
        elif self._start_dt < range_start and self._stop_dt < range_start:
            is_clipped = True
            self._start_dt = None
            self._stop_dt = None

        if self._stop_dt > range_stop and self._start_dt < range_stop:
            is_clipped = True
            self._stop_dt = range_stop
        elif self._stop_dt > range_stop and self._start_dt > range_stop:
            is_clipped = True
            self._start_dt = None
            self._stop_dt = None

        return is_clipped

    def set_period(self, period):
        """
        Set the period (monthly, weekly, etc) for the generation of
        iterations for the time range
        """

        if period in self._PERIODS:
            self._period = period
        else:
            raise ValueError("Invalid TimeRangeRequest period: %s" % period)

    def set_exclude_month(self, exclude_month_list):
        """
        Set a list of month, that shall be ignored during the generation of
        iterations for the time range
        """
        if exclude_month_list is None:
            exclude_month_list = []
        self._exclude_month = exclude_month_list

    def raise_on_error(self):

        if self.error.status:
            print self.error.message
            sys.exit(1)

    def get_iterations(self):
        """
        Return a list of iterations for the number of periods in the
        time range
        """

        iterations = []
        if self._start_dt is None or self._stop_dt is None:
            return iterations

        # monthly periods: return a list of time ranges that cover the full
        # month from the first to the last month
        if self._period == "monthly":

            # Get an iterator for integer year and month
            month_list = month_iterator(
                self._start_dt.year, self._start_dt.month,
                self._stop_dt.year, self._stop_dt.month)

            # Filter month that are excluded from processing
            month_list = [entry for entry in month_list if (
                entry[1] not in self._exclude_month)]

            # Create Iterations
            n_iterations = len(month_list)
            index = 1
            for year, month in month_list:

                # iteration will be a of type TimeRangeIteration
                time_range = TimeRangeIteration()

                # Per default get the full month
                month_start, month_stop = get_month_time_range(year, month)

                # limit the time range for first and last iteration
                # (only if the first and the last month are not in the
                #  exclude_month list)

                first_month = self._start_dt.month
                first_month_excluded = first_month in self._exclude_month
                if index == 1 and not first_month_excluded:
                    month_start = self._start_dt

                last_month = self._stop_dt.month
                last_month_excluded = last_month in self._exclude_month
                if index == n_iterations and not last_month_excluded:
                    month_stop = self.stop_dt

                # set final time range
                time_range.set_range(month_start, month_stop)
                time_range.set_indices(index, n_iterations)
                iterations.append(time_range)
                index += 1

        elif self._period == "custom":
            time_range = TimeRangeIteration(base_period="custom")
            time_range.set_range(self.start_dt, self.stop_dt)
            time_range.set_indices(1, 1)
            iterations = [time_range]

        return iterations

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
            self.error.add_error(self.__class__.__name__, error_message)
            return None

        # if stop time: add one period
        if start_or_stop == "stop":
            if n_entries == 2:
                extra_period = relativedelta(months=1, microseconds=-1)
            else:
                extra_period = relativedelta(days=1, microseconds=-1)
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

    @property
    def label(self):
        return str(self.start_dt)+" till "+str(self.stop_dt)


class TimeRangeIteration(object):

    def __init__(self, base_period="monthly"):
        self._index = 0
        self._num_iterations = 0
        self._start = None
        self._stop = None
        self._base_period = base_period

    def __repr__(self):
        output = "pysiral TimeRangeIteration Object:\n"
        output += "%6s: %s\n" % ("start", str(self.start))
        output += "%6s: %s\n" % ("stop", str(self.stop))
        output += "Iteration: %s of %s\n" % (
            str(self.index), str(self.num_iterations))
        return output

    def set_range(self, start, stop):
        self._start = start
        self._stop = stop

    def set_indices(self, i, n):
        self._index = i
        self._num_iterations = n

    @property
    def start(self):
        return self._start

    @property
    def stop(self):
        return self._stop

    @property
    def index(self):
        return self._index

    @property
    def num_iterations(self):
        return self._num_iterations

    @property
    def is_full_month(self):
        test_diff = (self.start+relativedelta(months=1))-self.stop
        return test_diff.days < 1

    @property
    def base_period(self):
        return self._base_period

    @property
    def label(self):
        return str(self.start)+" till "+str(self.stop)


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
                "help": 'list of months to be excluded from processing'},

            # Flag that indicates if previous versions shall be removed
            # before processing / plotting etc.
            "remove-old": {
                "action": "store_true",
                "dest": "remove_old",
                "default": False,
                "required": False,
                "help": 'remove all existing product in target directory'},

            # version tag of input data
            "input-version": {
                "action": "store",
                "dest": "input_version",
                "default": "default",
                "required": False,
                "help": 'input version name (see documentation)'},

            # override any critical prompts for cronjobs etc
            "no-critical-prompt": {
                "action": "store_true",
                "dest": "no_critical_prompt",
                "default": False,
                "required": False,
                "help": 'set to skip any required command line inputs'},

            # preset for level-1b (l1bdata) fiels
            "l1b_files": {
                "action": "store",
                "dest": "l1b_files_preset",
                "default": None,
                "required": False,
                "help": 'Path to one or many l1bdata files (e.g.: path/*.nc)'},

            # fetch the level-2 settings file
            "l2-settings": {
                "action": "store",
                "dest": "l2_settings",
                "default": None,
                "required": True,
                "help": 'id or path to Level-2 settings file'},

            # set the run tag for the Level-2 Processor
            "run-tag": {
                "action": "store",
                "dest": "run_tag",
                "default": None,
                "required": False,
                "help": 'tag for the Level-2 output'},

            # no overwrite protection for level-2 outputs
            "no-overwrite-protection": {
                "action": "store_false",
                "dest": "overwrite_protection",
                "default": False,
                "required": False,
                "help": 'disable writing Level-2 output to unique directory'},

            # no overwrite protection for level-2 outputs
            "overwrite-protection": {
                "action": "store_true",
                "dest": "overwrite_protection",
                "default": True,
                "required": False,
                "help": 'enable writing Level-2 output to unique directory ' +
                        '(default)'}}

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

def get_parameter_attributes(target):
    config = ConfigInfo()
    return config.parameter[target]
