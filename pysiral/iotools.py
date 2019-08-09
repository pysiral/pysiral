# -*- coding: utf-8 -*-
"""
Created on Sat Aug 01 17:33:02 2015

@author: Stefan
"""

from pysiral.config import ConfigInfo, TimeRangeIteration
from pysiral.errorhandler import ErrorStatus
from pysiral.output import NCDateNumDef, PysiralOutputFilenaming
from pysiral.path import file_basename
from netCDF4 import Dataset, num2date

import os
import glob
import tempfile
import uuid
import numpy as np


class ReadNC(object):
    """
    Quick & dirty method to parse content of netCDF file into a python object
    with attributes from file variables
    """
    def __init__(self, filename, verbose=False, autoscale=True,
                 nan_fill_value=False, global_attrs_only=False):
        self.error = ErrorStatus()
        self.time_def = NCDateNumDef()
        self.parameters = []
        self.attributes = []
        self.verbose = verbose
        self.autoscale = autoscale
        self.global_attrs_only = global_attrs_only
        self.nan_fill_value = nan_fill_value
        self.filename = filename
        self.parameters = []
        self.read_globals()
        self.read_content()

    def read_globals(self):
        pass
#        self.gobal_attributes = {}
#        f = Dataset(self.filename)
#        print f.ncattrs()
#        f.close()

    def read_content(self):

        self.keys = []

        # Open the file
        try:
            f = Dataset(self.filename)
        except RuntimeError:
            msg = "Cannot read netCDF file: %s" % self.filename
            self.error.add_error("nc-runtime-error", msg)
            self.error.raise_on_error()

        f.set_auto_scale(self.autoscale)

        # Get the global attributes
        for attribute_name in f.ncattrs():

            self.attributes.append(attribute_name)
            attribute_value = getattr(f, attribute_name)

            # Convert timestamps back to datetime objects
            # TODO: This needs to be handled better
            if attribute_name in ["start_time", "stop_time"]:
                attribute_value = num2date(
                    attribute_value, self.time_def.units,
                    calendar=self.time_def.calendar)
            setattr(self, attribute_name, attribute_value)

        # Get the variables
        if not self.global_attrs_only:
            for key in f.variables.keys():

                try:
                    variable = f.variables[key][:]
                except ValueError:
                    continue

                try:
                    is_float = variable.dtype in ["float32", "float64"]
                    has_mask = hasattr(variable, "mask")
                except:
                    is_float, has_mask = False, False

                if self.nan_fill_value and has_mask and is_float:
                    is_fill_value = np.where(variable.mask)
                    variable[is_fill_value] = np.nan

                setattr(self, key, variable)
                self.keys.append(key)
                self.parameters.append(key)
                if self.verbose:
                    print key
            self.parameters = f.variables.keys()
        f.close()


class NCMaskedGridData(object):

    def __init__(self, filename, squeeze=True):
        self.squeeze = squeeze
        self.filename = filename
        self.parse()

    def parse(self):
        from pysiral.iotools import ReadNC

        nc = ReadNC(self.filename)

        self.parameters = nc.parameters
        for parameter in nc.parameters:
            data = np.ma.array(getattr(nc, parameter))
            if self.squeeze:
                data = np.squeeze(data)
            data.mask = np.isnan(data)
            setattr(self, parameter, data)

        self.attributes = nc.attributes
        for attribute in nc.attributes:
            setattr(self, attribute, getattr(nc, attribute))

    def get_by_name(self, parameter_name):
        try:
            return getattr(self, parameter_name)
        except:
            return None


def get_temp_png_filename():
    return os.path.join(tempfile.gettempdir(), str(uuid.uuid4())+".png")


def get_l1bdata_files(mission_id, hemisphere, year, month, config=None,
                      version="default"):
    import glob
    if config is None:
        config = ConfigInfo()
    l1b_repo = config.local_machine.l1b_repository[mission_id][version].l1bdata
    directory = os.path.join(
        l1b_repo, hemisphere, "%04g" % year, "%02g" % month)
    l1bdata_files = sorted(glob.glob(os.path.join(directory, "*.nc")))
    return l1bdata_files


def get_local_l1bdata_files(mission_id, time_range, hemisphere, config=None,
                            version="default", allow_multiple_baselines=True):
    """
    Returns a list of l1bdata files for a given mission, hemisphere, version
    and time range
    XXX: Note: this function will slowly replace `get_l1bdata_files`, which
         is limited to full month
    """

    # parse config data (if not provided)
    if config is None or not isinstance(config, ConfigInfo):
        config = ConfigInfo()

    # Validate time_range (needs to be of type TimeRangeIteration)
    try:
        time_range_is_correct_object = time_range.base_period == "monthly"
    except:
        time_range_is_correct_object = False
    if not time_range_is_correct_object:
        error = ErrorStatus()
        msg = "Invalid type of time_range, required: %s, was %s" % (
            type(time_range), type(TimeRangeIteration))
        error.add_error("invalid-timerange-type", msg)
        error.raise_on_error()

    # 1) get list of all files for monthly folders
    yyyy, mm = "%04g" % time_range.start.year, "%02g" % time_range.start.month

    # NOTE: Check if folder is l1p instead of l1bdata
    #       -> l1p indicated new l1 data version with different filenaming. File structure is still the same
    repo_branch = config.local_machine.l1b_repository[mission_id][version]
    if repo_branch.has_key("l1bdata"):
        l1b_repo = repo_branch.l1bdata
    else:
        l1b_repo = repo_branch.l1p
    directory = os.path.join(l1b_repo, hemisphere, yyyy, mm)
    all_l1bdata_files = sorted(glob.glob(os.path.join(directory, "*.nc")))

    # 2) First filtering step: Check if different algorithm baseline values
    # exist in the list of l1bdata files
    algorithm_baselines = [l1bdata_get_baseline(f) for f in all_l1bdata_files]
    baselines = np.unique(np.array(algorithm_baselines))
    n_baselines = len(baselines)
    if not allow_multiple_baselines and n_baselines > 1:
        error = ErrorStatus()
        baseline_str_list = ", ".join(baselines)
        msg = "Multiple l1bdata baselines (%g) [%s] found in directory: %s" % (
                n_baselines, baseline_str_list, directory)
        error.add_error("multiple-l1b-baselines", msg)
        error.raise_on_error()

    # 3) Check if files are in requested time range
    # This serves two purporses: a) filter out files with timestamps that do
    # not belong in the directory. b) get a subset if required
    l1bdata_files_checked = [l1bdata_file for l1bdata_file in all_l1bdata_files
                             if l1bdata_in_trange(l1bdata_file, time_range)]

    # Done return list (empty or not)
    return l1bdata_files_checked, directory


def l1bdata_in_trange(fn, tr):
    """ Returns flag if filename is within time range """
    # Parse infos from l1bdata filename
    fnattr = PysiralOutputFilenaming()
    fnattr.parse_filename(fn)
    # Compute overlap between two start/stop pairs
    is_overlap = fnattr.start <= tr.stop and fnattr.stop >= tr.start
    return is_overlap


def l1bdata_get_baseline(filename):
    """ Returns version string in l1bdata filename  """
    # Parse infos from l1bdata filename
    fnattr = PysiralOutputFilenaming()
    fnattr.parse_filename(filename)
    return fnattr.version
