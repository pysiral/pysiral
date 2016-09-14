# -*- coding: utf-8 -*-
"""
Created on Sat Aug 01 17:33:02 2015

@author: Stefan
"""

from pysiral.config import ConfigInfo, TimeRangeIteration
from pysiral.errorhandler import ErrorStatus
from pysiral.output import NCDateNumDef, PysiralOutputFileNaming
from pysiral.path import file_basename
from netCDF4 import Dataset, num2date

import os
import glob
import tempfile
import uuid
import numpy as np


class ReadNC():
    """
    Quick & dirty method to parse content of netCDF file into a python object
    with attributes from file variables
    """
    def __init__(self, filename, verbose=False, autoscale=True,
                 nan_fill_value=False):
        self.time_def = NCDateNumDef()
        self.parameters = []
        self.attributes = []
        self.verbose = verbose
        self.autoscale = autoscale
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
        f = Dataset(self.filename)
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
        for key in f.variables.keys():

            variable = f.variables[key][:]

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

    def __init__(self, filename):
        self.filename = filename
        self.parse()

    def parse(self):
        from pysiral.iotools import ReadNC

        nc = ReadNC(self.filename)

        self.parameters = nc.parameters
        for parameter in nc.parameters:
            data = np.ma.array(getattr(nc, parameter))
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
                      version="default"):
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
        error.add("invalid-timerange-type", msg)
        error.raise_on_error()

    # 1) get list of full month
    yyyy, mm = "%04g" % time_range.start.year, "%02g" % time_range.start.month
    l1b_repo = config.local_machine.l1b_repository[mission_id][version].l1bdata
    directory = os.path.join(l1b_repo, hemisphere, yyyy, mm)
    l1bdata_files = sorted(glob.glob(os.path.join(directory, "*.nc")))

    # 2) check if month subset is requested
    tr = time_range
    if time_range.is_full_month:
        return l1bdata_files
    else:
        subset = [f for f in l1bdata_files if l1bdata_in_trange(f, tr)]
        return subset


def l1bdata_in_trange(fn, tr):
    """ Returns flag if filename is within time range """
    # Parse infos from l1bdata filename
    fnattr = PysiralOutputFileNaming()
    fnattr.parse_filename(fn)
    # Compute overlap between two start/stop pairs
    is_overlap = fnattr.start <= tr.stop and fnattr.stop >= tr.start
    return is_overlap
