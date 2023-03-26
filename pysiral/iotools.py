# -*- coding: utf-8 -*-
"""
Created on Sat Aug 01 17:33:02 2015

@author: Stefan

TODO: Evaluate usefulness (or move to internal module)

"""

import tempfile
import uuid
from pathlib import Path

import numpy as np
from cftime import num2pydate
from dateperiods import DatePeriod
from netCDF4 import Dataset

from pysiral import psrlcfg
from pysiral.core.errorhandler import ErrorStatus
from pysiral.output import NCDateNumDef, PysiralOutputFilenaming


# TODO: Replace by xarray
class ReadNC(object):
    """
    Quick & dirty method to parse content of netCDF file into a python object
    with attributes from file variables
    """
    def __init__(self, filename, verbose=False, autoscale=True,
                 nan_fill_value=False, global_attrs_only=False):
        self.error = ErrorStatus()
        self.time_def = NCDateNumDef()
        self.keys = []
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

        # Open the file
        try:
            f = Dataset(self.filename)
            f.set_auto_scale(self.autoscale)
        except RuntimeError:
            msg = "Cannot read netCDF file: %s" % self.filename
            self.error.add_error("nc-runtime-error", msg)
            self.error.raise_on_error()

        # Try to update the time units
        # NOTE: This has become necessary with the use of
        #       variable epochs
        try:
            time_units = f.variables["time"].units
            self.time_def.units = time_units
        except (KeyError, AttributeError):
            pass

        # Get the global attributes
        for attribute_name in f.ncattrs():

            self.attributes.append(attribute_name)
            attribute_value = getattr(f, attribute_name)

            # Convert timestamps back to datetime objects
            # TODO: This needs to be handled better
            if attribute_name in ["start_time", "stop_time"]:
                attribute_value = num2pydate(attribute_value, self.time_def.units,
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
                    print(key)
            self.parameters = f.variables.keys()
        f.close()


class NCMaskedGridData(object):

    def __init__(self, filename, squeeze=True):
        self.squeeze = squeeze
        self.filename = filename
        self.parse()

    def parse(self):

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
    return Path(tempfile.gettempdir()) / str(uuid.uuid4())+".png"


def get_l1bdata_files(mission_id, hemisphere, year, month, config=None, version="default"):
    if config is None:
        config = psrlcfg
    l1b_repo = config.local_machine.l1b_repository[mission_id][version].l1bdata
    directory = Path(l1b_repo) / hemisphere / "%04g" % year / "%02g" % month
    l1bdata_files = sorted(directory.glob("*.nc"))
    return l1bdata_files




def l1bdata_get_baseline(filename):
    """ Returns version string in l1bdata filename  """
    # Parse infos from l1bdata filename
    fnattr = PysiralOutputFilenaming()
    fnattr.parse_filename(filename)
    return fnattr.version
