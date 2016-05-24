# -*- coding: utf-8 -*-
"""
Created on Sat Aug 01 17:33:02 2015

@author: Stefan
"""

from pysiral.config import ConfigInfo
from pysiral.output import NCDateNumDef
from netCDF4 import Dataset, num2date

import os
import tempfile
import uuid
import numpy as np


class ReadNC():
    """
    Quick & dirty method to parse content of netCDF file into a python object
    with attributes from file variables
    """
    def __init__(self, filename, verbose=False, autoscale=True):
        self.time_def = NCDateNumDef()
        self.parameters = []
        self.attributes = []
        self.verbose = verbose
        self.autoscale = autoscale
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
        for attribute_name in f.ncattrs():
            self.attributes.append(attribute_name)
            attribute_value = getattr(f, attribute_name)
            # Convert timestamps back to datetime objects
            if attribute_name in ["start_time", "stop_time"]:
                attribute_value = num2date(
                    attribute_value, self.time_def.units,
                    calendar=self.time_def.calendar)
            setattr(self, attribute_name, attribute_value)

        for key in f.variables.keys():
            setattr(self, key, np.array(f.variables[key][:]))
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


def get_l1bdata_files(mission_id, hemisphere, year, month, config=None):
    import glob
    if config is None:
        config = ConfigInfo()
    l1b_repo = config.local_machine.l1b_repository[mission_id].l1bdata
    directory = os.path.join(
        l1b_repo, hemisphere, "%04g" % year, "%02g" % month)
    print directory
    l1bdata_files = sorted(glob.glob(os.path.join(directory, "*.nc")))
    return l1bdata_files
