# -*- coding: utf-8 -*-
"""
Created on Sat Aug 01 17:33:02 2015

@author: Stefan
"""

from netCDF4 import Dataset

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
        self.parameters = []
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
        for parameter in nc.parameters:
            data = np.ma.array(getattr(nc, parameter))
            data.mask = np.isnan(data)
            setattr(self, parameter, data)

    def get_by_name(self, parameter_name):
        try:
            return getattr(self, parameter_name)
        except:
            return None


def get_temp_png_filename():
    return os.path.join(tempfile.tempdir, str(uuid.uuid4())+".png")
