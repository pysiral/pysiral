# -*- coding: utf-8 -*-
"""
Created on Sat Aug 01 17:33:02 2015

@author: Stefan
"""

from netCDF4 import Dataset, date2num
from datetime import datetime

import numpy as np
import os


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


class NCDateNumDef(object):
    """
    Holds definition for datetime conversion to numbers and vice versa
    for netCDF operations
    """

    def __init__(self):
        self.units = "seconds since 1970-01-01"
        self.calendar = "standard"


class L1bNCfile(object):
    """
    Class to export a L1bdata object into a netcdf file
    """

    def __init__(self):
        self.datagroups = []
        self.output_folder = None
        self.l1b = None
        self.filename = None
        self.zlib = False
        self.time_def = NCDateNumDef()
        self._rootgrp = None

    def export(self):
        self._validate()
        self._create_filename()
        self._create_root_group()
        self._populate_data_groups()
        self._write_to_file()

    def _validate(self):
        pass

    def _create_filename(self):
        basename = file_basename(self.l1b.filename)
        self.filename = os.path.join(self.output_folder, basename+".nc")

    def _create_root_group(self):
        """
        Create the root group and add l1b metadata as global attributes
        """
        self._rootgrp = Dataset(self.filename, "w")
        # Save the l1b info data group as global attributes
        attdict = self.l1b.info.attdict
        self._convert_datetime_attributes(attdict)
        self._convert_bool_attributes(attdict)
        self._convert_nonetype_attributes(attdict)
        self._set_global_attributes(attdict)

    def _populate_data_groups(self):
        for datagroup in self.datagroups:
            print datagroup.upper()
            # Create the datagroup
            dgroup = self._rootgrp.createGroup(datagroup)
            content = getattr(self.l1b, datagroup)
            # Create the dimensions
            # (must be available as OrderedDict in Datagroup Container
            dims = content.dimdict.keys()
            for key in dims:
                dgroup.createDimension(key, content.dimdict[key])
            # Now add variables for each parameter in datagroup
            for parameter in content.parameter_list:
                data = getattr(content, parameter)
                # Convert datetime objects to number
                if type(data[0]) is datetime:
                    data = date2num(data, self.time_def.units,
                                    self.time_def.calendar)
                # Convert bool objects to integer
                if data.dtype.str == "|b1":
                    data = np.int8(data)
                dimensions = tuple(dims[0:len(data.shape)])
                print " "+parameter, dimensions, data.dtype.str
                var = dgroup.createVariable(
                    parameter, data.dtype.str, dimensions, zlib=self.zlib)
                var[:] = data

    def _write_to_file(self):
        self._rootgrp.close()

    def _convert_datetime_attributes(self, attdict):
        """
        Replace l1b info parameters of type datetime.datetime by a double
        representation to match requirements for netCDF attribute data type
        rules
        """
        for key in attdict.keys():
            content = attdict[key]
            if type(content) is datetime:
                attdict[key] = date2num(
                    content, self.time_def.units, self.time_def.calendar)

    def _convert_bool_attributes(self, attdict):
        """
        Replace l1b info parameters of type bool ['b1'] by a integer
        representation to match requirements for netCDF attribute data type
        rules
        """
        for key in attdict.keys():
            content = attdict[key]
            if type(content) is bool:
                attdict[key] = int(content)

    def _convert_nonetype_attributes(self, attdict):
        """
        Replace l1b info parameters of type bool ['b1'] by a integer
        representation to match requirements for netCDF attribute data type
        rules
        """
        for key in attdict.keys():
            content = attdict[key]
            if content is None:
                attdict[key] = ""

    def _set_global_attributes(self, attdict):
        """ Save l1b.info dictionary as global attributes """
        for key in attdict.keys():
            print key, type(attdict[key])
            self._rootgrp.setncattr(key, attdict[key])


def file_basename(filename, fullpath=False):
    """
    Returns the filename without file extension of a give filename (or path)
    """
    strarr = os.path.split(filename)
    file_name = strarr[-1]
    basename = file_name.split(".")[0]
    if fullpath:
        basename = os.path.join(strarr[0], basename)
    return basename
