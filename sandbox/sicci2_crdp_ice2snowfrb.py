# -*- coding: utf-8 -*-
"""
Takes l2i repository and creates l2i files changes freeboard and thickness
computation from ice freeboard assumption to snow freeboard assumption

(e.g. fast tracking the computation from ccicdr-southern-*-P002a
to ccicdr-southern-*-P002a)

Created on Fri Feb 17 13:19:28 2017

@author: shendric
"""

from pysiral.logging import DefaultLoggingClass
from pysiral.errorhandler import ErrorStatus
from pysiral.l2data import L2iNCFileImport
from pysiral.path import validate_directory, filename_from_path
from pysiral.sit import snowfreeboard2thickness
from pysiral.output import NCDataFile

from pysiral.config import PYSIRAL_VERSION
from pysiral.config import get_parameter_attributes


from netCDF4 import date2num
from datetime import datetime

from collections import OrderedDict

import numpy as np
import argparse
import glob
import sys
import os
import re


class CRDPConvertP002aToP002b(DefaultLoggingClass):

    def __init__(self):
        super(CRDPConvertP002aToP002b, self).__init__(self.__class__.__name__)
        self.config = CRDPConverterConfig()
        self.l2irepo = CRDPConverterL2iList()

    def collect_l2i_files(self):
        self.l2irepo.get_file_stack(self.config.l2i_source_repo)

    def convert_files(self):
        for year, month, l2i_files in self.l2irepo.iteration:
            self.log.info("+ Convert l2i files %s-%s" % (year, month))
            output_folder = self.config.get_output_folder(year, month)
            self.log.info("- Output folder: %s" % output_folder)
            validate_directory(output_folder)
            for l2i_file in l2i_files:
                self.convert(l2i_file, output_folder)

    def convert(self, l2i_file, output_folder):

        # Read the data file
        l2i = L2iNCFileImport(l2i_file)

        # step one: Replace freeboard by radar freeboard
        l2i.freeboard = l2i.radar_freeboard

        # filter freeboard
        valid_min, valid_max = -0.25, 2.25
        invalid = np.where(np.logical_or(l2i.freeboard < valid_min,
                                         l2i.freeboard > valid_max))[0]
        l2i.freeboard[invalid] = np.nan



        # step 2: recalculate sea ice thickness
        water_density = np.full(l2i.n_records, self.config.water_density)
        l2i.sea_ice_thickness = snowfreeboard2thickness(
                l2i.freeboard, l2i.snow_depth, water_density,
                l2i.ice_density, l2i.snow_density)

        valid_min, valid_max = -0.5, 10.5
        invalid = np.where(np.logical_or(l2i.sea_ice_thickness < valid_min,
                                         l2i.sea_ice_thickness > valid_max))[0]
        l2i.sea_ice_thickness[invalid] = np.nan

        output_filename = os.path.join(output_folder,
                                       filename_from_path(l2i_file))
        output = L2iNCFileWrapper()
        output.write_to_file(output_filename, l2i)


class L2iNCFileWrapper(NCDataFile):

    def __init__(self):
        super(L2iNCFileWrapper, self).__init__()
        self.parameter_list = [
                "timestamp", "longitude", "latitude",  "surface_type",
                "elevation", "mean_sea_surface", "sea_surface_anomaly",
                "radar_freeboard", "freeboard", "sea_ice_type",
                "sea_ice_concentration", "snow_depth", "snow_density",
                "ice_density", "sea_ice_thickness"]
        self.parameter = []
        self.l2 = None
        self._missing_parameters = []
        self.parameter_attributes = get_parameter_attributes("l2i")

    def set_filename(self, path):
        self.base_export_path = path

    def get_full_export_path(self, startdt):
        self._get_full_export_path(startdt)
        return self.export_path

    def write_to_file(self, filename, l2):
        dimdict = OrderedDict([("n_records", l2.n_records)])
        self.path = filename
        self._open_file()
        self._create_root_group(self._get_attr_dict(l2))
        self._populate_data_groups(l2, dimdict)
        self._write_to_file()

    def _get_attr_dict(self, l2):
        attr_dict = OrderedDict()
        for attribute_name in l2.attribute_list:
            attr_dict[attribute_name] = getattr(l2.info, attribute_name)
        return attr_dict

    def _populate_data_groups(self, l2, dimdict):
        dims = dimdict.keys()
        for key in dims:
                self._rootgrp.createDimension(key, dimdict[key])
        for parameter_name in self.parameter_list:
            data = getattr(l2, parameter_name)
            # Convert datetime objects to number
            if type(data[0]) is datetime:
                data = date2num(data, self.time_def.units,
                                self.time_def.calendar)
            # Convert bool objects to integer
            if data.dtype.str == "|b1":
                data = np.int8(data)
            dimensions = tuple(dims[0:len(data.shape)])
            var = self._rootgrp.createVariable(
                    parameter_name, data.dtype.str, dimensions, zlib=self.zlib)
            var[:] = data
            # Add Parameter Attributes
            attribute_dict = self._get_variable_attr_dict(parameter_name)
            for key in attribute_dict.keys():
                setattr(var, key, attribute_dict[key])

        # Report mission variable attributes (not in master release)
        not_master = "master" not in PYSIRAL_VERSION
        if not_master and len(self._missing_parameters) > 0:
            print "Warning: Missing parameter attributes for "+"; ".join(
                self._missing_parameters)


class CRDPConverterConfig(object):

    def __init__(self):
        self.l2i_source_repo = None
        self.water_density = 1024.0

    def parse_args(self):
        args = self.argparser.parse_args()
        self.l2i_source_repo = args.l2i_repo

    def get_output_folder(self, year, month):
        output_folder = self.l2i_source_repo.replace("p002a", "p002b")
        output_folder = os.path.join(output_folder, year, month)
        if output_folder == self.l2i_source_repo:
            sys.exit("something went wrong, files may be overwritten")
        return output_folder

    @property
    def argparser(self):
        parser = argparse.ArgumentParser()
        parser.add_argument(
                action='store',
                dest='l2i_repo',
                type=str,
                help='l2i repository (..\l2i)')
        return parser


class CRDPConverterL2iList(DefaultLoggingClass):

    def __init__(self):
        super(CRDPConverterL2iList, self).__init__(self.__class__.__name__)
        self.error = ErrorStatus()
        self.l2i_search = r"l2i*.nc"
        self.stack = dict()
        self.n_files = 0

    def get_file_stack(self, repo):
        self.validate_repo_path(repo)
        if self.error.status:
            return
        msg = "List input files (%s) in %s" % (self.l2i_search, repo)
        self.log.info(msg)
        self.get_directory_tree(repo)
        self.list_l2i_files(repo)
        msg = "Found %g l2i files" % self.n_files
        self.log.info(msg)

    def validate_repo_path(self, repo):
        """ Check if valid path """
        is_p002a = re.match(".*p002a.*", repo)
        path_exists = os.path.isdir(repo)
        ends_on_l2i = os.path.split(repo)[1] == "l2i"
        if not path_exists or not ends_on_l2i or not is_p002a:
            msg = "Invalid argument (path to l2i folder (p002a) required)"
            self.error.add_error("invalid-l2i-path", msg)

    def get_directory_tree(self, repo):
        """ Make a dict tree of the directory structure """
        for root, dirs, files in walklevel(repo, level=1):
            # skip top level folder
            if root == repo:
                continue
            # this can be done as walk level is limited to 1
            year = os.path.split(root)[1]
            if not re.match(r'([1-3][0-9]{3})', year):
                continue
            # weed out any subdirectories that are not 01, 02, ... , 12
            month_list = [d for d in dirs if re.match(r'([0-1][0-9])', d)]
            self.add_l2i_branch(year, month_list)

    def add_l2i_branch(self, year, month_list):
        if year not in self.stack.keys():
            self.stack[year] = dict()
            for month in month_list:
                self.stack[year][month] = []
        else:
            self.log.warning("Overwriting directory tree %s" % year)

    def list_l2i_files(self, repo):
        self.n_files = 0
        for year in self.years:
            for month in self.get_month_list(year):
                search = os.path.join(repo, year, month, self.l2i_search)
                l2i_list = glob.glob(search)
                self.stack[year][month] = l2i_list
                self.n_files += len(l2i_list)

    def get_month_list(self, year):
        return sorted(self.stack[year].keys())

    def get_l2i_list(year, month):
        pass

    @property
    def years(self):
        return sorted([year for year in self.stack.keys()])

    @property
    def iteration(self):
        output = []
        for year in self.years:
            for month in self.get_month_list(year):
                iteritem = (year, month, self.stack[year][month])
                output.append(iteritem)
        return output


def walklevel(some_dir, level=1):
    some_dir = some_dir.rstrip(os.path.sep)
    assert os.path.isdir(some_dir)
    num_sep = some_dir.count(os.path.sep)
    for root, dirs, files in os.walk(some_dir):
        yield root, dirs, files
        num_sep_this = root.count(os.path.sep)
        if num_sep + level <= num_sep_this:
            del dirs[:]

if __name__ == "__main__":
    job = CRDPConvertP002aToP002b()
    job.config.parse_args()
    job.collect_l2i_files()
    if job.l2irepo.error.status:
        job.l2irepo.error.raise_on_error()
    job.convert_files()
