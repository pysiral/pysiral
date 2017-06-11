# -*- coding: utf-8 -*-


from pysiral.config import (PYSIRAL_VERSION, PYSIRAL_VERSION_FILENAME,
                            ConfigInfo, get_yaml_config)
from pysiral.path import filename_from_path, file_basename
from pysiral.errorhandler import ErrorStatus
from pysiral.logging import DefaultLoggingClass
from pysiral.config import options_from_dictionary, get_parameter_attributes
from pysiral.path import validate_directory

from glob import glob
from netCDF4 import Dataset, date2num
from datetime import datetime
from dateutil import parser as dtparser
import numpy as np
import parse
import os
import re


class OutputHandlerBase(DefaultLoggingClass):

    subfolder_format = {"month": "%02g", "year": "%04g", "day": "%02g"}

    def __init__(self, output_def):
        super(OutputHandlerBase, self).__init__(self.__class__.__name__)
        self.pysiral_config = ConfigInfo()
        self.error = ErrorStatus()
        self._basedir = "n/a"
        self._init_from_output_def(output_def)

    def fill_template_string(self, template, dataset):
        """ Fill an template string with information of a dataset
        object (in this case Level2Data) """
        attributes = self.get_template_attrs(template)
        result = str(template)
        for attribute in attributes:
            attribute_name, option, placeholder = attribute
            attribute = dataset.get_attribute(attribute_name, *option)
            result = result.replace(placeholder, attribute)
        return result

    def get_dt_subfolders(self, dt, subfolder_tags):
        """ Returns a list of subdirectories based on a datetime object
        (usually the start time of data collection) """
        subfolders = []
        for subfolder_tag in subfolder_tags:
            parameter = getattr(dt, subfolder_tag)
            subfolder = self.subfolder_format[subfolder_tag] % parameter
            subfolders.append(subfolder)
        return subfolders

    def get_template_attrs(self, template):
        """ Extract attribute names and options (if defined) for a
        give template string """
        attr_defs = re.findall("{.*?}", str(template))
        attrs, options = [], []
        for attr_def in attr_defs:
            attr_name, _, optstr = attr_def[1:-1].partition(":")
            attrs.append(attr_name)
            options.append(optstr.split(","))
        return zip(attrs, options, attr_defs)

    def _init_from_output_def(self, output_def):
        """ Adds the information for the output def yaml files (either
        full filename or treedict structure) """
        if os.path.isfile(output_def):
            try:
                self._output_def = get_yaml_config(output_def)
            except Exception, msg:
                self.error.add_error("outputdef-parsers-error", msg)
                self.error.raise_on_error()
        else:
            self._output_def = output_def

    def _set_basedir(self, basedir, create=True):
        """ Sets and and (per default) creates the main output directory """
        self._basedir = basedir
        if create:
            self._create_directory(self._basedir)

    def _create_directory(self, directory):
        """ Convinience method to create a directory and add an error
        when failed """
        status = validate_directory(directory)
        if not status:
            msg = "Unable to create directory: %s" % str(directory)
            self.error.add_error("directory-error", msg)

    def _get_subdirectories(self, dt):
        directory = self.basedir
        for subfolder_tag in self.subfolders:
            parameter = getattr(dt, subfolder_tag)
            subfolder = self.subfolder_format[subfolder_tag] % parameter
            directory = os.path.join(directory, subfolder)

    def _get_directory_from_dt(self, dt):
        subfolders = self.get_dt_subfolders(dt, self.subfolder_tags)
        return os.path.join(self.basedir, *subfolders)

    @property
    def id(self):
        try:
            return self._output_def.metadata.output_id
        except:
            return None

    @property
    def product_level_subfolder(self):
        subfolder = self._output_def.product_level_subfolder
        if type(subfolder) is not str:
            msg = "root.product_level_subfolder missing"
            self.error.add_error("outputdef-invalid", msg)
            self.error.raise_on_error()
        return subfolder

    @property
    def basedir(self):
        return self._basedir

    @property
    def output_def(self):
        return self._output_def

    @property
    def now_directory(self):
        """ Returns a directory suitable string with the current time """
        return datetime.now().strftime("%Y%m%dT%H%M%S")

    @property
    def variable_def(self):
        t = self.output_def.variables
        variables = list(t.iterkeys(recursive=False, branch_mode='only'))
        attribute_dicts = [self.output_def.variables[a] for a in variables]
        return zip(variables, attribute_dicts)


class DefaultLevel2OutputHandler(OutputHandlerBase):
    """ Default output handler with pysiral conventions. Uses product
    directory from local_machine_def.yaml as standard repository """

    # Some fixed parameters for this class
    default_file_location = ["settings", "outputdef", "l2i_default.yaml"]
    subfolder_tags = ["year", "month"]
    applicable_data_level = 2

    def __init__(self, output_def="default", subdirectory="default_output",
                 overwrite_protection=True):
        # Fall back to default output if no output_def is given
        # (allows default initialization for the Level2 processor)
        if output_def == "default":
            output_def = self.default_output_def_filename
        super(DefaultLevel2OutputHandler, self).__init__(output_def)
        self.error.caller_id = self.__class__.__name__
        self.log.name = self.__class__.__name__
        self.subdirectory = subdirectory
        self.overwrite_protection = overwrite_protection
        self._init_product_directory()

    def get_filename_from_l2(self, l2):
        """ Return the filename for a defined level-2 data object
        based on tag filenaming in output definition file """
        filename_template = self.output_def.filenaming
        return self.fill_template_string(filename_template, l2)

    def get_directory_from_l2(self, l2, create=True):
        """ Return the output directory based on information provided
        in an l2 data object """
        directory = self._get_directory_from_dt(l2.info.start_time)
        if create:
            self._create_directory(directory)
        return directory

    def get_fullpath_from_l2(self, l2):
        """ Return export path and filename based on information
        provided in the l2 data object """
        export_directory = self.get_directory_from_l2(l2)
        export_filename = self.get_filename_from_l2(l2)
        return os.path.join(export_directory, export_filename)

    def get_global_attribute_dict(self, l2):
        attr_dict = {}
        for attr_name in self.output_def.global_attributes.iterkeys():
            attr_template = self.output_def.global_attributes[attr_name]
            attribute = self.fill_template_string(attr_template, l2)
            attr_dict[attr_name] = attribute
        return attr_dict

    def remove_old(self, time_range):
        """ This method will erase all files in the target orbit for a
        given time range. Use with care """

        # Get the target directory
        # XXX: Assumption time_range is monthly
        directory = self._get_directory_from_dt(time_range.start)
        # Get list of output files
        search_pattern = os.path.join(directory, "*.*")
        l2output_files = glob(search_pattern)

        # Delete files
        self.log.info("Removing %g l2 product files [ %s ] in %s" % (
                len(l2output_files), self.id, directory))
        for l2output_file in l2output_files:
                os.remove(l2output_file)

    def _init_product_directory(self):
        """ Get main product directory from local_machine_def, add mandatory
        runtag subdirectory, optional second subdirectory for overwrite
        protection and product level id subfolder"""
        pysiral_config = ConfigInfo()
        basedir = pysiral_config.local_machine.product_repository
        basedir = os.path.join(basedir, self.subdirectory)
        if self.overwrite_protection:
            basedir = os.path.join(basedir, self.now_directory)
        basedir = os.path.join(basedir, self.product_level_subfolder)
        self._set_basedir(basedir)



    @property
    def default_output_def_filename(self):
        pysiral_config = ConfigInfo()
        local_settings_path = pysiral_config.pysiral_local_path
        return os.path.join(local_settings_path, *self.default_file_location)


class NCDateNumDef(object):
    """
    Holds definition for datetime conversion to numbers and vice versa
    for netCDF operations
    """

    def __init__(self):
        self.units = "seconds since 1970-01-01"
        self.calendar = "standard"


class NCDataFile(object):

    def __init__(self):
        self.filename = None
        self.time_def = NCDateNumDef()
        self.zlib = True
        self._rootgrp = None
        self._options = None
        self._proc_settings = None
        self.verbose = False

    def set_options(self, **opt_dict):
        self._options = options_from_dictionary(**opt_dict)

    def set_processor_settings(self, proc_settings):
        self._proc_settings = proc_settings

    def _create_root_group(self, attdict, **global_attr_keyw):
        """
        Create the root group and add l1b metadata as global attributes
        """
        self._convert_datetime_attributes(attdict)
        self._convert_bool_attributes(attdict)
        self._convert_nonetype_attributes(attdict)
        self._set_global_attributes(attdict, **global_attr_keyw)

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

    def _set_global_attributes(self, attdict, prefix=""):
        """ Save l1b.info dictionary as global attributes """
        for key in sorted(attdict.keys()):
            self._rootgrp.setncattr(prefix+key, attdict[key])

    def _get_variable_attr_dict(self, parameter):
        """ Retrieve the parameter attributes """
        default_attrs = {
            "long_name": parameter,
            "standard_name": parameter,
            "scale_factor": 1.0,
            "add_offset": 0.0}
        if parameter not in self.parameter_attributes:
            # self._missing_parameters.append(parameter)
            return default_attrs
        else:
            return dict(self.parameter_attributes[parameter])

    def _write_processor_settings(self):
        if self._proc_settings is None:
            pass
        settings = self._proc_settings
        for item in settings.iterkeys():
            self._rootgrp.setncattr(item, str(settings[item]))

    def _open_file(self):
        self._rootgrp = Dataset(self.path, "w")

    def _write_to_file(self):
        self._rootgrp.close()


class L1bDataNC(NCDataFile):
    """
    Class to export a L1bdata object into a netcdf file
    """

    def __init__(self):
        super(L1bDataNC, self).__init__()

        self.datagroups = ["waveform", "surface_type", "time_orbit",
                           "classifier", "correction"]
        self.output_folder = None
        self.l1b = None
        self.parameter_attributes = get_parameter_attributes("l1b")

    def export(self):
        self._validate()
        self._open_file()
        # Save the l1b info data group as global attributes
        attdict = self.l1b.info.attdict
        self._create_root_group(attdict)
        self._populate_data_groups()
        self._write_to_file()

    def _validate(self):
        if self.filename is None:
            self._create_filename()
        self.path = os.path.join(self.output_folder, self.filename)

    def _create_filename(self):
        self.filename = file_basename(self.l1b.filename)+".nc"

    def _populate_data_groups(self):
        self._missing_parameters = []
        for datagroup in self.datagroups:
            if self.verbose:
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
                if self.verbose:
                    print " "+parameter, dimensions, data.dtype.str, data.shape
                var = dgroup.createVariable(
                    parameter, data.dtype.str, dimensions, zlib=self.zlib)
                var[:] = data
                # Add Parameter Attributes
                attribute_dict = self._get_variable_attr_dict(parameter)
                for key in attribute_dict.keys():
                    setattr(var, key, attribute_dict[key])

        # Report mission variable attributes (not in master release)
        not_master = "master" not in PYSIRAL_VERSION
        if not_master:
            print "Warning: Missing parameter attributes for "+"; ".join(
                self._missing_parameters)


class Level2Output(NCDataFile):
    """
    Class to export a l2data object into a netcdf file
    """

    def __init__(self, l2, output_handler):
        super(Level2Output, self).__init__()
        self.l2 = l2
        self.output_handler = output_handler
        self._export_l2()

    def set_base_export_path(self, path):
        self.base_export_path = path

    def get_full_export_path(self, startdt):
        self._get_full_export_path(startdt)
        return self.export_path

    def _export_l2(self):
        self.path = self.full_path
        self._open_file()
        self._write_global_attributes()
        self._populate_data_groups()
        self._write_to_file()

    def _write_global_attributes(self):
        attr_dict = self.output_handler.get_global_attribute_dict(self.l2)
        self._set_global_attributes(attr_dict)

    def _populate_data_groups(self):

        dimdict = self.l2.dimdict
        dims = dimdict.keys()

        for key in dims:
                self._rootgrp.createDimension(key, dimdict[key])

        for parameter_name, attribute_dict in self.output_handler.variable_def:

            data = self.l2.get_parameter_by_name(parameter_name)

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
            for key in attribute_dict.keys():
                setattr(var, key, attribute_dict[key])

    @property
    def export_path(self):
        """ Evoking this property will also create the directory if it
        does not already exists """
        return self.output_handler.get_directory_from_l2(self.l2, create=True)

    @property
    def export_filename(self):
        """ Returns the filename for the level2 output file """
        return self.output_handler.get_filename_from_l2(self.l2)

    @property
    def full_path(self):
        return os.path.join(self.export_path, self.export_filename)


class L3SDataNC(NCDataFile):
    """
    Class to export a l2data object into a netcdf file
    """

    def __init__(self):
        super(L3SDataNC, self).__init__()
        self.parameter = []
        self.export_path = None
        self.metadata = None
        self.l2 = None
        self.parameter_attributes = get_parameter_attributes("l3c")

    def set_export_folder(self, path):
        self.export_path = path

    def set_metadata(self, metadata):
        self.metadata = metadata

    def export(self, l3):
        self._validate()
        self._open_file()
        self._create_root_group(self.metadata.attdict)
        self._populate_data_groups(l3)
        self._add_time_dummy_variable()
        self._write_to_file()

    def export_parameter_dict(self, pardict, dimdict=None):
        self._validate()
        self._open_file()
        self._create_root_group(self.metadata.attdict)
        self._add_time_dummy_variable()
        self._populate_data_groups_from_dict(pardict, dimdict)
        self._write_to_file()

    def _validate(self):
        # Validate the export directory
        path = self.export_path
        validate_directory(path)
        # get full output filename
        filenaming = PysiralOutputFilenaming()
        filename = filenaming.from_l3s(self.metadata)
        self.path = os.path.join(path, filename)

    def _populate_data_groups(self, l3):
        dimdict = l3.dimdict
        dims = dimdict.keys()
        for key in dims:
                self._rootgrp.createDimension(key, dimdict[key])
        for parameter_name in l3.parameter_list:
            data = l3.get_parameter_by_name(parameter_name)
            # Add time axis
            if parameter_name not in ["lon", "lat"]:
                data = np.array([data])
                dimensions = tuple(dims[0:len(data.shape)])
            else:
                dimensions = tuple(dims[1:len(data.shape)+1])
            var = self._rootgrp.createVariable(
                    parameter_name, data.dtype.str, dimensions, zlib=self.zlib)
            var[:] = data
            # Add Parameter Attributes
            attribute_dict = self._get_variable_attr_dict(parameter_name)
            for key in attribute_dict.keys():
                setattr(var, key, attribute_dict[key])

    def _populate_data_groups_from_dict(self, pardict, dimdict):
        dims = dimdict.keys()
        for key in dims:
                self._rootgrp.createDimension(key, dimdict[key])
        for parameter_name in sorted(pardict.keys()):
            data = pardict[parameter_name]
            dimensions = tuple(dims[0:len(data.shape)])
            var = self._rootgrp.createVariable(
                    parameter_name, data.dtype.str, dimensions, zlib=self.zlib)
            var[:] = np.array([data])
            # Add Parameter Attributes
            attribute_dict = self._get_variable_attr_dict(parameter_name)
            for key in attribute_dict.keys():
                setattr(var, key, attribute_dict[key])

    def _add_time_dummy_variable(self):
        var = self._rootgrp.createVariable("time", "f8", ('time'),
                                           zlib=self.zlib)
        var.standard_name = "time"
        var.long_name = self.time_def.units
        var.units = "seconds"


class PysiralOutputFilenaming(object):
    """
    Class for generating and parsing of pysiral output
    filenames for all data levels
    """

    def __init__(self):
        self.error = ErrorStatus()
        self.data_level = None
        self.version = None
        self.hemisphere = None
        self.mission_id = None
        self.orbit = None
        self.start = None
        self.stop = None
        self.resolution = None
        self.grid = None

        self._registered_parsers = {
            "l1bdata": "l1bdata_{version}_{mission_id}_{hemisphere}_{start}_{stop}.nc",
            "l2i": "l2i_{version}_{mission_id}_{hemisphere}_{start}_{stop}.nc",
            "l3s": "l3s_{version}_{mission_id}_{grid}_{resolution}_{start}_{stop}.nc"}

    def from_l1b(self, l1b):
        """ Level-1b preprocessed filename """
        export_filename = self._registered_parsers["l1bdata"]
        export_filename = export_filename.format(
            version=PYSIRAL_VERSION_FILENAME,
            hemisphere=l1b.info.hemisphere,
            mission_id=l1b.mission,
            start=self._datetime_format(l1b.info.start_time),
            stop=self._datetime_format(l1b.info.stop_time))
        return export_filename

    def from_l2i(self, l2i):
        """ Level-2 Intermediate filename """
        export_filename = self._registered_parsers["l2i"]
        export_filename = export_filename.format(
            version=PYSIRAL_VERSION_FILENAME,
            hemisphere=l2i.hemisphere,
            mission_id=l2i.info.mission,
            start=self._datetime_format(l2i.info.start_time),
            stop=self._datetime_format(l2i.info.stop_time))
        return export_filename

    def from_l3s(self, l3s):
        """ Level-3 super-collocated filename """
        export_filename = self._registered_parsers["l3s"]
        export_filename = export_filename.format(
            version=PYSIRAL_VERSION_FILENAME,
            mission_id=l3s.mission,
            grid=l3s.grid_tag,
            resolution=l3s.resolution_tag,
            start=self._datetime_format(l3s.start_period),
            stop=self._datetime_format(l3s.stop_period))
        return export_filename

    def parse_filename(self, fn):
        """ Parse info from pysiral output filename """
        filename = filename_from_path(fn)
        match_found = False
        for data_level in self._registered_parsers.keys():
            parser = parse.compile(self._registered_parsers[data_level])
            match = parser.parse(filename)
            if match:
                match_found = True
                self.data_level = data_level
                for parameter in match.named.keys():
                    value = match[parameter]
                    if parameter in ["start", "stop"]:
                        try:
                            value = dtparser.parse(value)
                        except:
                            match_found = False
                            break
                    setattr(self, parameter, value)
                break

        if not match_found:
            print "Unrecognized filename: %s" % filename

    def _datetime_format(self, datetime):
        return "{dt:%Y%m%dT%H%M%S}".format(dt=datetime)


class PysiralOutputFolder(object):
    """
    Class for generating and retrieving output folders
    """

    def __init__(self, config=None, load_config=True):
        self.error = ErrorStatus()
        self.data_level = None
        self.path = None
        self.version = "default"
        self.mission_id = None
        self.year = None
        self.month = None
        if not load_config:
            return
        if config is None or not isinstance(config, ConfigInfo):
            self.config = ConfigInfo()
        else:
            self.config = config

    def l1bdata_from_list(self, mission_id, version, hemisphere, year, month):
        self.mission_id = mission_id
        self.version = version
        self.hemisphere = hemisphere
        self.year = year
        self.month = month
        self._set_folder_as_l1bdata()

    def l1bdata_from_l1b(self, l1b, version="default"):
        self.mission_id = l1b.mission
        self.version = version
        self.hemisphere = l1b.info.hemisphere
        self.year = l1b.info.start_time.year
        self.month = l1b.info.start_time.month
        self._set_folder_as_l1bdata()

    def l2i_from_startdt(self, startdt, base_path, subfolders):
        stringify = {"month": "%02g", "year": "%04g", "day": "%02g"}
        self.path = base_path
        for subfolder_tag in subfolders:
            parameter = getattr(startdt, subfolder_tag)
            subfolder = stringify[subfolder_tag] % parameter
            self.path = os.path.join(self.path, subfolder)

    def create(self):
        validate_directory(self.path)

    def _set_folder_as_l1bdata(self):
        self.data_level = "l1b"
        local_repository = self.config.local_machine.l1b_repository
        export_folder = local_repository[self.mission_id][self.version].l1bdata
        yyyy = "%04g" % self.year
        mm = "%02g" % self.month
        self.path = os.path.join(export_folder, self.hemisphere, yyyy, mm)


def get_output_class(name):
    return globals()[name]()
