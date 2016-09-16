# -*- coding: utf-8 -*-


from pysiral.config import PYSIRAL_VERSION_FILENAME, ConfigInfo
from pysiral.path import filename_from_path, file_basename
from pysiral.errorhandler import ErrorStatus
from pysiral.config import options_from_dictionary, get_parameter_attributes
from pysiral.path import validate_directory


from netCDF4 import Dataset, date2num
from datetime import datetime
from dateutil import parser as dtparser
import numpy as np
import parse
import os


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
        self.verbose = False

    def set_options(self, **opt_dict):
        self._options = options_from_dictionary(**opt_dict)

    def _create_root_group(self, attdict):
        """
        Create the root group and add l1b metadata as global attributes
        """
        self._convert_datetime_attributes(attdict)
        self._convert_bool_attributes(attdict)
        self._convert_nonetype_attributes(attdict)
        self._set_global_attributes(attdict)

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
            self._rootgrp.setncattr(key, attdict[key])

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
        print "Warning: Missing parameter attributes for "+"; ".join(
            self._missing_parameters)

    def _get_variable_attr_dict(self, parameter):
        """ Retrieve the parameter attributes """
        default_attrs = {
            "long_name": parameter,
            "standard_name": parameter,
            "scale_factor": 1.0,
            "add_offset": 0.0}
        if not self.parameter_attributes.has_key(parameter):
            self._missing_parameters.append(parameter)
            return default_attrs
        else:
            return dict(self.parameter_attributes[parameter])


class L2iDataNC(NCDataFile):
    """
    Class to export a l2data object into a netcdf file
    """

    def __init__(self):
        super(L2iDataNC, self).__init__()
        self.parameter = []
        self.base_export_path = None
        self.l2 = None

    def set_base_export_path(self, path):
        self.base_export_path = path

    def get_full_export_path(self, startdt):
        self._get_full_export_path(startdt)
        return self.export_path

    def write_to_file(self, l2):
        self._get_full_export_path(l2.info.start_time)
        self._get_export_filename(l2)
        self._open_file()
        self._create_root_group(l2.info.attdict)
        self._populate_data_groups(l2)
        self._write_to_file()

    def _get_full_export_path(self, startdt):
        # Comput und create the export directory
        base_path = self.base_export_path
        sub_folders = self._options.subfolders
        folder = PysiralOutputFolder(load_config=False)
        folder.l2i_from_startdt(startdt, base_path, sub_folders)
        folder.create()
        self.export_path = folder.path

    def _get_export_filename(self, l2):
        # get full output filename
        filenaming = PysiralOutputFilenaming()
        self.filename = filenaming.from_l2i(l2)
        self.path = os.path.join(self.export_path, self.filename)

    def _populate_data_groups(self, l2):
        dimdict = l2.dimdict
        dims = dimdict.keys()
        for key in dims:
                self._rootgrp.createDimension(key, dimdict[key])
        for parameter_name in self._options.parameter:
            data = l2.get_parameter_by_name(parameter_name)
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

    def set_export_folder(self, path):
        self.export_path = path

    def set_metadata(self, metadata):
        self.metadata = metadata

    def export(self, l3):
        self._validate()
        self._open_file()
        self._create_root_group(self.metadata.attdict)
        self._populate_data_groups(l3)
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
            dimensions = tuple(dims[0:len(data.shape)])
            var = self._rootgrp.createVariable(
                    parameter_name, data.dtype.str, dimensions, zlib=self.zlib)
            var[:] = data


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
            mission=l3s.mission,
            gri=l3s.grid_tag,
            resolution_tag=l3s.resolution_tag,
            start=self._datetime_format(l3s.start_period),
            stop_period=self._datetime_format(l3s.stop_period))
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
                        value = dtparser.parse(value)
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



#def get_l1bdata_export_folder(l1b, config=None, version="default"):
#    """ Returns the l1bdata export folder for a l1b data object """
#    if config is None or not isinstance(config, ConfigInfo):
#            config = ConfigInfo()
#
#    local_repository = config.local_machine.l1b_repository
#    export_folder = local_repository[mission_id][version].l1bdata
#    yyyy = "%04g" % year
#    mm = "%02g" % month
#    export_folder = os.path.join(export_folder, hemisphere, yyyy, mm)
#    return export_folder

def get_output_class(name):
    return globals()[name]()
