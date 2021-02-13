# -*- coding: utf-8 -*-


from pysiral import psrlcfg
from pysiral.config import get_yaml_config
from pysiral.errorhandler import ErrorStatus
from pysiral.logging import DefaultLoggingClass
from loguru import logger
from netCDF4 import Dataset, date2num
from datetime import datetime
from dateutil import parser as dtparser
from collections import OrderedDict
from pathlib import Path
from attrdict import AttrDict
import numpy as np
import cftime
import parse
import re


class OutputHandlerBase(DefaultLoggingClass):
    """
    A class that defines properties of output files (content, location, format)
    based on the output definition, data container and other processor settings
    """

    subfolder_format = {"month": "%02g", "year": "%04g", "day": "%02g"}

    def __init__(self, output_def, applicable_data_level=None, subfolder_tags=None,
                 default_file_location=None):
        """
        Init the output handler with the content of the output definition file and
        keywords specific for the data processing levels. These keywords have to be
        set during the initialization of the parent (this) class.

        TODO: Move applicable data level and subfolder tags to output definition file
        TODO: Remove default_file_location

        :param output_def: (str or pathlib.Path): The full file path to the output definition file
        :param applicable_data_level: (int) Intended processing level for validation of output definition
        :param subfolder_tags: (str list) A list of intended sub-folders and their meaning
            (example ["year", "month"])
        :param default_file_location: A list relative to the pysiral resource directory of sub-directories
            and filenames that links to the default output definition for the respectice
            data processing level (deprecated)
        """

        # Init the parent
        super(OutputHandlerBase, self).__init__(self.__class__.__name__)
        self.pysiral_config = psrlcfg
        self.error = ErrorStatus()
        self._basedir = "n/a"

        # Attributes
        self._doi = None
        self.subfolders = None
        self.subfolder_tags = subfolder_tags
        self.applicable_data_level = applicable_data_level
        self.default_file_location = default_file_location

        self._init_from_output_def(output_def)
        self.output_def_filename = output_def

    def fill_template_string(self, template, dataset):
        """ Fill an template string with information of a dataset
        object (in this case Level2Data) """
        attributes = self.get_template_attrs(template)
        result = str(template)
        for attribute in attributes:
            attribute_name, option, placeholder = attribute
            attribute = dataset.get_attribute(attribute_name, *option)
            if attribute is None:
                attribute = "unknown"
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

    @staticmethod
    def get_template_attrs(template):
        """ Extract attribute names and options (if defined) for a
        give template string """
        try:
            template = template.encode('utf-8').strip()
        except AttributeError:
            template = str(template)
        attr_defs = re.findall("{.*?}", str(template))
        attrs, options = [], []
        for attr_def in attr_defs:
            attr_name, _, optstr = attr_def[1:-1].partition(":")
            attrs.append(attr_name)
            options.append(optstr.split(";"))
        return zip(attrs, options, attr_defs)

    def _init_from_output_def(self, output_def):
        """ Adds the information for the output def yaml files (either
        full filename or treedict structure) """
        if Path(output_def).is_file():
            try:
                self._output_def = get_yaml_config(output_def)
            except Exception as ex:
                self.error.add_error("outputdef-parser-error", ex)
                self.error.raise_on_error()
        else:
            self._output_def = output_def
        self._validate_outputdef()

    def _set_basedir(self, basedir, create=True):
        """ Sets and and (per default) creates the main output directory """
        self._basedir = basedir
        if create:
            self._create_directory(self._basedir)

    def _create_directory(self, directory):
        """ Convinience method to create a directory and add an error
        when failed """
        Path(directory).mkdir(exist_ok=True, parents=True)
        if not Path(directory).is_dir():
            msg = "Unable to create directory: %s" % str(directory)
            self.error.add_error("directory-error", msg)

    def _get_subdirectories(self, dt):
        directory = Path(self.basedir)
        for subfolder_tag in self.subfolders:
            parameter = getattr(dt, subfolder_tag)
            subfolder = self.subfolder_format[subfolder_tag] % parameter
            directory = directory / subfolder

    def _get_directory_from_dt(self, dt):
        subfolders = self.get_dt_subfolders(dt, self.subfolder_tags)
        return Path(self.basedir) / Path(*subfolders)

    def _validate_outputdef(self):
        """ Run a series of tests to check if a valid output definition
        has been passed. Note: theses tests will only check existing
        items of the output definition. If the requested item is missing
        a separate exception will be evoked """
        # Test 1: Applicable data level needs
        if self.applicable_data_level != self.data_level:
            msg = "outputdef data level (%g) does not match %s reqirement (%g)"
            msg = msg % (self.data_level, self.__class__.__name__,
                         self.applicable_data_level)
            self.error.add_error("datalevel-mismatch", msg)
            self.error.raise_on_error()

    @property
    def has_doi(self):
        try:
            return self._doi is not None
        except AttributeError:
            return False

    @property
    def id(self):
        try:
            return self._output_def.metadata.output_id
        except (AttributeError, KeyError):
            return None

    @property
    def product_level_subfolder(self):
        subfolder = self._output_def.product_level_subfolder
        if type(subfolder) is not str:
            msg = "root.product_level_subfolder (str) missing or wrong dtype"
            self.error.add_error("outputdef-invalid", msg)
            self.error.raise_on_error()
        return subfolder

    @property
    def data_level(self):
        data_level = self._output_def.metadata.data_level
        if type(data_level) is not int:
            msg = "root.metadata.data_level (int) missing or wrong dtype"
            self.error.add_error("outputdef-invalid", msg)
            self.error.raise_on_error()
        return data_level

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
        variables = sorted(list(self.output_def.variables.keys()))
        attribute_dicts = [self.output_def.variables[a] for a in variables]
        return zip(variables, attribute_dicts)


class DefaultLevel2OutputHandler(OutputHandlerBase):
    """
    Default output handler with pysiral conventions. Uses product directory from
    local_machine_def.yaml as standard repository
    """

    def __init__(self, output_def="default", subdirectory="default_output",
                 overwrite_protection=True, period="default"):
        """
        Init the Level-2 output handler
        :param output_def: (str) The full ilepath of the output definition file
        :param subdirectory: (str) The subdirectory relativ to the standard pysiral output path defined
            in the local_machine_def.yaml file
        :param overwrite_protection: (bool) If true an additional sub-directory will be added to the
            output path with the exact date and time of this run. The intention is to avoid over-writing
            existing files
        :param period: (str) An identifier string for the period. This string will be used to choose
            the filename template from the output definition file
        """

        # Fall back to default output if no output_def is given
        # (allows default initialization for the Level2 processor)
        if output_def == "default":
            output_def = self.default_output_def_filename

        # Init the parent class with Level-2 data output parameters
        super(DefaultLevel2OutputHandler, self).__init__(
            output_def, applicable_data_level=2, subfolder_tags=["year", "month"],
            default_file_location=["settings", "outputdef", "l2i_default.yaml"])

        self.error.caller_id = self.__class__.__name__
        logger.name = self.__class__.__name__

        self.subdirectory = subdirectory
        self.overwrite_protection = overwrite_protection
        self._period = period

        # Set parameters for this class
        self._init_product_directory()

    def get_filename_from_data(self, l2):
        """
        Return the filename for a defined level-2 data object based on tag filenaming
        in output definition file
        :param l2: l2data.Level2Data
        :return: str
        """

        # Get the filename definition (depending on period definition)
        filename_template = ""
        try:
            template_ids = self.output_def.filenaming.keys()
            period_id = self._period
            # Fall back to default if no filename convention for given
            # data period
            if period_id not in template_ids:
                period_id = "default"
            filename_template = self.output_def.filenaming[period_id]
        except AttributeError:
            filename_template = self.output_def.filenaming
        except KeyError:
            msg = "Missing filenaming convention for period [%s] in [%s]"
            msg = msg % (str(self._period), self.output_def_filename)
            self.error.add_error("invalid-outputdef", msg)
            self.error.raise_on_error()

        # Fill the tags in the filename template
        filename = self.fill_template_string(filename_template, l2)
        return filename

    def get_directory_from_data(self, l2, create=True):
        """
        Return the output directory based on information provided in an l2 data object
        :param l2:
        :param create:
        :return:
        """
        directory = self._get_directory_from_dt(l2.info.start_time)
        if create:
            self._create_directory(directory)
        return directory

    def get_fullpath_from_data(self, l2):
        """
        Return export path and filename based on information provided in the l2 data object
        :param l2:
        :return:
        """
        return Path(self.get_directory_from_data(l2)) / self.get_filename_from_data(l2)

    def get_global_attribute_dict(self, l2):
        attr_dict = OrderedDict()
        for attr_entry in self.output_def.global_attributes:
            attr_name, attr_template = zip(*attr_entry.items())
            attribute = self.fill_template_string(attr_template[0], l2)
            attr_dict[attr_name[0]] = attribute
        return attr_dict

    def remove_old(self, time_range):
        """ This method will erase all files in the target orbit for a
        given time range. Use with care """

        # Get the target directory
        # XXX: Assumption time_range is monthly
        directory = Path(self._get_directory_from_dt(time_range.start))
        # Get list of output files
        l2output_files = list(directory.glob("*.*"))

        # Delete files
        logger.info("Removing %g l2 product files [ %s ] in %s" % (len(l2output_files), self.id, directory))
        for l2output_file in l2output_files:
            Path(l2output_file).unlink()

    def _init_product_directory(self):
        """ Get main product directory from local_machine_def, add mandatory
        runtag subdirectory, optional second subdirectory for overwrite
        protection and product level id subfolder"""
        basedir = Path(psrlcfg.local_machine.product_repository)
        if not isinstance(self.subdirectory, list):
            basedir = basedir / self.subdirectory
        else:
            basedir = basedir / Path(*self.subdirectory)
        if self.overwrite_protection:
            basedir = basedir / self.now_directory
        basedir = basedir / self.product_level_subfolder
        self._set_basedir(basedir)

    @property
    def default_output_def_filename(self):
        local_settings_path = psrlcfg.pysiral_local_path
        return Path(local_settings_path) / Path(*self.default_file_location)


class NCDateNumDef(object):
    """
    Holds definition for datetime conversion to numbers and vice versa
    for netCDF operations
    """

    def __init__(self):
        self.units = "seconds since 1970-01-01"
        self.calendar = "standard"


class NCDataFile(DefaultLoggingClass):

    def __init__(self, output_handler):
        """
        Init the netCDF output parent class.
        NOTE: This class should only be used as a parent class.
        :param output_handler: An output handler class for the different processing level
        """

        # Init parent
        class_name = self.__class__.__name__
        super(NCDataFile, self).__init__(class_name)
        self.error = ErrorStatus(caller_id=class_name)

        # Output handler property
        self.output_handler = output_handler

        # Class attributes
        self.data = None
        self.filename = None
        self.base_export_path = None
        self.parameter_attributes = None

        self.time_def = NCDateNumDef()

        # TODO: Make this an option?
        self.zlib = True

        self._rootgrp = None
        self._options = None
        self._proc_settings = None
        self.verbose = False

    def set_options(self, **opt_dict):
        self._options = AttrDict(**opt_dict)

    def set_processor_settings(self, proc_settings):
        self._proc_settings = proc_settings

    def set_base_export_path(self, path):
        self.base_export_path = path

    def _set_doi(self):
        if self.output_handler.has_doi:
            self.data.set_doi(self.output_handler.doi)

    def _set_data_record_type(self):
        if self.output_handler.has_doi:
            self.data.set_data_record_type(self.output_handler.data_record_type)

    def _write_global_attributes(self):
        attr_dict = self.output_handler.get_global_attribute_dict(self.data)
        self._set_global_attributes(attr_dict)

    def _populate_data_groups(self, level3=False, flip_yc=False):

        lonlat_parameter_names = ["lon", "lat", "longitude", "latitude"]

        dimdict = self.data.dimdict
        dims = dimdict.keys()

        for key in dims:
            self._rootgrp.createDimension(key, dimdict[key])

        for parameter_name, attribute_dict in self.output_handler.variable_def:
            # Check if parameter name is also the the name or the source
            # parameter

            if "var_source_name" in attribute_dict.keys():
                attribute_dict = dict(attribute_dict)
                var_source_name = attribute_dict.pop("var_source_name")
            else:
                var_source_name = parameter_name

            data = self.data.get_parameter_by_name(var_source_name, raise_on_error=False)

            if data is None:
                msg = "Invalid parameter name for data object: %s"
                msg = msg % parameter_name
                logger.error(msg)
                self.error.add_error("invalid-paramater", msg)
                self.error.raise_on_error()

            # Convert datetime objects to number
            if isinstance(data[0], (datetime, cftime.datetime, cftime.real_datetime)):
                data = date2num(data, self.time_def.units, self.time_def.calendar)

            # Convert bool objects to integer
            if data.dtype.str == "|b1":
                data = np.int8(data)

            # Set dimensions (dependent on product level)
            if level3:
                if flip_yc:
                    data = np.flipud(data)
                if parameter_name not in lonlat_parameter_names:
                    data = np.array([data])
                    dimensions = tuple(list(dims)[0:len(data.shape)])
                else:
                    dimensions = tuple(list(dims)[1:len(data.shape)+1])
            else:
                dimensions = tuple(list(dims)[0:len(data.shape)])

            # Create and set the variable
            var = self._rootgrp.createVariable(parameter_name, data.dtype.str, dimensions, zlib=self.zlib)
            var[:] = data

            # Add Parameter Attributes
            for key in sorted(attribute_dict.keys()):
                attribute = attribute_dict[key]
                attribute = self.output_handler.fill_template_string(attribute, self.data)
                setattr(var, key, attribute)

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
            if isinstance(content, (datetime, cftime.datetime, cftime.real_datetime)):
                attdict[key] = date2num(content, self.time_def.units, self.time_def.calendar)

    @staticmethod
    def _convert_bool_attributes(attdict):
        """
        Replace l1b info parameters of type bool ['b1'] by a integer
        representation to match requirements for netCDF attribute data type
        rules
        """
        for key in attdict.keys():
            content = attdict[key]
            if type(content) is bool:
                attdict[key] = int(content)

    @staticmethod
    def _convert_nonetype_attributes(attdict):
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
        for key in attdict.keys():
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
        try:
            self._rootgrp = Dataset(self.full_path, "w")
        except RuntimeError:
            msg = "Unable to create netCDF file: %s" % self.full_path
            self.error.add_error("nc-runtime-error", msg)
            self.error.raise_on_error()

    def _write_to_file(self):
        self._rootgrp.close()

    @property
    def export_path(self):
        """ Evoking this property will also create the directory if it
        does not already exists """
        return self.output_handler.get_directory_from_data(self.data, create=True)

    @property
    def export_filename(self):
        """ Returns the filename for the level2 output file """
        return self.output_handler.get_filename_from_data(self.data)

    @property
    def full_path(self):
        return Path(self.export_path) / self.export_filename


class L1bDataNC(DefaultLoggingClass):
    """
    Class to export a L1bdata object into a netcdf file
    NOTE: This class is different from the other netCDF output classes as
          it uses data groups.
    """

    def __init__(self):

        super(L1bDataNC, self).__init__(self.__class__.__name__)

        self.datagroups = ["waveform", "surface_type", "time_orbit", "classifier", "correction"]
        self.filename = None
        self.time_def = NCDateNumDef()
        self.zlib = True
        self._rootgrp = None
        self._options = None
        self._proc_settings = None
        self.verbose = False

        self.output_folder = None
        self.l1b = None
        # TODO: Remove parameter attributes altogether with l1 setting files
        self.parameter_attributes = []

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
        self.path = Path(self.output_folder) / self.filename

    def _create_filename(self):
        self.filename = Path(self.l1b.filename.stem+".nc")

    def _set_global_attributes(self, attdict, prefix=""):
        """ Save l1b.info dictionary as global attributes """
        for key in sorted(attdict.keys()):
            self._rootgrp.setncattr(prefix+key, attdict[key])

    def _create_root_group(self, attdict, **global_attr_keyw):
        """
        Create the root group and add l1b metadata as global attributes
        :param attdict: A dictionary containing the global attributes
        :param global_attr_keyw:
        :return:
        """
        self._convert_datetime_attributes(attdict)
        self._convert_bool_attributes(attdict)
        self._convert_nonetype_attributes(attdict)
        self._set_global_attributes(attdict, **global_attr_keyw)

    def _populate_data_groups(self):

        self._missing_parameters = []

        for datagroup in self.datagroups:

            if self.verbose:
                print(datagroup.upper())

            # Create the data group
            dgroup = self._rootgrp.createGroup(datagroup)
            content = getattr(self.l1b, datagroup)

            # Create the dimensions
            # (must be available as OrderedDict in data group container
            dims = content.dimdict.keys()
            for key in dims:
                dgroup.createDimension(key, content.dimdict[key])

            # Now add variables for each parameter in data group
            for parameter in content.parameter_list:

                data = getattr(content, parameter)

                # Convert datetime objects to number
                if isinstance(data[0], (datetime, cftime.datetime, cftime.real_datetime)):
                    data = date2num(data, self.time_def.units, self.time_def.calendar)

                # Convert bool objects to integer
                if data.dtype.str == "|b1":
                    data = np.int8(data)
                dimensions = tuple(list(dims)[0:len(data.shape)])
                if self.verbose:
                    print(" "+parameter, dimensions, data.dtype.str, data.shape)

                var = dgroup.createVariable(parameter, data.dtype.str, dimensions, zlib=self.zlib)
                var[:] = data

                # Add Parameter Attributes
                attribute_dict = self._get_variable_attr_dict(parameter)
                for key in attribute_dict.keys():
                    setattr(var, key, attribute_dict[key])

    def _convert_datetime_attributes(self, attdict):
        """
        Change any datetime-like object to a netCDF compatible numerical representation
        :param attdict: dictionary-like
        :return:
        """

        for key in attdict.keys():
            content = attdict[key]
            if isinstance(content, (datetime, cftime.datetime, cftime.real_datetime)):
                attdict[key] = date2num(content, self.time_def.units, self.time_def.calendar)

    @staticmethod
    def _convert_bool_attributes(attdict):
        """
        Replace l1b info parameters of type bool ['b1'] by a integer representation
        to match requirements for netCDF attribute data type rules
        :param attdict: dictionary-like
        :return:
        """
        for key in attdict.keys():
            content = attdict[key]
            if type(content) is bool:
                attdict[key] = int(content)

    @staticmethod
    def _convert_nonetype_attributes(attdict):
        """
        Change empty (None) values in an attrdict to an empty string
        :param attdict: dictionary like
        :return: None, will be change in place
        """
        for key in attdict.keys():
            content = attdict[key]
            if content is None:
                attdict[key] = ""

    def _get_variable_attr_dict(self, parameter):
        """ Retrieve the parameter attributes """
        default_attrs = {"long_name": parameter}
        if parameter not in self.parameter_attributes:
            # self._missing_parameters.append(parameter)
            return default_attrs
        else:
            return dict(self.parameter_attributes[parameter])

    def _open_file(self):
        self._rootgrp = Dataset(self.path, "w")

    def _write_to_file(self):
        self._rootgrp.close()


class Level1POutput(NCDataFile):
    """ Class to export a l2data object into a netcdf file """

    def __init__(self, data, output_handler):
        super(Level1POutput, self).__init__(output_handler)
        self.data = data


class Level2Output(NCDataFile):
    """
    Class to export a l2data object into a netcdf file
    """

    def __init__(self, data, output_handler):
        """
        A class for data of processing level-2 (l2i, l2p)
        :param data: The data container (l2data.Level2Data)
        :param output_handler: The output handler (pysiral.output.DefaultLevel2OutputHandler)
        """

        # Init the parent
        super(Level2Output, self).__init__(output_handler)

        # Store the data container
        # FIXME: The data container does not need to be stored here, can be piped to _export_content()
        self.data = data
        self._set_doi()
        try:
            self._set_data_record_type()
        except AttributeError:
            pass
        self._export_content()

    def _export_content(self):
        """
        :return:
        """
        self.path = self.full_path
        self._open_file()
        self._write_global_attributes()
        self._populate_data_groups()
        self._write_to_file()


class Level3Output(NCDataFile):
    """
    Class to export a Level-3 data object into a netcdf file.
    Differences to Level2Output are small but substantial (e.g.
    with the additional time dimension)
    """

    def __init__(self, data, output_handler):
        super(Level3Output, self).__init__(output_handler)
        self.data = data
        self._set_doi()
        self._set_data_record_type()
        self._export_content()

    def _export_content(self):
        self._open_file()
        self._write_global_attributes()
        self._populate_data_groups(level3=True, flip_yc=self.output_handler.flip_yc)
        self._add_time_variables()
        self._add_grid_mapping_variables()
        self._write_to_file()

    def _add_time_variables(self):
        """
        Add the required variables describing the temporal coverage and the reference time of
        the Level-3 data product.
        :return:
        """

        rgrp = self._rootgrp

        # Set Time Variable
        var = rgrp.createVariable("time", "f8", tuple(['time']), zlib=self.zlib)
        var.standard_name = "time"
        var.units = self.time_def.units
        var.long_name = "reference time of product"
        var.axis = "T"
        var.calendar = self.time_def.calendar
        var.bounds = "time_bnds"
        var[:] = date2num(self.data.metadata.time_coverage_start,
                          self.time_def.units, self.time_def.calendar)

        # Set Time Bounds
        rgrp.createDimension("nv", 2)
        time_bounds_dt = self.data.time_bounds
        td_units, td_cal = self.time_def.units, self.time_def.calendar
        time_bnds = [[date2num(dt, td_units, td_cal) for dt in time_bounds_dt]]
        dims = ("time", "nv")
        var = rgrp.createVariable("time_bnds", "f8", dims, zlib=self.zlib)
        var.units = self.time_def.units
        var[:] = time_bnds

    def _add_grid_mapping_variables(self):
        """
        This method adds a variable to the netCDF file that contains the grid definition and is named
        after the projection. E.g:

        byte Lambert_Azimuthal_Grid;
            :latitude_of_projection_origin = 90.0; // double
            :longitude_of_projection_origin = 0.0; // double
            :proj4_string = "+proj=laea +lon_0=0 +datum=WGS84 +ellps=WGS84 +lat_0=90.0";
            :semi_major_axis = 6378137.0; // double
            :false_easting = 0.0; // double
            :false_northing = 0.0; // double
            :grid_mapping_name = "lambert_azimuthal_equal_area";
            :inverse_flattening = 298.257223563; // double

        The information is taken 1:1 from the `netcdf_grid_description` entry in the grid definition file.

        :return:
        """

        rgrp = self._rootgrp

        # Set x coordinate
        var = rgrp.createVariable("xc", "f8", tuple(['xc']), zlib=self.zlib)
        var.standard_name = "projection_x_coordinate"
        var.units = "km"
        var.long_name = "x coordinate of projection (eastings)"
        var[:] = self.data.griddef.xc_km

        # Set y coordinate
        yc_km = self.data.griddef.yc_km
        if self.output_handler.flip_yc:
            yc_km = np.flip(yc_km, 0)
        var = rgrp.createVariable("yc", "f8", tuple(['yc']), zlib=self.zlib)
        var.standard_name = "projection_y_coordinate"
        var.units = "km"
        var.long_name = "y coordinate of projection (eastings)"
        var[:] = yc_km

        # Set grid definition
        grid_nc_cfg = self.data.griddef.netcdf_vardef
        name = list(grid_nc_cfg.keys())[0]
        attrs = grid_nc_cfg[name]
        var = rgrp.createVariable(name, "i1", ())
        for key in sorted(attrs.keys()):
            setattr(var, key, attrs[key])


class PysiralOutputFilenaming(object):
    """
    Class for generating and parsing of pysiral output  filenames for all data levels
    TODO: Find out if this is still needed (probably for l1p) and remove if obsolet
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
            "l1p": "pysiral-l1p-{mission_id}-{source}-{timeliness}-{hemisphere}-{start}-{stop}-{version}.nc",
            "l2i": "l2i_{version}_{mission_id}_{hemisphere}_{start}_{stop}.nc",
            "l3s": "l3s_{version}_{mission_id}_{grid}_{resolution}_{start}_{stop}.nc"}

    def parse_filename(self, fn):
        """ Parse info from pysiral output filename """
        filename = Path(fn).name
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
                        except dtparser.ParserError:
                            match_found = False
                            break
                    setattr(self, parameter, value)
                break

        if not match_found:
            print("Unrecognized filename: %s" % filename)

    @staticmethod
    def _datetime_format(dt):
        return "{dt:%Y%m%dT%H%M%S}".format(dt=dt)


class PysiralOutputFolder(object):
    """
    Class for generating and retrieving output folders
    """

    def __init__(self):
        self.error = ErrorStatus()
        self.data_level = None
        self.path = None
        self.version = "default"
        self.mission_id = None
        self.year = None
        self.month = None
        self.hemisphere = None
        self.config = psrlcfg

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
            self.path = Path(self.path) / subfolder

    def create(self):
        Path(self.path).mkdir(exist_ok=True, parents=True)

    def _set_folder_as_l1bdata(self):
        self.data_level = "l1b"
        local_repository = self.config.local_machine.l1b_repository
        export_folder = local_repository[self.mission_id][self.version].l1bdata
        yyyy = "{:04g}".format(self.year)
        mm = "{:02g}".format(self.month)
        self.path = Path(export_folder) / self.hemisphere / yyyy / mm


def get_output_class(name):
    return globals()[name]()
