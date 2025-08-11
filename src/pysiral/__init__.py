# -*- coding: utf-8 -*-

"""
pysiral is the PYthon Sea Ice Radar ALtimetry toolbox
"""

__all__ = ["_logger", "auxdata", "cryosat2", "envisat", "ers", "sentinel3",
           "filter", "frb", "grid",
           "l1data", "l1preproc", "l2data", "l2preproc", "l2proc", "l3proc",
           "mask", "proj", "retracker",
           "sit", "surface", "waveform", "psrlcfg", "import_submodules", "get_cls",
           "set_psrl_cpu_count", "__version__", "__software_version__"]

import importlib

import multiprocessing
import pkgutil
import shutil
import socket
import sys
from datetime import datetime, UTC
from pathlib import Path
from typing import Iterable, Union

import yaml
from dateperiods import DatePeriod
from loguru import logger

from pysiral.core.legacy_classes import AttrDict


# Get version from VERSION in package root
PACKAGE_ROOT_DIR = Path(__file__).absolute().parent
VERSION_FILE_PATH = PACKAGE_ROOT_DIR / "VERSION"
try:
    version_file = open(str(VERSION_FILE_PATH))
    with version_file as f:
        version = f.read().strip()
except IOError:
    sys.exit(f'Cannot find VERSION file in package (expected: {PACKAGE_ROOT_DIR / "VERSION"}')


# Package Metadata
__version__ = version
__author__ = "Stefan Hendricks"
__author_email__ = "stefan.hendricks@awi.de"

# Get git version (allows tracing of the exact commit)
# TODO: Implement git version retrieval with git hooks (server-side post receive hook)
__software_version__ = None


class _MissionDefinitionCatalogue(object):
    """
    Container for storing and querying information from mission_def.yaml
    """

    def __init__(self, filepath):
        """
        Create a catalogue for all altimeter missions definitions
        :param filepath:
        """

        # Store Argument
        self._filepath = filepath

        # Read the file and store the content
        self._content = None
        with open(str(self._filepath)) as fh:
            self._content = AttrDict(yaml.safe_load(fh))

    def get_platform_info(self, platform_id) -> Union[AttrDict, None]:
        """
        Return the full configuration attr dict for a given platform id

        :param platform_id:

        :return:
        """
        platform_info = self._content.platforms.get(platform_id, None)
        return platform_info if platform_info is None else AttrDict(**platform_info)

    def get_platform_id(self, platform_name: str) -> Union[str, None]:
        """
        Return the name of a platform.
        :param platform_name:
        :return:
        """

        # Query the source dictionary
        platforms = [entry for entry in self._content.platforms.items() if entry[1]["long_name"] == platform_name]

        # No valid entry found -> Warning and returning None
        if not platforms:
            logger.warning(f"Did not find entry for {platform_name} in {self._filepath}")
            return None

        # Multiple Entries -> Error in configuration: Raise Exception
        elif len(platforms) > 1:
            msg = f"Multitple entries found for {platform_name} in {self._filepath}"
            logger.error(msg)
            raise ValueError(msg)

        platform_id, _ = platforms[0]
        return platform_id

    def get_name(self, platform_id):
        """
        Return the name of a platform.
        :param platform_id:
        :return:
        """
        platform_info = self.get_platform_info(platform_id)
        return None if platform_info is None else platform_info.long_name

    def get_sensor(self, platform_id):
        """
        Return the sensor name of a platform
        :param platform_id:
        :return:
        """
        platform_info = self.get_platform_info(platform_id)
        return None if platform_info is None else platform_info.sensor

    def get_orbit_inclination(self, platform_id):
        """
        Return the orbit inclination of a platform
        :param platform_id:
        :return:
        """
        platform_info = self.get_platform_info(platform_id)
        return None if platform_info is None else platform_info.orbit_max_latitude

    def get_time_coverage(self, platform_id):
        """
        Get the time coverage (start and end of data coverage) of the requested plaform.
        If the end data is not defined because the platform is still active, the current
        date is returned.
        :param platform_id:
        :return: time coverage start & time coverage end
        """
        platform_info = self.get_platform_info(platform_id)
        if platform_info is None:
            return None, None
        tcs = platform_info.time_coverage.start
        tce = platform_info.time_coverage.end
        if tce is None:
            tce = datetime.now(UTC)
        return tcs, tce

    @property
    def content(self) -> AttrDict:
        """
        The content of the definition file as an attribute-enabled dictionary.
        :return:
        """
        return self._content

    @property
    def ids(self):
        """
        A list of id's for each platform.

        :return: list with platform ids
        """
        return list(self.content.platforms.keys())


class _AuxdataCatalogueItem(object):
    """
    Container for an auxiliary data item
    """

    def __init__(self, category, auxid, config_dict):
        """
        Data class to manage an auxiliary data set definition
        :param category:
        :param auxid:
        :param config_dict:
        """

        # Arguments
        self._category = category
        self._id = auxid
        self._config_dict = config_dict

    @property
    def id(self):
        return str(self._id)

    @property
    def category(self):
        return str(self._category)

    @property
    def keys(self):
        return self._config_dict.keys()

    @property
    def attrdict(self):
        return AttrDict(**self._config_dict)


class _AuxdataCatalogue(object):
    """
    Container for the content of the auxdata_def.yaml definition file
    for auxiliary data sets
    """

    def __init__(self, filepath):
        """
        Data container with query functionality for auxiliary data.
        :param filepath:
        """

        # Arguments
        self.filepath = filepath

        # Read contents
        with open(str(self.filepath)) as fh:
            self._yaml_content = AttrDict(**yaml.safe_load(fh))

        # Create a catalogue of the content
        self.ctlg = {}
        for category, auxdata_items in self._yaml_content.items():
            self.ctlg[category] = {}
            for auxdata_id in auxdata_items:
                entry_dict = self._yaml_content[category][auxdata_id]
                item = _AuxdataCatalogueItem(category, auxdata_id, entry_dict)
                self.ctlg[category][auxdata_id] = item

    def get_category_items(self, category):
        """
        List all id's in a given category
        :param category:
        :return: list of ids
        """

        # Sanity check
        if category not in self.categories:
            raise ValueError(f'Invalid category: {str(category)} [{", ".join(self.categories)}]')

        # Return a sorted str list
        return sorted(self.ctlg[category].keys())

    def get_definition(self, category, auxid):
        """
        Retrieve the auxiliary data definition for a category and auxiliary data set id
        :param category: (str) Auxiliary data category (must be in self.categories)
        :param auxid: (str) The ID of the auxililary data set
        :return: AttrDict or None (if auxiliary dataset does not exist
        """

        # Check if valid category
        if category not in self.categories:
            return None

        # Extract & return the definition
        return self.ctlg[category].get(auxid, None)

    @property
    def categories(self):
        return self.ctlg.keys()

    @property
    def iter_keys(self):
        """
        List with two items per entry: (category, id)
        :return:
        """
        keys = []
        for category in self.categories:
            ids = self.get_category_items(category)
            keys.extend((category, auxid) for auxid in ids)
        return keys

    @property
    def items(self):
        """
        List with three items per entry: (category, id, catalogue_entry)
        :return:
        """
        keys = self.iter_keys
        return [(category, auxid, self.ctlg[category][auxid]) for category, auxid in keys]


class _PysiralPackageConfiguration(object):
    """
    Container for the content of the pysiral definition files
    (in pysiral/configuration) and the local machine definition file
    (local_machine_definition.yaml)
    """

    # --- Global variables of the pysiral package configuration ---

    # Filenames of definitions files
    _DEFINITION_FILES = {
        "platforms": "mission_def.yaml",
        "auxdata": "auxdata_def.yaml",
    }

    # name of the file containing the data path on the local machine
    _LOCAL_MACHINE_DEF_FILE = "local_machine_def.yaml"

    # valid settings types, processor levels and data level ids.
    # NOTE: These tie into the naming and content of definition files
    VALID_SETTING_TYPES = ["proc", "output", "grid"]
    VALID_PROCESSOR_LEVELS = ["l1", "l2", "l3"]
    VALID_DATA_LEVEL_IDS = ["l1", "l2", "l2i", "l2p", "l3", None]
    VALID_CONFIG_TARGETS = ["PACKAGE", "USER_HOME"]

    # Multiprocessing properties
    # Allow to package-wide specification of number of CPU's. Default value
    # is the CPU count from the python multiprocessing package.
    #
    # NOTE: This is intended when `multiprocessing.cpu_count()` is unreliable, e.g.
    #       when using slurm workload managers, or when the number of CPU's should be
    #       limited due to other concerns.
    CPU_COUNT = multiprocessing.cpu_count()

    def __init__(self):
        """
        Collect package configuration data from the various definition files and provide an interface
        to pysiral processor, output and grid definition files.
        This class is intended to be only called inside the init module of pysiral and to store the
        pysiral package configuration in the global variable `psrlcfg`
        """

        # --- Establish the path information ---
        # This step gets the default path (user home, set path for the resources)
        # NOTE: The current path to the active pysiral package is already set in the global
        #       variable `pysiral.PACKAGE_ROOT_DIR`, since this is required for reading the
        #       version file
        self._path = AttrDict()
        self._get_pysiral_path_information()

        # --- Check the dedicated pysiral config path ---
        # As per default, the pysiral config path is set in the user home directory.
        self._check_pysiral_config_path()

        # --- Read the configuration files ---
        self.local_machine = None
        self._read_config_files()

    def _get_pysiral_path_information(self):
        """
        Get the different path information for pysiral. This method will add the following
        attributes to self.path:
            1. package_root_path: The root directory of this package
            2. package_config_path: The directory of the pysiral config items in this package
            3. userhome_config_dir: The intended configuration directory in the user home
            4. config_target: The value given in the `PYSIRAL-CFG-LOC` file
        :return: None
        """

        # Store the root dir of this pysiral package
        self._path["package_root_path"] = PACKAGE_ROOT_DIR

        # Get the config directory of the package
        # NOTE: This approach should work for a local script location or an installed package
        self._path["package_config_path"] = self._path["package_root_path"] / "resources" / "pysiral-cfg"

        # Get an indication of the location for the pysiral configuration path
        # NOTE: In its default version, the text file `PYSIRAL-CFG-LOC` does only contain the
        #       string `USER_HOME`. In this case, pysiral will expect the a .pysiral-cfg sub-folder
        #       in the user home. The only other valid option is an absolute path to a specific
        #       directory with the same content as .pysiral-cfg. This was introduced to enable
        #       fully encapsulated pysiral installation in virtual environments

        # Get the home directory of the current user
        self._path["userhome_config_path"] = Path.home() / ".pysiral-cfg"

        # Read pysiral config location indicator file
        cfg_loc_file = PACKAGE_ROOT_DIR / "PYSIRAL-CFG-LOC"
        try:
            with open(str(cfg_loc_file)) as fh:
                self._path["config_target"] = fh.read().strip()
        except IOError:
            sys.exit(f"Cannot find PYSIRAL-CFG-LOC file in package (expected: {cfg_loc_file})")

    def _check_pysiral_config_path(self):
        """
        This class ensures that the pysiral configuration files are in the chosen
        configuration directory
        :return:
        """

        # Make alias of
        config_path = Path(self.config_path)
        package_config_path = Path(self.path.package_config_path)

        # Check if current config dir is package config dir
        # if yes -> nothing to do (files are either there or aren't)
        if config_path == package_config_path:
            return

        # current config dir is not package dir and does not exist
        # -> must be populated with content from the package config dir
        if not config_path.is_dir():
            print(f"Creating pysiral config directory: {config_path}")
            shutil.copytree(str(self.path.package_config_path), str(config_path), dirs_exist_ok=True)
            print("Init local machine def")
            template_filename = package_config_path / "templates" / "local_machine_def.yaml"
            target_filename = config_path / "local_machine_def.yaml"
            shutil.copy(str(template_filename), str(target_filename))

    def _read_config_files(self):
        """
        Read the three main configuration files for
            1. supported platforms
            2. supported auxiliary datasets
            3. path on local machine
        and create the necessary catalogues
        :return:
        """

        # --- Get information of supported platforms ---
        # The general information for supported radar altimeter missions (mission_def.yaml)
        # provides general metadata for each altimeter missions that can be used to sanity checks
        # and queries for sensor names etc.
        #
        # NOTE: This is just general information on altimeter platform and not to be confused with
        #       settings for actual primary data files. These are located in each l1p processor
        #       definition file.
        self.mission_def_filepath = self.config_path / Path(self._DEFINITION_FILES["platforms"])
        if not self.mission_def_filepath.is_file():
            error_msg = "Cannot load pysiral package files: \n %s" % self.mission_def_filepath
            print(error_msg)
            sys.exit(1)
        self.platforms = _MissionDefinitionCatalogue(self.mission_def_filepath)

        # --- Get information on supported auxiliary data sets ---
        # The auxdata_def.yaml config file contains the central definition of the properties
        # of supported auxiliary data sets. Each auxiliary data set is uniquely defined by
        # the type of auxiliary data set and a name id.
        # The central definition allows accessing auxiliary data by its id in processor definition files
        self.auxdata_def_filepath = self.config_path / self._DEFINITION_FILES["auxdata"]
        if not self.auxdata_def_filepath.is_file():
            error_msg = "Cannot load pysiral package files: \n %s" % self.auxdata_def_filepath
            print(error_msg)
            sys.exit(1)
        self.auxdef = _AuxdataCatalogue(self.auxdata_def_filepath)

        # read the local machine definition file
        self._read_local_machine_file()

    @staticmethod
    def get_yaml_config(filename):
        """
        Read a yaml file and return it content as an attribute-enabled dictionary
        :param filename: path to the yaml file
        :return: attrdict.AttrDict
        """
        with open(str(filename)) as fileobj:
            settings = AttrDict(yaml.safe_load(fileobj))
        return settings

    def get_setting_ids(self, settings_type, data_level=None):
        lookup_directory = self.get_local_setting_path(settings_type, data_level)
        ids, files = self.get_yaml_setting_filelist(lookup_directory)
        return ids

    def get_platform_period(self, platform_id):
        """
        Get a period definition for a given platform ID
        :param platform_id:
        :return: dateperiods.DatePeriod
        """
        tcs, tce = self.platforms.get_time_coverage(platform_id)
        return DatePeriod(tcs, tce)

    def get_processor_definition_ids(self, processor_level):
        """
        Returns a list of available processor definitions ids for a given processor
        level (see self.VALID_PROCESSOR_LEVELS)
        :param processor_level:
        :return:
        """
        lookup_directory = self.get_local_setting_path("proc", processor_level)
        return self.get_yaml_setting_filelist(lookup_directory, return_value="ids")

    def get_settings_files(self, settings_type: str, data_level: str) -> Iterable[Path]:
        """
        Returns all processor settings or output definitions files for a given data level.
        :param settings_type:
        :param data_level:
        :return:
        """

        if settings_type not in self.VALID_SETTING_TYPES:
            return []

        if data_level not in self.VALID_DATA_LEVEL_IDS:
            return []

        # Get all settings files in settings/{data_level} and its
        # subdirectories
        lookup_directory = self.get_local_setting_path(settings_type, data_level)
        _, files = self.get_yaml_setting_filelist(lookup_directory)

        # Test if ids are unique and return error for the moment
        return files

    def get_settings_file(self, settings_type, data_level, setting_id_or_filename):
        """ Returns a processor settings file for a given data level.
        (data level: l2 or l3). The second argument can either be a
        direct filename (which validity will be checked) or an id, for
        which the corresponding file (id.yaml) will be looked up in
        the default directory """

        if settings_type not in self.VALID_SETTING_TYPES:
            return None

        if data_level not in self.VALID_DATA_LEVEL_IDS:
            return None

        # Check if filename
        if Path(setting_id_or_filename).is_file():
            return setting_id_or_filename

        # Get all settings files in settings/{data_level} and its
        # subdirectories
        lookup_directory = self.get_local_setting_path(settings_type, data_level)
        ids, files = self.get_yaml_setting_filelist(lookup_directory)

        # Test if ids are unique and return error for the moment
        if len(set(ids)) != len(ids):
            msg = f"Non-unique {settings_type}-{str(data_level)} setting filename"
            print(f"ambiguous-setting-files: {msg}")
            sys.exit(1)

        # Find filename to setting_id
        try:
            index = ids.index(setting_id_or_filename)
            return Path(files[index])
        except (IOError, ValueError):
            return None

    @staticmethod
    def get_yaml_setting_filelist(directory, return_value="both"):
        """ Retrieve all yaml files from a given directory (including
        subdirectories). Directories named "obsolete" are ignored if
        ignore_obsolete=True (default) """
        setting_ids = []
        setting_files = []
        for filepath in directory.rglob("*.yaml"):
            setting_ids.append(filepath.name.replace(".yaml", ""))
            setting_files.append(filepath)
        if return_value == "both":
            return setting_ids, setting_files
        elif return_value == "ids":
            return setting_ids
        elif return_value == "files":
            return setting_files
        else:
            raise ValueError(f"Unknown return value {str(return_value)} [`both`, `ids`, `files`]")

    def get_local_setting_path(self, settings_type, data_level=None):
        """
        Return the absolute path on the local productions system to the configuration file. The
        returned path depends on the fixed structure below the `resources` directory in the pysiral
        package and the choice in the config file "PYSIRAL-CFG-LOC"
        :param settings_type:
        :param data_level:
        :return:
        """
        if settings_type in self.VALID_SETTING_TYPES and data_level in self.VALID_DATA_LEVEL_IDS:
            args = [settings_type]
            if data_level is not None:
                args.append(data_level)
            return Path(self.config_path) / Path(*args)
        else:
            return None

    def reload(self):
        """
        Method to trigger reading the configuration files again, e.g. after changing the config target
        :return:
        """
        self._read_config_files()
        self._check_pysiral_config_path()

    def set_config_target(self, config_target, permanent=False):
        """
        Set the configuration target
        :param config_target:
        :param permanent:
        :return:
        """

        # Input validation
        if config_target in self.VALID_CONFIG_TARGETS or Path(config_target).is_dir():
            self._path["config_target"] = config_target
        else:
            msg = "Invalid config_target: {} must be {} or valid path"
            msg = msg.format(str(config_target), ", ".join(self.VALID_CONFIG_TARGETS))
            raise ValueError(msg)

        if permanent:
            raise NotImplementedError()

    def _read_local_machine_file(self):
        """
        :return:
        """
        filename = self.local_machine_def_filepath
        try:
            local_machine_def = self.get_yaml_config(filename)
        except IOError:
            msg = f"local_machine_def.yaml not found (expected: {filename})"
            print(f"local-machine-def-missing: {msg}")
            local_machine_def = None
        self.local_machine = local_machine_def

    @property
    def platform_ids(self):
        return self.platforms.ids

    @property
    def path(self):
        return AttrDict(**self._path)

    @property
    def userhome_config_path(self):
        return Path(self.path.userhome_config_path)

    @property
    def package_config_path(self):
        return Path(self.path.package_config_path)

    @property
    def package_path(self):
        return Path(PACKAGE_ROOT_DIR)

    @property
    def current_config_target(self):
        return str(self._path["config_target"])

    @property
    def config_target(self):
        return str(self._path["config_target"])

    @property
    def config_path(self):
        """
        nstruct the target config path based on the value in `PYSIRAL-CFG-LOC`
        :return:
        """
        # Case 1 (default): pysiral config path is in user home
        if self._path["config_target"] == "USER_HOME":
            return Path(self._path["userhome_config_path"])

        # Case 2: pysiral config path is the package itself
        elif self._path["config_target"] == "PACKAGE":
            return Path(self._path["package_config_path"])

        # Case 3: package specific config path
        else:
            # This should be an existing path, but in the case it is not, it is created
            return Path(self._path["config_target"])

    @property
    def local_machine_def_filepath(self):
        if self.current_config_target != "PACKAGE":
            return self.config_path / self._LOCAL_MACHINE_DEF_FILE
        # TODO: Disable warnings that are run on simple import (as long as everything runs)
        # msg = "Current config path is `PACKAGE`, lookup directory for local_machine_def.yaml changed to `USERHOME`"
        # logger.warning(msg)
        return self.userhome_config_path / self._LOCAL_MACHINE_DEF_FILE

    @property
    def processor_levels(self):
        return list(self.VALID_PROCESSOR_LEVELS)

    @property
    def hostname(self):
        return socket.gethostname()

    @property
    def version(self):
        return str(__version__)


# Create a package configuration object as global variable
psrlcfg = _PysiralPackageConfiguration()


def get_cls(module_name, class_name, relaxed=True):
    """ Small helper function to dynamically load classes"""
    try:
        module = importlib.import_module(module_name)
    except ImportError as e:
        if relaxed:
            return None, e
        else:
            raise ImportError(f"Cannot load module: {module_name}") from e
    try:
        return getattr(module, class_name), None
    except AttributeError as e:
        if relaxed:
            return None, e
        else:
            raise NotImplementedError(f"Cannot load class: {module_name}.{class_name}") from e


def import_submodules(package, recursive=True):
    """ Import all submodules of a module, recursively, including subpackages

    :param package: package (name or actual module)
    :param recursive: Flag if package is a submodule
    :type package: str | module
    :rtype: dict[str, types.ModuleType]
    """
    if isinstance(package, str):
        package = importlib.import_module(package)
    results = {}
    for loader, name, is_pkg in pkgutil.walk_packages(package.__path__):
        full_name = f'{package.__name__}.{name}'
        results[full_name] = importlib.import_module(full_name)
        if recursive and is_pkg:
            results.update(import_submodules(full_name))
    return results


def set_psrl_cpu_count(cpu_count: int) -> None:
    """
    Set the pysiral-wide CPU count for multiprocessing to the pysiral package
    configuration

    :param cpu_count: The number of CPU's to use

    :raises ValueError: cpu_count is not a positive integer
    """

    try:
        assert isinstance(cpu_count, int)
        assert cpu_count > 0
    except AssertionError as e:
        raise ValueError(
            f"specified number of CPU's ({cpu_count}) not a positive integer"
        ) from e
    cpu_count_mp = multiprocessing.cpu_count()
    if cpu_count > cpu_count_mp:
        logger.warning(f"Specified number of CPU's ({cpu_count}) > number of CPU's ({cpu_count_mp})")
    psrlcfg.CPU_COUNT = cpu_count
