# -*- coding: utf-8 -*-

""" """

__all__ = ["auxdata", "bnfunc", "cryosat2", "envisat", "ers", "esa", "icesat", "sentinel3", "classifier", "clocks",
           "config", "datahandler", "errorhandler", "filter", "flag", "frb", "grid", "io_adapter",
           "iotools", "l1bdata", "l1bpreproc", "l2data", "l2preproc", "l2proc", "l3proc", "legacy",
           "logging", "maptools", "mask", "orbit", "output", "path", "proj", "retracker", "roi",
           "sit", "surface_type", "units", "validator", "waveform"]


import warnings
warnings.filterwarnings("ignore")

import os
import sys
from distutils import log, dir_util
log.set_verbosity(log.INFO)
log.set_threshold(log.INFO)
import importlib

# Get version from VERSION in package root
PACKAGE_ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
try:
    version_file = open(os.path.abspath(os.path.join(PACKAGE_ROOT_DIR, "VERSION")))
    with version_file as f:
        version = f.read().strip()
except IOError:
    sys.exit("Cannot find VERSION file in package (expected: %s" % version_file)

# Package Metadata
__version__ = version
__author__ = "Stefan Hendricks"
__author_email__ = "stefan.hendricks@awi.de"

# Get the config directory of the package
# NOTE: This approach should work for a local script location of distributed package
PACKAGE_CONFIG_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "pysiral-cfg")


# Get an indication of the location for the pysiral configuration path
# NOTE: In its default version, the text file `PYSIRAL-CFG-LOC` does only contain the
#       string `USER_HOME`. In this case, pysiral will expect the a .pysiral-cfg subfolder
#       in the user home. The only other valid option is an absolute path to a specific
#       directory with the same content as .pysiral-cfg. This was introduced to enable
#       fully encapsulated pysiral installation in virtual environments


# Get the home directory of the current user
CURRENT_USER_HOME_DIR = os.path.expanduser("~")


# Read pysiral config location indicator file
try:
    cfg_loc_file = open(os.path.abspath(os.path.join(PACKAGE_ROOT_DIR, "PYSIRAL-CFG-LOC")))
    with cfg_loc_file as f:
        cfg_loc = f.read().strip()
except IOError:
    sys.exit("Cannot find PYSIRAL-CFG-LOC file in package (expected: %s)" % cfg_loc_file)


# Case 1 (default): pysiral config path is in user home
if cfg_loc == "USER_HOME":

    # NOTE: This is where to expect the pysiral configuration files
    USER_CONFIG_PATH = os.path.join(CURRENT_USER_HOME_DIR, ".pysiral-cfg")

    # Check if pysiral configuration exists in user home directory
    # if not: create the user configuration directory
    # Also add an initial version of local_machine_def, so that pysiral does not raise an exception
    if not os.path.isdir(USER_CONFIG_PATH):
        print "Creating pysiral config directory: %s" % USER_CONFIG_PATH
        dir_util.copy_tree(PACKAGE_CONFIG_PATH, USER_CONFIG_PATH, verbose=1)
        print "Init local machine def"
        template_filename = os.path.join(PACKAGE_CONFIG_PATH, "templates", "local_machine_def.yaml.template")
        target_filename = os.path.join(USER_CONFIG_PATH, "local_machine_def.yaml")
        dir_util.copy(template_filename, USER_CONFIG_PATH, verbose=1)

# Case 2: package specific config path
else:
    USER_CONFIG_PATH = cfg_loc
    # This must be an existing path
    # NOTE: The idea is that this case needs to set up specifically (i.e. calling a script that will alter the value
    #       in the config location indicator file
    if not os.path.isdir(USER_CONFIG_PATH):
        sys.exit("Invalid pysiral config location in PYSIRAL-CFG-LOC (%s)" % USER_CONFIG_PATH)


def get_cls(module_name, class_name):
    """ Small helper function to dynamically load classes"""
    module = importlib.import_module(module_name)
    return getattr(module, class_name, None)