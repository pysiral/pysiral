# -*- coding: utf-8 -*-

""" """

__version__ = "0.6.1"

__all__ = ["bnfunc", "cryosat2", "envisat", "ers", "esa", "icesat", "sentinel3", "auxdata", "classifier", "clocks",
           "config", "datahandler", "errorhandler", "filter", "flag", "frb", "grid", "io_adapter",
           "iotools", "l1bdata", "l1bpreproc", "l2data", "l2preproc", "l2proc", "l3proc", "legacy",
           "logging", "maptools", "mask", "mss", "orbit", "output", "path", "proj", "retracker", "roi",
           "sic", "sit", "sitype", "snow", "surface_type", "units", "validator", "waveform"]

import warnings
warnings.filterwarnings("ignore")

import os
import sys
import shutil
import pkg_resources

# First get the home directory of the current user
# NOTE: This is where to expect the pysiral configuration files
CURRENT_USER_HOME_DIR = os.path.expanduser("~")

# Get the config directory of the package
# NOTE: This approach should work for a local script location of distributed package

PACKAGE_CONFIG_PATH = pkg_resources.resource_filename("pysiral", "resources/pysiral-cfg")

# Check if pysiral configuration exists in user home directory
# if not: create the user configuration directory
USER_CONFIG_PATH = os.path.join(CURRENT_USER_HOME_DIR, ".pysiral-cfg")
if not os.path.isdir(USER_CONFIG_PATH):
    print "Creating pysiral config directory: %s" % USER_CONFIG_PATH
    shutil.copytree(PACKAGE_CONFIG_PATH, USER_CONFIG_PATH)


