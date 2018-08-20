# -*- coding: utf-8 -*-

""" """

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

# First get the home directory of the current user
# NOTE: This is where to expect the pysiral configuration files
CURRENT_USER_HOME_DIR = os.path.expanduser("~")

# Get the config directory of the package
# NOTE: This approach should work for a local script location of distributed package
PACKAGE_CONFIG_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "pysiral-cfg")

# Check if pysiral configuration exists in user home directory
# if not: create the user configuration directory
USER_CONFIG_PATH = os.path.join(CURRENT_USER_HOME_DIR, ".pysiral-cfg")
if not os.path.isdir(USER_CONFIG_PATH):
    print "Creating pysiral config directory: %s" % USER_CONFIG_PATH
    shutil.copytree(PACKAGE_CONFIG_PATH, USER_CONFIG_PATH)


