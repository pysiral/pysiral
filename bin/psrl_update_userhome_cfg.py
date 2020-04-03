"""
A short script that updates the settings files in the user home
"""

import sys

from distutils import log, dir_util
log.set_verbosity(log.INFO)
log.set_threshold(log.INFO)

from pysiral import PACKAGE_CONFIG_PATH, USER_CONFIG_PATH


def main():
    """ Copy all files from the package resources to the user home"""

    # Get confirmation from user since all settings file will be overwritten
    var = raw_input("All config files (except `local_machine_def.yaml` will be overwritten, type YES to proceed: ")
    if var != "YES":
        print " abort"
        sys.exit()

    # Copy the entire tree
    print "copy pysiral config files from package to user home: "
    dir_util.copy_tree(PACKAGE_CONFIG_PATH, USER_CONFIG_PATH, verbose=1)

if __name__ == "__main__":
    main()