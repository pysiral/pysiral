#!/usr/bin/env python
"""
A short script that updates the settings files in the user home
"""

import shutil
import sys

from loguru import logger

from pysiral import psrlcfg


def main():
    """ Copy all files from the package resources to the user home"""

    # Get confirmation from user since all settings file will be overwritten
    var = input("All config files (except `local_machine_def.yaml` will be overwritten, type YES to proceed: ")
    if var != "YES":
        print(" abort")
        sys.exit()

    # Check if pysiral config path is package
    # -> in this case copying the files is pointless
    if psrlcfg.config_target == "PACKAGE":
        logger.warning("Config target is `PACKAGE`, aborting ...")
        sys.exit()

    # Copy the entire tree
    logger.info("copy pysiral config files from package to config dir: {}".format(psrlcfg.config_path))
    shutil.copytree(str(psrlcfg.package_config_path), str(psrlcfg.config_path), dirs_exist_ok=True)


if __name__ == "__main__":
    main()
