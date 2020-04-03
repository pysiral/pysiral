"""
A short script that should be used to set the location of the pysiral config file for the
"""

import sys
import argparse
import shutil
import os
from distutils import log, dir_util
log.set_verbosity(log.INFO)
log.set_threshold(log.INFO)

from pysiral import PACKAGE_CONFIG_PATH, USER_CONFIG_PATH, PACKAGE_ROOT_DIR

def main(args):
    """
    Handle the actions of
     1. create a new specific config directory
     2. activating an specific config directory
     4. reset default config directory in user home
    NOTE: The mutually exclusive selection of the 3 actions is ensured by argparse
    """

    print("Create, activate or reset the pysiral configuration directory (handle with care!)")

    # Action 1:
    if args.target_create is not None:

        var = input("Create and activate config dir %s? [YES/NO]: " % args.target_create)
        if var != "YES":
            print(" abort")
            sys.exit(2)

        # Check if directory exists
        if os.path.isdir(args.target_create):
            print("Error: Directory exists: %s" % str(args.target_create))
            print("Creating new config dir [FAILED]")
            sys.exit(1)

        # a. Create the user home directory
        try:
            dir_util.copy_tree(PACKAGE_CONFIG_PATH, args.target_create, verbose=1)
        except:
            print("Error: Could not create directory: %s" % str(args.target_create))
            print("Creating new config dir [FAILED]")
            sys.exit(1)

        # b. Copy the local_machine_def.yaml

        # b1: Specific files
        if args.lmd != "":
            if not os.path.isfile(args.lmd):
                print("Error: Could copy local_machine_def.yaml: %s" % str(args.lmd))
                print("Creating new config dir [FAILED]")
                sys.exit(1)

            head, tail = os.path.split(args.lmd)
            if tail != "local_machine_def.yaml":
                print("Error: argument not named `local_machine_def.yaml`: %s" % str(args.lmd))
                print("Creating new config dir [FAILED]")
                sys.exit(1)
                shutil.copy(args.lmd, os.path.join(args.target_create, "local_machine_def.yaml"))

        # b2: from template
        else:
            template_filename = os.path.join(PACKAGE_CONFIG_PATH, "templates", "local_machine_def.yaml.template")
            target_filename = os.path.join(args.target_create, "local_machine_def.yaml")
            shutil.copy(template_filename, target_filename)

        # c. Change the value in PYSIRAL-CFG-LOC
        set_pysiral_cfg_loc(args.target_create)
        print("Creating new config dir [SUCCESS]")
        sys.exit(0)

    # Action 3:
    if args.target_activate is not None:

        var = input("Activate config dir %s? [YES/NO]: " % args.target_activate)
        if var != "YES":
            print(" abort")
            sys.exit(2)

        if os.path.isdir(args.target_activate):
            set_pysiral_cfg_loc(args.target_activate)
            print("Activating config dir [SUCCESS]")
            sys.exit(0)
        else:
            print("Error: Invalid directory: %s" % str(args.target_activate))
            print("Activating config dir [FAILED]")
            sys.exit(1)

    # Action 3:
    if args.reset:

        var = input("Reset config dir to user home? [YES/NO]: ")
        if var != "YES":
            print(" abort")
            sys.exit(2)

        set_pysiral_cfg_loc("USER_HOME")
        print("Resetting config dir to user home [SUCCESS]")
        sys.exit(0)


def set_pysiral_cfg_loc(target):
    """
    Write the location of the pysiral configuration for the current package
    NOTE: If you don't know what this means: Don't!
    """
    cfg_loc_file = os.path.join(PACKAGE_ROOT_DIR, "PYSIRAL-CFG-LOC")
    with open(cfg_loc_file, 'w') as fh:
        fh.write(target)


def get_args():

    # create the parser
    parser = argparse.ArgumentParser()

    action_group = parser.add_mutually_exclusive_group(required=True)

    # Action 1: Create a new pysiral config directory
    action_group.add_argument('-create-config-dir', dest='target_create', action='store', default=None, type=str,
                               help='Target directory for specific config path')

    # Action 2: Activate an existing config directory
    action_group.add_argument('-activate-target-dir', dest='target_activate', action='store', default=None, type=str,
                               help='Target directory for specific config path')

    # Action 3: Reset to user home default
    action_group.add_argument('--reset', dest='reset', action='store_true', default=False,
                              help='Set this keyword to reset the config path to $user_home')

    # (Optional) copy a defined local machine def (will be ignored for Actions 2 & 3)
    parser.add_argument('-use-local-machine-def', dest='lmd', action='store', default="", type=str,
                        help='Copy the following local_machine_def.yaml into config dir')

    # Get the Arguments
    args = parser.parse_args()

    return args

if __name__ == "__main__":
    args = get_args()
    main(args)