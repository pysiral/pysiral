# -*- coding: utf-8 -*-

"""

"""

import shutil
import argparse
from enum import StrEnum
from pathlib import Path
from typing import List, Optional, Union
from loguru import logger

from pysiral import psrlcfg
from pysiral.scripts._argparse_types import config_target_type



__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"


class ConfigActions(StrEnum):
    INIT = "init"
    SET = "set"


def config(
    action: ConfigActions = ConfigActions.INIT,
    target_or_directory: Optional[Union[Path, str]] = None,
    yes: bool = False,
) -> None:
    """
    Main entry point for the config script. Can either initialize the pysiral configuration
    in the default configuration directory or set/update the configuration to a specific
    directory.

    :param action:
    :param target_or_directory:
    :param yes:

    :return: None
    """
    match action:
        case ConfigActions.INIT:
            config_init()
        case ConfigActions.SET:
            config_set(target_or_directory, yes)
        case _:
            raise ValueError(f"Unknown action: {action}")


def config_init() -> None:
    """
    Initial pysiral configuration in the default configuration directory.

    :return: None
    """

    logger.info("Initializing pysiral configuration in the default configuration directory.")

    # Configuration Directory
    if psrlcfg.config_path.is_dir():
        qualifier = "exists - no action taken"
    else:
        qualifier = "created"
    logger.info(f"Configuration directory: {psrlcfg.config_path} [{qualifier}]")

    # Local Machine Definition
    if psrlcfg.local_machine_def_filepath.is_file():
        qualifier = "exists - no action taken"
    else:
        local_machine_def_copy_template(psrlcfg.local_machine_def_filepath)
        qualifier = "Template copied"
    logger.info(f"Local machine def path: {psrlcfg.local_machine_def_filepath} [{qualifier}]")


def config_set(target: Union[Path, str], yes: bool) -> None:
    """

    :param target:
    :param yes:

    :return: None
    """

    logger.info(f"Update the pysiral configuration directory with target {target}")

    # Copy the configuration files (if needed)
    match target:

        # Set to user home directory (copy needed)
        case "USER_HOME":
            logger.info(f"Copy configuration files to user home directory. [{psrlcfg.userhome_config_path}]")
            copy_config_files(psrlcfg.userhome_config_path, yes)

        case "PACKAGE":
            logger.info(f"Package configuration files will be used. [{psrlcfg.package_config_path}]")

        case Path():
            logger.info(f"Copy configuration files to dedicated directory. [{target}]")
            copy_config_files(target, yes)

        case _:
            raise ValueError(f"Unknown target: {target}")

    # Update the pysiral configuration
    psrlcfg.set_config_target(target)

    # Copy the local machine definition template (if needed)

    if not psrlcfg.local_machine_def_filepath.is_file():
        action = "[missing] -> copy template"
        local_machine_def_copy_template(psrlcfg.local_machine_def_filepath)
    else:
        action = "[exists] -> no action taken"
    logger.info(f"Expect local_machine_def.yaml at {psrlcfg.local_machine_def_filepath} {action}")

    # Write the location of the pysiral configuration to file for the next start
    logger.info(f"Write new configuration target to file: ({target})")
    set_pysiral_cfg_loc(target)


def local_machine_def_copy_template(target_filepath: Path) -> None:
    """
    Copy the local machine definition template to the target filepath.

    :param target_filepath:

    :return: None
    """
    target_filepath.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(
        psrlcfg.package_config_path / "templates" / "local_machine_def.yaml",
        target_filepath
    )


def copy_config_files(target_directory: Path, yes: bool) -> None:
    """
    Copy the configuration files to the target directory.

    :param target_directory: The
    :param yes:

    :return:
    """

    # Check if the target directory exists and is not empty
    # If so, ask the user for confirmation to overwrite the files.
    # The `yes` flag can be used to skip the confirmation.
    target_directory_not_empty =  target_directory.is_dir() and any(target_directory.iterdir())
    if target_directory_not_empty and not yes:
        var = input(f"Overwrite configuration files to {target_directory}? [YES/NO]: ")
        if var != "YES":
            logger.info("Abort")
            raise KeyboardInterrupt("User aborted the operation.")

    # Create the target directory if it does not exist
    if not target_directory.is_file():
        try:
            target_directory.mkdir(parents=True, exist_ok=True)
        except IOError:
            msg = f"Cannot create directory: {target_directory}"
            raise IOError(msg)

    # Copy the configuration files from the pysiral package to the target directory
    try:
        shutil.copytree(psrlcfg.package_config_path, target_directory, dirs_exist_ok=True)
    except:
        msg = f"Could not create directory: {target_directory}. Copying configuration files failed."
        raise IOError(msg)


def set_pysiral_cfg_loc(target):
    """
    Write the location of the pysiral configuration for the current package
    NOTE: If you don't know what this means: Please Don't! An incorrect setting
    can break your pysiral installation!

    :param target: The target directory or identifier (e.g., "USER_HOME", "PACKAGE", or a specific path)

    :return: None
    """
    cfg_loc_file = psrlcfg.package_path / "PYSIRAL-CFG-LOC"
    with open(str(cfg_loc_file), 'w') as fh:
        fh.write(str(target))


class ConfigScriptArguments(object):

    def __init__(self) -> None:
        self.parser = self.get_argument_parser()

    def get(self, args_list: List[str] = None) -> "argparse.Namespace":
        args = self.parser.parse_args() if args_list is None else self.parser.parse_args(args_list)
        if args.action == ConfigActions.SET and args.target_or_directory is None:
            raise argparse.ArgumentError(f"The --dir argument is required for action '{args.action.value}'")
        return args

    @staticmethod
    def get_argument_parser() -> argparse.ArgumentParser:
        """
            Set up the command line argument parser for the Level-2 Processor.

            :return: The argument parser object.
            """

        # create the parser
        parser = argparse.ArgumentParser(
            prog="pysiral config",
            description="""
                        Initialize pysiral configuration and set or update the configuration directory.
                        """,
            epilog="For more information, see: https://pysiral.readthedocs.io",
            formatter_class=lambda prog: argparse.HelpFormatter(prog, width=96, indent_increment=4)  # noqa: E501
        )

        parser.add_argument(
            "action",
            type=str,
            choices=[v.value for v in ConfigActions],
            metavar="<action>",
            help=f"Configuration action to perform: {[v.value for v in ConfigActions]}"
        )

        parser.add_argument(
            "-t", "--target",
            type=config_target_type,
            dest="target_or_directory",
            metavar="<target_directory>",
            help="""
            Target directory for the configuration set or update actions.
            """
        )

        parser.add_argument(
            "--yes",
            action="store_true",
            dest="yes",
            default=False,
            required=False,
            help="Skips manual confirmation of actions that modify or overwrite existing configuration files."
        )

        return parser



