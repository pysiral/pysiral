# -*- coding: utf-8 -*-

"""
"""

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"

import re
import argparse
from datetime import date, datetime
from pathlib import Path
from typing import Callable, Union, List, Literal

from pysiral import psrlcfg
from pysiral.core.flags import PysiralProcessingLevels


def dir_type(ends_with: Union[str, List[str]] = None, must_exist: bool = True) -> Callable:
    """
    Factory function that provides additional option to the argparse directory type validation.

    :param ends_with: Last subdirectory name that the directory path must end with.
    :param must_exist: If True, the directory must exist; if False, it will be created if it does not exist.

    :return: Function to validiates input for argparse
    """

    def dir_type_func(value: str) -> Path:

        if ends_with is not None:
            # Optional check: Directory path needs to end with certain directory name
            # NOTE: You cannot re-assign a variable from the outer scope in a nested function,
            # so we need to check the type of `ends_with` and create a list of valid last directories.
            valid_last_dirs = [ends_with] if isinstance(ends_with, str) else ends_with
            if not any(Path(value).name.endswith(suffix) for suffix in valid_last_dirs):
                raise argparse.ArgumentTypeError(f"{value} must end with: {', '.join(valid_last_dirs)}")

        # First Check: Directory must exist
        if Path(value).is_dir():
            value = Path(value).absolute()
        else:
            if must_exist:
                raise argparse.ArgumentTypeError(f"Not a directory: {value}")
            else:
                try:
                    value = Path(value).absolute()
                    # If the directory does not exist, we create it
                    value.mkdir(parents=True, exist_ok=True)
                except IOError:
                    raise argparse.ArgumentTypeError(f"Cannot create directory: {value}")

        return value

    return dir_type_func


def date_type(value: str) -> date:
    """
    Small helper function to convert input automatically to Path.
    Also raises an exception if the input is not a valid directory.

    :param value: Input argument

    :raises ValueError:

    :return: Input as `datetime.date`
    """
    try:
        return date.fromisoformat(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Not a date definition: {value}")


def date_datetime_type(value: str) -> datetime:
    try:
        return datetime.fromisoformat(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Not a date definition: {value}")


def file_type(suffix: str = None) -> Callable:
    """
    Factory function that provides additional option to the argparse directory type validation.

    :param suffix: Ensures that the file has the correct suffix

    :return: Function to validiates input for argparse
    """

    def file_type_func(string: str) -> Path:
        """
        Small helper function to convert input automatically to Path.
        Also raises an exception if the input is not a valid directory.

        :param string: Input argument

        :raises argparse.ArgumentTypeError: raised if invalid input

        :return: Input as pathlib.Path
        """

        # First Check: Directory must exist
        if Path(string).is_file():
            file_path = Path(string).absolute()
        else:
            raise argparse.ArgumentTypeError(f"Not a file: {string}")

        # Optional check: Directory path needs to end with certain directory name
        if file_path.suffix != suffix:
            raise argparse.ArgumentTypeError(f"{file_path} is not correct suffix: `{suffix}`")
        return file_path

    return file_type_func


def dirpath_type(string: str) -> Path:
    """
    Small helper function to convert input automatically to Path.
    Also raises an exception if the input is not a valid directory.

    :param string: Input argument

    :raises argparse.ArgumentTypeError:

    :return: Input as pathlib.Path
    """
    if Path(string).is_dir():
        return Path(string).absolute()
    else:
        raise argparse.ArgumentTypeError(f"Not a directory {string}")


def pysiral_settings_type(
        target: Literal["output", "proc"],
        level: str
) -> Callable:
    """
    Factory function that provides additional option to the argparse directory type validation.

    :param target: Target for the settings file, e.g. "output" or "proc"
    :param level: Level of processing definition, e.g. "l1", "l2", "l3"

    :return: Function to validiates input for argparse
    """

    def pysiral_settings_type_func(string: str) -> Path:
        """
        Small helper function to convert input automatically to Path.
        Also raises an exception if the input is not a valid directory.

        :param string: Input argument

        :raises argparse.ArgumentTypeError:

        :return: Input as pathlib.Path
        """
        if level not in PysiralProcessingLevels:
            raise argparse.ArgumentTypeError(
                f"Invalid processing level: {level} [{PysiralProcessingLevels.__members__}"
            )

        settings_filepath = psrlcfg.get_settings_file(target, level, string)
        if settings_filepath is None:
            msg = f"Invalid {target}:{level} settings filename or id: {string}\n"
            msg += f"Recognized {level} processor setting ids:\n"
            for procdef_id in psrlcfg.get_setting_ids(target, level):
                msg += f"  {procdef_id}\n"
            raise argparse.ArgumentTypeError(msg)
        return psrlcfg.get_settings_file(target, level, string)

    return pysiral_settings_type_func


def positive_int_type(value: str) -> int:
    """
    Convert a string to a positive integer.

    :param value: Input argument

    :raises argparse.ArgumentTypeError: if the input is not a valid positive integer

    :return: Positive integer
    """
    try:
        int_value = int(value)
        if int_value <= 0:
            raise argparse.ArgumentTypeError(f"Value must be a positive integer: {value}")
        return int_value
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid integer value: {value}")


def doi_type(value: str) -> str:
    """
    Convert a string to a valid DOI format.

    :param value: Input argument

    :raises argparse.ArgumentTypeError: if the input is not a valid DOI

    :return: Valid DOI string
    """
    regex_doi = re.compile(r"\b10\.\d{4,9}/[-.;()/:\w]+")
    match = regex_doi.findall(value)
    if not match:
        raise argparse.ArgumentTypeError(f"Invalid DOI format: {value}")
    return value


def pysiral_grid_id_type(value: str) -> str:
    """
    Convert a string to a valid pysiral grid ID.

    :param value: Input argument

    :raises argparse.ArgumentTypeError: if the input is not a valid grid ID

    :return: Valid grid ID string
    """
    valid_grid_ids = psrlcfg.get_setting_ids("grid")
    if value not in valid_grid_ids:
        raise argparse.ArgumentTypeError(f"Invalid grid ID: {value}. Valid IDs are: {', '.join(valid_grid_ids)}")
    return value
