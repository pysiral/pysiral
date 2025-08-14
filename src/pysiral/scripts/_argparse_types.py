# -*- coding: utf-8 -*-

"""
"""

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"

import argparse
from dateperiods import DatePeriod
from datetime import date, datetime
from pathlib import Path
from typing import Callable, List

from pysiral import psrlcfg
from pysiral.core.flags import BasicProcessingLevels


def dir_type(value: str) -> Path:
    # First Check: Directory must exist
    if Path(value).is_dir():
        return Path(value).absolute()
    else:
        raise argparse.ArgumentTypeError(f"Not a directory: {value}")


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


def pysiral_procdef_type(level: BasicProcessingLevels) -> Callable:
    """
    Factory function that provides additional option to the argparse directory type validation.

    :param level: Level of processing definition, e.g. "l1", "l2", "l3"

    :return: Function to validiates input for argparse
    """

    def procdef_type_func(string: str) -> Path:
        """
        Small helper function to convert input automatically to Path.
        Also raises an exception if the input is not a valid directory.

        :param string: Input argument

        :raises argparse.ArgumentTypeError:

        :return: Input as pathlib.Path
        """
        if not level in BasicProcessingLevels:
            raise argparse.ArgumentTypeError(f"Invalid processing level: {level} [{BasicProcessingLevels.__members__}")

        settings_filepath = psrlcfg.get_settings_file("proc", level, string)
        if settings_filepath is None:
            msg = f"Invalid {level} settings filename or id: {string}\n"
            msg += f"Recognized {level} processor setting ids:\n"
            for procdef_id in psrlcfg.get_setting_ids("proc", level):
                msg += f"  {procdef_id}\n"
            raise argparse.ArgumentTypeError(msg)
        return psrlcfg.get_settings_file("proc", level, string)

    return procdef_type_func


def proc_period_type(arg_string: str) -> DatePeriod:
    """
    Checks that a string is a valid date period start/end definitions [YYYY-MM-DD or YYYY-MM].

    :param arg_string: Input argument in the format "YYYY-MM[-DD]"

    :raises argparse.ArgumentTypeError: if the input is not a valid period definition

    :return: string object
    """
    # Note: Number of arguments should already be handled by `required_length` action in the parser.
    #       This we can safely split at whitespace and expect two or one date arguments.
    args = arg_string.split()

    # If only one date is provided, use it as both start and stop date.
    if len(args) == 1:
        args.append(args[0])

    def parse_date(date_str: str) -> List[int]:
        try:
            return [int(part) for part in date_str.split("-")]
        except ValueError:
            raise argparse.ArgumentTypeError(f"Invalid date format: {date_str}. Expected YYYY-MM-DD.")

    date_parts = [parse_date(date_str) for date_str in args]
    try:
        return DatePeriod(*date_parts)
    except ValueError as e:
        raise argparse.ArgumentTypeError(f"Invalid date period format: {arg_string}. Error: {e}")


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
