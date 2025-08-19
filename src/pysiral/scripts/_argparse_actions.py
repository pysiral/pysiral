# -*- coding: utf-8 -*-

"""

"""

import argparse
from datetime import datetime
from pathlib import Path
from typing import Type, Literal, Union


from dateperiods import DatePeriod
from pysiral import psrlcfg
from pysiral.core.flags import PysiralProcessingLevels, ProductProcessingLevels

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"


def required_length(nmin: int, nmax: int) -> Type[argparse.Action]:
    """
    argparse action that defines the required number of arguments for an argument
    with nargs `+` or `*`. The action will raise an error if the number of
    arguments is not within the specified range.

    :param nmin: minimum number of arguments required
    :param nmax: maximum number of arguments allowed

    :return: a subclass of `argparse.Action` that checks the number of arguments
    """

    class RequiredLength(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if not nmin <= len(values) <= nmax:
                msg = f"argument {self.dest} requires between {nmin} and {nmax} arguments"
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)

    return RequiredLength


def is_valid_date(date_str: str) -> bool:
    """
    Check if a date string is valid in the format YYYY-MM[-DD].

    :param date_str: The string to check, e.g. "2020-01" or "2020-01-15".

    :return: Boolean indicating whether the date string is valid.
    """
    try:
        date_nums = [int(part) for part in date_str.split("-")]
        if len(date_nums) == 2:
            date_nums.append(1)  # Default to the first day of the month
        assert isinstance(datetime(*date_nums), datetime)
        return True
    except (AssertionError, TypeError):
        return False


def period_conversion() -> Type[argparse.Action]:
    """
    Argparse action to convert a list of datetime objects into a start and end time coverage.

    :return: argparse.Action subclass
    """

    class PeriodConversion(argparse.Action):
        """
        Converts a period definition to a `dateperiods.DatePeriod` object.

        E.g.

            `2020-01`             -> `DatePeriod(2020-01-01, 2020-01-31)`
            `2020-01 2020-03`     -> `DatePeriod(2020-01-01, 2020-03-31)`
            `2020-01-10 2020-03`  -> `DatePeriod(2020-01-10, 2020-01-31)`
            etc.
        """

        def __call__(self, parser, args, value: str, option_string=None) -> None:

            values = value.split(":")

            # Check if either 1 or 2 values are provided
            if not 1 <= len(values) <= 2:
                msg = f"argument {self.dest} requires between 1 (start only) and 2 (start & end) date definitions"
                raise argparse.ArgumentTypeError(msg)

            # Check if the values are valid date strings
            for value in values:
                if not is_valid_date(value):
                    msg = f"Invalid date format: {value}. Expected YYYY-MM[-DD]."
                    raise argparse.ArgumentTypeError(msg)

            # Autocomplete the end time if only start time is provided
            if len(values) == 1:
                # If only one value is provided, set the end time to the start time
                values.append(values[0])

            # Set the argument value as a DatePeriod object
            setattr(
                args,
                self.dest,
                DatePeriod(
                    [int(v) for v in values[0].split("-")],
                    [int(v) for v in values[1].split("-")]
                )
            )

    return PeriodConversion


def pysiral_settings_action(
        target: Literal["proc", "output"],
        level: Union[PysiralProcessingLevels, ProductProcessingLevels] = None
) -> Type[argparse.Action]:

    class PysiralSettingsLookup(argparse.Action):
        """
        Validates pysiral settings file input for argparse and resolves
        file path if required for one or multiple settings definitions.
        """

        def __call__(self, parser, args, value: str, option_string=None) -> None:
            """
            Ensure that the argument is the full file path to the configuration file.
            Valid inputs are either a file path (taken as is) or a settings definition id
            that will be resolved using pysiral package configuration.

            The input type is expected to be a string or a list of strings (if nargs='+' is used).

            :param parser:
            :param args:
            :param value:

            :param option_string:

            :return:
            """

            # Check if one or multiple settings definitions are provided
            if is_single_setting := isinstance(value, str):
                value = [value]  # Convert to list for uniform processing

            # Resolve the settings file path (if path)
            # NOTE: This will raise an exception if the input is not a valid settings file
            value = [self._resolve_file(v) if isinstance(v, str) else v for v in value]

            # Set attribute in args
            setattr(args, self.dest, value[0] if is_single_setting else value)

        @staticmethod
        def _resolve_file(string: str) -> Path:
            """
            Small helper function to convert input automatically to Path.
            Also raises an exception if the input is not a valid directory.

            :param string: Input argument

            :raises argparse.ArgumentTypeError:

            :return: Input as pathlib.Path
            """
            if level not in [*PysiralProcessingLevels, *ProductProcessingLevels]:
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

    return PysiralSettingsLookup
