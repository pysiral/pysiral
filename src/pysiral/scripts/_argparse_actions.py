# -*- coding: utf-8 -*-

"""

"""

import argparse
from datetime import datetime
from typing import Type, List


from dateperiods import DatePeriod

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

        def __call__(self, parser, args, values: List[str], option_string=None) -> None:

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
