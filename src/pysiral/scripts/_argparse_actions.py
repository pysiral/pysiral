# -*- coding: utf-8 -*-

"""

"""

import argparse
from typing import Type

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
