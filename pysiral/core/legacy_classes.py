# -*- coding: utf-8 -*-

"""
This module provides access to deprecated legacy classes and functions that are still used in some parts of the codebase.
These classes and functions are maintained for backward compatibility but should not be used in new development.
"""

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"


from collections import UserDict


class AttrDict(UserDict):
    """
    Short implementation of attrdict.AttrDict using UserDict. The code is based on the solutions shared here
    https://stackoverflow.com/a/76231823 and has been modified to allow nestesd AttrDict instances.
    """

    def __getattr__(self, key):
        item = self.__getitem__(key)
        return AttrDict(**item) if isinstance(item, dict) else item

    def __setattr__(self, key, value):
        if key == "data":
            return super().__setattr__(key, value)
        return self.__setitem__(key, value)