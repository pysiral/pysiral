# -*- coding: utf-8 -*-

"""
"""

__author__ = "Stefan Hendricks"

__all__ = ["class_template", "clocks", "errorhandler", "flags"]

from loguru import logger


class DefaultLoggingClass(object):
    # TODO: Remove all instances of this class

    def __init__(self, name):
        self.log = logger
