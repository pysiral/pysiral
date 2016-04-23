# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 21:25:26 2016

@author: Stefan
"""

from logbook import Logger, StreamHandler
import sys


def stdout_logger(name):
    StreamHandler(sys.stdout).push_application()
    log = Logger(name)
    return log


class DefaultLoggingClass(object):

    def __init__(self, name):
        self.log = stdout_logger(name)

    def info(self, *args, **kwargs):
        self.log.info(*args, **kwargs)

    def debug(self, *args, **kwargs):
        self.log.debug(*args, **kwargs)
