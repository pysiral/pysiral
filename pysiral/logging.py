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
