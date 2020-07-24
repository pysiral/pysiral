# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 21:25:26 2016

@author: Stefan

TODO: Move functionality to __init__

"""

from loguru import logger


class DefaultLoggingClass(object):

    def __init__(self, name):
        self.log = logger

    def info(self, *args, **kwargs):
        self.log.info(*args, **kwargs)

    def debug(self, *args, **kwargs):
        self.log.debug(*args, **kwargs)
