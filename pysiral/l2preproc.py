# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 17:07:02 2017

@author: shendric
"""

from pysiral.logging import DefaultLoggingClass


class Level2PreProcessor(DefaultLoggingClass):

    def __init__(self):
        super(Level2PreProcessor, self).__init__(self.__class__.__name__)


class Level2PreProcProductDefinition(DefaultLoggingClass):

    def __init__(self):
        class_name = self.__class__.__name__
        super(Level2PreProcProductDefinition, self).__init__(class_name)
