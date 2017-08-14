# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 17:07:02 2017

@author: shendric
"""

from pysiral.l2data import Level2PContainer, L2iNCFileImport
from pysiral.logging import DefaultLoggingClass
from pysiral.errorhandler import ErrorStatus


class Level2PreProcessor(DefaultLoggingClass):

    def __init__(self, product_def):
        super(Level2PreProcessor, self).__init__(self.__class__.__name__)
        self.error = ErrorStatus()

        # Sanity check of product definition object
        if not isinstance(product_def, Level2PreProcProductDefinition):
            msg = "Invalid Level-2 PreProcessor product definition: %s" % \
                type(product_def)
            self.error.add_error("invalid-l2preproc-def", msg)
            self.error.raise_on_error()
        self._job = product_def

    def process_l2i_files(self, l2i_files):
        """ Reads all l2i files and merges the valid data into a l2p
        summary file """

        # l2p: Container for storing l2i objects
        l2p = Level2PContainer()

        # Add all l2i objects to the l2p container.
        # NOTE: Only memory is the limit
        for l2i_file in l2i_files:
            l2i = L2iNCFileImport(l2i_file)
            l2p.append_l2i(l2i)

        # Merge the l2i object to a single L2Data object
        l2 = l2p.get_merged_l2()

    @property
    def job(self):
        return self._job


class Level2PreProcProductDefinition(DefaultLoggingClass):

    def __init__(self):
        class_name = self.__class__.__name__
        super(Level2PreProcProductDefinition, self).__init__(class_name)

    def add_output_definition(self, output_def_file,
                              overwrite_protection=True):
        pass
