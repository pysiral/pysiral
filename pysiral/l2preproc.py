# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 17:07:02 2017

@author: shendric
"""

from pysiral.config import ConfigInfo
from pysiral.errorhandler import ErrorStatus
from pysiral.l2data import Level2PContainer, L2iNCFileImport
from pysiral.logging import DefaultLoggingClass
from pysiral.output import Level2Output, OutputHandlerBase

import os


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

    def process_l2i_files(self, l2i_files, period):
        """ Reads all l2i files and merges the valid data into a l2p
        summary file """

        # l2p: Container for storing l2i objects
        l2p = Level2PContainer(period)

        # Add all l2i objects to the l2p container.
        # NOTE: Only memory is the limit
        for l2i_file in l2i_files:
            l2i = L2iNCFileImport(l2i_file)
            l2p.append_l2i(l2i)

        # Merge the l2i object to a single L2Data object
        l2 = l2p.get_merged_l2()

        # Write output
        output = Level2Output(l2, self.job.output_handler)
        self.log.info("- Wrote %s data file: %s" % (
                self.job.output_handler.id, output.export_filename))

    @property
    def job(self):
        return self._job


class Level2PreProcProductDefinition(DefaultLoggingClass):

    def __init__(self):
        class_name = self.__class__.__name__
        super(Level2PreProcProductDefinition, self).__init__(class_name)

    def add_output_definition(self, l2i_product_dir, output_def_file,
                              period="default", overwrite_protection=True):
                # Set given or default output handler
        self._output_handler = Level2POutputHandler(
            l2i_product_dir,
            output_def=output_def_file,
            period=period,
            overwrite_protection=overwrite_protection)

    @property
    def output_handler(self):
        return self._output_handler


class Level2POutputHandler(OutputHandlerBase):
    """ Default output handler with pysiral conventions. Uses product
    directory from local_machine_def.yaml as standard repository """

    # Some fixed parameters for this class
    default_file_location = ["settings", "outputdef", "l2p_default.yaml"]
    subfolder_tags = ["year", "month"]
    applicable_data_level = 2

    def __init__(self, l2i_product_dir, output_def="default",
                 subdirectory=None, period="default",
                 overwrite_protection=True):
        # Fall back to default output if no output_def is given
        # (allows default initialization for the Level2 processor)
        if output_def == "default":
            output_def = self.default_output_def_filename
        super(Level2POutputHandler, self).__init__(output_def)
        self.error.caller_id = self.__class__.__name__
        self.log.name = self.__class__.__name__
        self.l2i_product_dir = l2i_product_dir
        self._period = period
        self.subdirectory = subdirectory
        self.overwrite_protection = overwrite_protection
        self._init_product_directory()

    def get_filename_from_data(self, l2p):
        """ Return the filename for a defined level-2 data object
        based on tag filenaming in output definition file """

        try:
            template_ids = self.output_def.filenaming.keys()
            period_id = self._period
            # Fall back to default if no filenaming convention for given
            # data period
            if period_id not in template_ids:
                period_id = "default"
            filename_template = self.output_def.filenaming[period_id]
        except AttributeError:
            filename_template = self.output_def.filenaming
        except KeyError:
            msg = "Missing filenaming convention for period [%s] in [%s]"
            msg = msg % (str(self._period), self.output_def_filename)
            self.error.add_error("invalid-outputdef", msg)
            self.error.raise_on_error()
        filename = self.fill_template_string(filename_template, l2p)
        return filename

    def get_directory_from_data(self, l2p, create=True):
        """ Return the output directory based on information provided
        in an l2 data object """
        directory = self._get_directory_from_dt(l2p.info.start_time)
        if create:
            self._create_directory(directory)
        self.error.raise_on_error()
        return directory

    def get_fullpath_from_data(self, l2):
        """ Return export path and filename based on information
        provided in the l2 data object """
        export_directory = self.get_directory_from_l2(l2)
        export_filename = self.get_filename_from_l2(l2)
        return os.path.join(export_directory, export_filename)

    def get_global_attribute_dict(self, l2):
        attr_dict = {}
        for attr_name in self.output_def.global_attributes.iterkeys():
            attr_template = self.output_def.global_attributes[attr_name]
            attribute = self.fill_template_string(attr_template, l2)
            attr_dict[attr_name] = attribute
        return attr_dict

    def remove_old(self, time_range):
        """ This method will erase all files in the target orbit for a
        given time range. Use with care """
        # Get the target directory
        # XXX: Assumption time_range is monthly
        pass

    def _init_product_directory(self):
        """ Get main product directory from local_machine_def, add mandatory
        runtag subdirectory, optional second subdirectory for overwrite
        protection and product level id subfolder"""
        basedir = self.l2i_product_dir
        basedir, l2i_subdir = os.path.split(basedir)
        basedir = os.path.join(basedir, self.product_level_subfolder)
        if self.subdirectory is not None:
            basedir = os.path.join(basedir, self.subdirectory)
        if self.overwrite_protection:
            basedir = os.path.join(basedir, self.now_directory)
        self._set_basedir(basedir)

    @property
    def default_output_def_filename(self):
        pysiral_config = ConfigInfo()
        local_settings_path = pysiral_config.pysiral_local_path
        return os.path.join(local_settings_path, *self.default_file_location)
