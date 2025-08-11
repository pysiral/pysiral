# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 17:07:02 2017

@author: shendric
"""


from collections import OrderedDict
from pathlib import Path

from loguru import logger

from pysiral import psrlcfg
from pysiral.core import DefaultLoggingClass
from pysiral.core.errorhandler import ErrorStatus
from pysiral.core.output import Level2Output, OutputHandlerBase
from pysiral.l2data import L2iNCFileImport, Level2Data, Level2PContainer


class Level2PreProcessor(DefaultLoggingClass):

    def __init__(self, product_def):
        super(Level2PreProcessor, self).__init__(self.__class__.__name__)
        self.error = ErrorStatus()

        # Sanity check of product definition object
        if not isinstance(product_def, Level2PreProcProductDefinition):
            msg = f"Invalid Level-2 PreProcessor product definition: {type(product_def)}"
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
            try:
                l2i = L2iNCFileImport(l2i_file)
            except Exception as ex:
                msg = "Error (%s) in l2i file: %s"
                msg %= (ex, Path(l2i_file).name)
                logger.error(msg)
                continue
            l2p.append_l2i(l2i)

        # Merge the l2i object to a single L2Data object
        l2 = l2p.get_merged_l2()
        if l2 is None:
            logger.warning("- No valid freeboard data found for, skip day")
            return

        # Write output
        output = Level2Output(l2, self.job.output_handler)
        logger.info(f"- Wrote {self.job.output_handler.id} data file: {output.export_filename}")

    @property
    def job(self):
        return self._job


class Level2PreProcProductDefinition(DefaultLoggingClass):

    def __init__(self):
        self._output_handler = None
        class_name = self.__class__.__name__
        super(Level2PreProcProductDefinition, self).__init__(class_name)

    def add_output_definition(self, l2i_product_dir, output_def_file,
                              period="default", overwrite_protection=True,
                              doi=None):

        # Set given or default output handler
        self._output_handler = Level2POutputHandler(
            l2i_product_dir,
            output_def=output_def_file,
            period=period,
            doi=doi,
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
                 overwrite_protection=True, doi=None):
        # Fall back to default output if no output_def is given
        # (allows default initialization for the Level2 processor)
        if output_def == "default":
            output_def = self.default_output_def_filename
        super(Level2POutputHandler, self).__init__(output_def, applicable_data_level=2,
                                                   subfolder_tags=["year", "month"])
        self.error.caller_id = self.__class__.__name__
        logger.name = self.__class__.__name__
        self.l2i_product_dir = l2i_product_dir
        self._period = period
        self._doi = doi
        self.subdirectory = subdirectory
        self.overwrite_protection = overwrite_protection
        self._init_product_directory()

    def get_filename_from_data(self, l2p: Level2Data) -> str:
        """ Return the filename for a defined level-2 data object
        based on tag filenaming in output definition file """

        global filename_template
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
            msg %= (str(self._period), self.output_def_filename)
            self.error.add_error("invalid-outputdef", msg)
            self.error.raise_on_error()
        return self.fill_template_string(filename_template, l2p)

    def get_directory_from_data(self, l2p: Level2Data, create: bool = True) -> str:
        """
        Return the output directory based on l2 data object metadata

        :param l2p:
        :param create:
        :return:
        """
        l2p_period = l2p.period.tcs.dt
        directory = self._get_directory_from_dt(l2p_period)
        if create:
            self._create_directory(directory)
        self.error.raise_on_error()
        return directory

    def get_fullpath_from_data(self, l2):
        """ Return export path and filename based on information
        provided in the l2 data object """
        export_directory = self.get_directory_from_data(l2)
        export_filename = self.get_filename_from_data(l2)
        return Path(export_directory) / export_filename

    def get_global_attribute_dict(self, l2):
        attr_dict = OrderedDict()
        for attr_entry in self.output_def.global_attributes:
            attr_name, attr_template = zip(*attr_entry.items())
            attribute = self.fill_template_string(attr_template[0], l2)
            attr_dict[attr_name[0]] = attribute
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
        basedir = Path(self.l2i_product_dir).parent
        basedir = basedir / self.product_level_subfolder
        if self.subdirectory is not None:
            basedir = basedir / self.subdirectory
        if self.overwrite_protection:
            basedir = basedir / self.now_directory
        self._set_basedir(basedir)

    @property
    def has_doi(self):
        return self._doi is not None

    @property
    def doi(self):
        return self._doi

    @property
    def default_output_def_filename(self):
        local_settings_path = psrlcfg.pysiral_local_path
        return Path(local_settings_path) / Path(*self.default_file_location)
