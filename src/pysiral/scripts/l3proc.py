#!/usr/bin/env python
# -*- coding: utf-8 -*-


import argparse
import sys

from pathlib import Path
from typing import List, Union

from dateperiods import DatePeriod
from loguru import logger

from pysiral import psrlcfg, set_psrl_cpu_count
from pysiral.scripts.parser_items import DefaultCommandLineArguments
from pysiral.core.datahandler import L2iDataHandler
from pysiral.core.flags import DurationType, ProcessingLevels, DataRecordType
from pysiral.l3proc import (Level3GridDefinition, Level3OutputHandler,
                            Level3Processor, Level3ProductDefinition)
from pysiral.scripts.parser_items import (
    L2iDirectory, L3Settings, L3Grid, L3Output, L3Directory, ExcludeMonths,
    ProcessingPeriod, DOI, Duration, DataRecord, ProcessingLevel
)


def l3proc(
        processing_period: DatePeriod = None,
        l2i_product_directories: List[Path] = None,
        l3_product_directory: Union[str, Path] = None,
        l3_settings_file: Union[str, Path] = None,
        l3_grid_settings: Union[str, Path] = None,
        l3_output_definitions: List[Union[str, Path]] = None,
        duration: DurationType = DurationType.P1M,
        doi: str = None,
        processing_level: ProcessingLevels = ProcessingLevels.LEVEL3_COLLATED,
        data_record_type: DataRecordType = None,
        **kwargs: dict

) -> None:
    """

    :param processing_period:
    :param l2i_product_directories:
    :param l3_product_directory:
    :param l3_settings_file:
    :param l3_grid_settings:
    :param l3_output_definitions:
    :param duration:
    :param doi:
    :param processing_level:
    :param data_record_type:
    :return:
    """

    # --- Get the period segments for the Level-3 processor ---
    # NOTE: These depend on the chosen total time range and the duration period for the grid.
    if duration == "custom":
        period_segments = [processing_period]
        n_periods = 1
    else:
        period_segments = processing_period.get_segments(duration)
        n_periods = period_segments.n_periods

    # Get the output grid
    grid = Level3GridDefinition(l3_grid_settings)

    # Initialize the interface to the l2i products
    l2i_handler = L2iDataHandler(l2i_product_directories, search_str="l2")

    # Initialize the output handler
    # Currently,  overwrite protection is disabled per default
    output = []
    for l3_output_file in l3_output_definitions:
        output_handler = Level3OutputHandler(
            output_def=l3_output_file,
            base_directory=l3_product_directory,
            period=duration,
            doi=doi,
            data_record_type=processing_level,
            overwrite_protection=False
        )
        output.append(output_handler)

    # Compile the product def
    product_def = Level3ProductDefinition(l3_settings_file, grid, output, processing_period)

    # Initialize the Processor
    l3_processor = Level3Processor(product_def)

    # Loop over all iterations
    for i, time_range in enumerate(period_segments):

        # Report processing period
        msg = "# Processing %s period (%g of %g): %s"
        msg %= (processing_period, i+1, n_periods, time_range.date_label)
        logger.info(msg)

        # Retrieve files
        l2i_files = l2i_handler.get_files_from_time_range(time_range)
        logger.info("Num l2i files: %g" % len(l2i_files))
        if len(l2i_files) == 0:
            logger.info("Skip data period")
            continue

        # Start the Level-3 processing
        l3_processor.process_l2i_files(l2i_files, time_range)


class L3ProcScriptArguments(object):

    def __init__(self):
        self.parser = self.get_argument_parser()

    def get(self, args_list: List[str] = None) -> "argparse.Namespace":
        args = self.parser.parse_args() if args_list is None else self.parser.parse_args(args_list)
        if args.multiprocessing_num_cores is not None:
            set_psrl_cpu_count(args.multiprocessing_num_cores)
        return args

    @staticmethod
    def get_argument_parser() -> argparse.ArgumentParser:
        """
            Set up the command line argument parser for the Level-2 Processor.

            :return: The argument parser object.
            """

        # List of command line option required for the Level-1 pre-processor
        arg_item_list = [
            # Positional arguments
            L3Settings(),
            L2iDirectory(required=True),
            ProcessingPeriod(),
            # Mandatory arguments
            L3Output(required=True),
            L3Grid(required=True),
            # Optional arguments
            L3Directory(),
            Duration(),
            DataRecord(),
            ProcessingLevel(
                choices=[ProcessingLevels.LEVEL3_COLLATED, ProcessingLevels.LEVEL3_SUPERCOLLATED],
                default=ProcessingLevels.LEVEL3_COLLATED
            ),
            ExcludeMonths(),
            DOI(),
        ]

        # create the parser
        parser = argparse.ArgumentParser(
            prog="pysiral l2proc",
            description="""
                    The Level-2 Processor (l2proc) generates Level-2 files (l2/l2i) from l1p input files.
                    Level-2 files contain geophysical information and auxiliary data along the 
                    orbit at full sensor resolution. The processor uses a Level-2 product definition
                    file to define the product metadata, the list of auxiliary data files and the
                    algorithm steps to be applied to the Level-1P data. The output can be written
                    into multiple files.
                    Input Level-1 files are automatically selected based on the the source dataset ID 
                    and l1p version. (see also: `pysiral l2procfiles --help` for running the 
                    Level-2 processor on individual l1p files).
                    """,
            epilog="For more information, see: https://pysiral.readthedocs.io",
            formatter_class=lambda prog: argparse.HelpFormatter(prog, width=96, indent_increment=4)  # noqa: E501
        )
        for arg_item in arg_item_list:
            arg_flags, args_dict = arg_item.get()
            parser.add_argument(*arg_flags, **args_dict)

        return parser


class Level3ProcArgParser(object):

    def __init__(self):
        super(Level3ProcArgParser, self).__init__(self.__class__.__name__)

        self._args = None

    def parse_command_line_arguments(self):
        # use python module argparse to parse the command line arguments
        # (first validation of required options and data types)
        self._args = self.parser.parse_args()

        # Add additional check to make sure either `l1b-files` or
        # `start ` and `stop` are set

    #        l1b_file_preset_is_set = self._args.l1b_files_preset is not None
    #        start_and_stop_is_set = self._args.start_date is not None and \
    #            self._args.stop_date is not None
    #
    #        if l1b_file_preset_is_set and start_and_stop_is_set:
    #            self.parser.error("-start & -stop and -l1b-files are exclusive")
    #
    #        if not l1b_file_preset_is_set and not start_and_stop_is_set:
    #            self.parser.error("either -start & -stop or -l1b-files required")


    @property
    def parser(self):
        # XXX: Move back to caller

        # Take the command line options from default settings
        # -> see config module for data types, destination variables, etc.
        clargs = DefaultCommandLineArguments()

        # List of command line option required for pre-processor
        # (argname, argtype (see config module), destination, required flag)
        options = [
            ("-l2i-product-dir", "l2i-product-dir", "l2i_basedir", True),
            ("-l3-product-dir", "l3-product-dir", "l3_product_dir", False),
            ("-l3-settings", "l3-settings", "l3_settings", False),
            ("-l3-griddef", "l3-griddef", "l3_griddef", True),
            ("-l3-output", "l3-output", "l3_output", True),
            ("-start", "date", "start_date", True),
            ("-stop", "date", "stop_date", True),
            ("-period", "period", "period", False),
            ("-doi", "doi", "doi", False),
            ("-data-record-type", "data_record_type", "data_record_type", False),
            ("--remove-old", "remove-old", "remove_old", False),
            ("--no-critical-prompt", "no-critical-prompt",
             "no_critical_prompt", False)]

        # create the parser
        parser = argparse.ArgumentParser()
        for option in options:
            argname, argtype, destination, required = option
            argparse_dict = clargs.get_argparse_dict(
                argtype, destination, required)
            parser.add_argument(argname, **argparse_dict)

        return parser
