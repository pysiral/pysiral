#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys
import time
from datetime import timedelta
from pathlib import Path
from typing import List, Union

from dateperiods import DatePeriod
from loguru import logger


from pysiral import psrlcfg, set_psrl_cpu_count
from pysiral.scripts.parser_items import DefaultCommandLineArguments
from pysiral.core.datahandler import L2iDataHandler
from pysiral.core.legacy_classes import DefaultLoggingClass, ErrorStatus
from pysiral.l2preproc import (Level2PreProcessor, Level2PreProcProductDefinition)
from pysiral.scripts.parser_items import (
    ProcessingPeriod, ExcludeMonths, L2POutputs, L2iDirectory
)



def l2preproc(
    processing_period: DatePeriod = None,
    l2i_product_dir: Union[str, Path] = None,
    l2p_outputs: List[Union[str, Path]] = None,
    doi: str = None,
    exclude_month: List[int] = None,
    mp_cpu_count: int = None,
    force_l2def_record_type: bool = False,
    **kwargs: dict
):
    """ Caller for converting Level-2 Intermediate (l2i) into
    Level-2 Pre-Processed (l2p) data products.
    NOTE: At the moment that only means summary of valid freeboard/thickness
          data points into daily summary files. """

    # Collect job settings from pysiral configuration data and


    # Get the product definition
    product_def = Level2PreProcProductDefinition()

    # Specifically add an output handler
    product_def.add_output_definition(
            l2i_product_dir,
            l2p_outputs,
            period="daily",
            doi=doi
    )

    # Prepare DataHandler
    # The l2 pre-processor requires l2i input files
    l2i_handler = L2iDataHandler(l2i_product_dir)

    # Get list of days for processing
    # start and/or stop can be omitted. In this case fall back to the
    # start and/or stop of l2i product availability
    breakpoint()
    days = processing_period.get_segments("day")
    if exclude_month is not None:
        days.filter_month(exclude_month)

    # Processor Initialization
    # NOTE: This is only for later cases. Not much is done here at this
    #       point
    l2_pre_processor = Level2PreProcessor(product_def)

#    # Loop over iterations (one per day)
    for day in days:

        # Do some extra logging
        logger.info("Processing Day [%s]" % day.label)

        # Get input files
        l2i_daily_files = l2i_handler.get_files_for_day(day.tcs.dt)
        if len(l2i_daily_files) == 0:
            logger.info("- no l2i products, skip day")
            continue
        logger.info("- Found %g l2i product files" % len(l2i_daily_files))

        # Process the orbits
        l2_pre_processor.process_l2i_files(l2i_daily_files, day)


class L2PreProcScriptArguments(DefaultLoggingClass):

    def __init__(self):
        super(L2PreProcScriptArguments, self).__init__(self.__class__.__name__)
        self.error = ErrorStatus()
        self._args = None

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
            L2iDirectory(),
            ProcessingPeriod(),
            L2POutputs(required=True),
            # Optional arguments
            ExcludeMonths(),
        ]

        # create the parser
        parser = argparse.ArgumentParser(
            prog="pysiral l2preproc",
            description="""
                       The Level-2 Pre-Processor (l2preproc) generates Level-2P files (l2p) from l2(i) 
                       input files. Level-2 files contain geophysical information and auxiliary data along the 
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

    @property
    def parser(self):
        # XXX: Move back to caller

        # Take the command line options from default settings
        # -> see config module for data types, destination variables, etc.
        clargs = DefaultCommandLineArguments()

        # List of command line option required for pre-processor
        # (argname, argtype (see config module), destination, required flag)
        options = [
            ("-start", "date", "start_date", False),
            ("-stop", "date", "stop_date", False),
            ("-l2i-product-dir", "l2i-product-dir", "l2i_product_dir", True),
            ("-l2p-output", "l2p-output", "l2p_output", False),
            ("-exclude-month", "exclude-month", "exclude_month", False),
            ("-doi", "doi", "doi", False),
            ("--remove-old", "remove-old", "remove_old", False),
            ("--no-critical-prompt", "no-critical-prompt",
             "no_critical_prompt", False),
            ("--no-overwrite-protection", "no-overwrite-protection",
             "overwrite_protection", False),
            ("--overwrite-protection", "overwrite-protection",
             "overwrite_protection", False)]

        # create the parser
        parser = argparse.ArgumentParser()
        for option in options:
            argname, argtype, destination, required = option
            argparse_dict = clargs.get_argparse_dict(
                argtype, destination, required)
            parser.add_argument(argname, **argparse_dict)

        parser.set_defaults(overwrite_protection=False)

        return parser
