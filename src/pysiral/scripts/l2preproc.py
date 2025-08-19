#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path
from typing import List, Union

from dateperiods import DatePeriod
from loguru import logger

from pysiral.core.datahandler import L2iDataHandler
from pysiral.l2preproc import Level2PreProcessor, Level2PreProcProductDefinition
from pysiral.scripts.parser_items import (
    ProcessingPeriod, ExcludeMonths, L2POutputs, L2iDirectory,
    DOI
)


def l2preproc(
    processing_period: DatePeriod = None,
    l2i_product_dir: Union[str, Path] = None,
    l2p_outputs: List[Union[str, Path]] = None,
    doi: str = None,
    exclude_month: List[int] = None,
) -> None:
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


class L2PreProcScriptArguments(object):

    def __init__(self) -> None:
        self.parser = self.get_argument_parser()

    def get(self, args_list: List[str] = None) -> "argparse.Namespace":
        return self.parser.parse_args() if args_list is None else self.parser.parse_args(args_list)

    @staticmethod
    def get_argument_parser() -> argparse.ArgumentParser:
        """
            Set up the command line argument parser for the Level-2 Processor.

            :return: The argument parser object.
            """

        # List of command line option required for the Level-1 pre-processor
        arg_item_list = [
            # Positional arguments
            L2iDirectory(nargs=None).as_positional("l2_product_directory"),
            L2POutputs(required=True).as_positional("l2p_output"),
            ProcessingPeriod(),
            # Optional arguments
            ExcludeMonths(),
            DOI(),
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
