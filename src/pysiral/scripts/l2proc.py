#!/usr/bin/env python
# -*- coding: utf-8 -*-


import argparse
import glob
import re
import time
from typing import List, Union
from pathlib import Path
from datetime import timedelta

from dateperiods import DatePeriod
from loguru import logger

from pysiral import psrlcfg, set_psrl_cpu_count
from pysiral.scripts.parser_items import DefaultCommandLineArguments
from pysiral.core.datahandler import L1PDataHandler
from pysiral.l2proc import Level2Processor, Level2ProductDefinition

from pysiral.scripts.parser_items import (
    ProcessingPeriod, ExcludeMonths, Hemisphere, PlatformID,
    L2Settings, L2Outputs, SourceDatasetID, MultiProcesssingNumCores,
    UseMultiProcesssing, ForceL2DefRecordType
)


def l2proc(
    processing_period: DatePeriod = None,
    l2_settings: Union[str, Path] = None,
    l2_outputs: List[Union[str, Path]] = None,
    exclude_month: List[int] = None,
    input_version: str = None,
    l1p_version: str = None,
    mp_cpu_count: int = None,
    force_l2def_record_type: bool = False,
    **kwargs: dict
) -> None:
    """
    Main entry point for the Level-2 Processor job.

    This function initializes the Level-2 processor with the provided settings
    and processes the Level-1P data files for the specified time range.

    :param processing_period: Start date for processing in [year, month, [day]] format.
    :param l2_settings: Path to the Level-2 settings file or its identifier.
    :param l2_outputs: List of output definitions for Level-2 products.
    :param exclude_month: List of months to exclude from processing (1-12).
    :param input_version: Version of the input data (optional).
    :param l1p_version: Version of the Level-1P data (optional).
    :param mp_cpu_count: Number of CPUs to use for multiprocessing (optional).
    :param force_l2def_record_type: If True, forces the use of a specificLevel-2 definition record type.
    """

    # Update pysiral multiprocessing settings
    if mp_cpu_count is not None:
        logger.info(f"Using {mp_cpu_count} CPU cores.")
        set_psrl_cpu_count(mp_cpu_count)

    # Get the Level-2 product definition containing the product metadata, list of
    # auxiliary data files and the algorithm steps. This information is saved in
    # yaml L2 product definition files (`pysiral/resources/pysiral-cfg/proc/l2/`).
    product_def = Level2ProductDefinition(
        l2_settings,
        force_l2def_record_type=force_l2def_record_type
    )
    platform = product_def.l2def.metadata.platform
    hemisphere = product_def.l2def.metadata.hemisphere

    # The Level-2 data can be written to one or multiple output files.
    for l2_output in l2_outputs:
        product_def.add_output_definition(l2_output)

    # Clip the time range to the valid time range of the target platform
    period = processing_period.intersect(psrlcfg.get_platform_period(platform))
    if period is None:
        msg = f"Invalid period definition ({processing_period.label}) for platform {platform}"
        raise ValueError(msg)

    # The Level-2 processor operates in monthly segments, and for a given period,
    # certain months can be excluded from processing.
    period_segments = period.get_segments("month", crop_to_period=True)
    if exclude_month is not None:
        period_segments.filter_month(exclude_month)

    # Initialize the Level-1P data handler
    l1b_data_handler = L1PDataHandler(
        platform,
        hemisphere,
        source_version=input_version,
        file_version=l1p_version
    )

    # Processor Initialization
    logger.info("Level-2 processor: initializing")
    l2_processor = Level2Processor(product_def)

    # Loop over monthly segments of the period
    for time_range in period_segments:

        # Do some extra logging
        logger.info(f"Processing period: {time_range.label}")

        # Get input files
        l1b_files = l1b_data_handler.get_files_from_time_range(time_range)
        logger.info("Found %g files in %s" % (len(l1b_files), l1b_data_handler.last_directory))

        # Process the orbits
        l2_processor.process_l1b_files(l1b_files)

    # All done
    logger.info("Level-2 processor: completed")


class L2ProcScriptArguments(object):

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
            L2Settings(),
            ProcessingPeriod(),
            L2Outputs(required=True),
            # Optional arguments
            ExcludeMonths(),
            SourceDatasetID(),
            UseMultiProcesssing(),
            MultiProcesssingNumCores(),
            ForceL2DefRecordType()
        ]

        # create the parser
        parser = argparse.ArgumentParser(
            prog="pysiral l2proc",
            description="""
                    The Level-2 Processor (l2roc) generates Level-2 files (l2/l2i) from l1p input files.
                    Level-2 files contain geophysical information and auxiliary data along the 
                    orbit at full sensor resolution. The processor uses a Level-2 product definition
                    file to define the product metadata, the list of auxiliary data files and the
                    algorithm steps to be applied to the Level-1P data. The output can be written
                    into multiple files. 
                    """,
            epilog="For more information, see: https://pysiral.readthedocs.io",
            formatter_class=lambda prog: argparse.HelpFormatter(prog, width=96, indent_increment=4)  # noqa: E501
        )
        for arg_item in arg_item_list:
            arg_flags, args_dict = arg_item.get()
            parser.add_argument(*arg_flags, **args_dict)

        return parser
