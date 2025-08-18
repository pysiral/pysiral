#!/usr/bin/env python
# -*- coding: utf-8 -*-


import argparse

from pathlib import Path
from typing import List, Union

from dateperiods import DatePeriod
from loguru import logger

from pysiral import set_psrl_cpu_count
from pysiral.scripts.parser_items import DefaultCommandLineArguments
from pysiral.core.datahandler import L2iDataHandler
from pysiral.core.flags import DurationType, ProductProcessingLevels, DataRecordType
from pysiral.l3proc import (Level3GridDefinition, Level3OutputHandler,
                            Level3Processor, Level3ProductDefinition)
from pysiral.scripts.parser_items import (
    L2iDirectory, L3Settings, L3Grid, L3Outputs, L3Directory, ExcludeMonths,
    ProcessingPeriod, DOI, Duration, DataRecord, ProductProcessingLevel
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
        processing_level: ProductProcessingLevels = ProductProcessingLevels.LEVEL3_COLLATED,
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

        # Set product processing level based on the number of l2i directories
        if args.processing_level is not None:
            args.processing_level = (
                ProductProcessingLevels.LEVEL3_COLLATED if len(args.l2i_directory) == 1
                else ProductProcessingLevels.LEVEL3_SUPERCOLLATED
            )
            logger.info("Using processing level %s" % args.processing_level)

        # Set the l3 product directory if not provided
        if args.l3_directory is None:
            dirs = args.l2i_directory.parts + args.processing_level
            args.l3_directory = Path.joinpath(dirs)
            logger.info(f"Using default L3 directory: {args.l3_directory}")

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
            ProcessingPeriod(),
            L2iDirectory(required=True),
            # Mandatory arguments
            L3Outputs(required=True),
            L3Grid(required=True),
            # Optional arguments
            L3Directory(),
            Duration(),
            DataRecord(),
            ProductProcessingLevel(
                choices=[ProductProcessingLevels.LEVEL3_COLLATED,
                         ProductProcessingLevels.LEVEL3_SUPERCOLLATED],
                default=None
            ),
            ExcludeMonths(),
            DOI(),
        ]

        # create the parser
        parser = argparse.ArgumentParser(
            prog="pysiral l3proc",
            description="""
                    The Level-3 Processor (l3proc) generates Level-3 files (l3c|l3s) from l2i input files.
                    Level-3 files contain geophysical information and auxiliary data on spatio-temporal
                    grids. The processor uses a Level-3 product definition file to define the product metadata, 
                    the required Level-2 input data and Level-3 processor algorithm steps 
                    to be applied to either trajectory or gridded data. The output can be written
                    into multiple files.
                    """,
            epilog="For more information, see: https://pysiral.readthedocs.io",
            formatter_class=lambda prog: argparse.HelpFormatter(prog, width=96, indent_increment=4)  # noqa: E501
        )
        for arg_item in arg_item_list:
            arg_flags, args_dict = arg_item.get()
            parser.add_argument(*arg_flags, **args_dict)

        return parser
