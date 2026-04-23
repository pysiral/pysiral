#!/usr/bin/env python
# -*- coding: utf-8 -*-


import argparse

from pathlib import Path
from typing import List, Union

from dateperiods import DatePeriod
from loguru import logger

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
        l2_product_directory: List[Path] = None,
        l3_product_directory: Union[str, Path] = None,
        l3_settings: Union[str, Path] = None,
        l3_grid_id: str = None,
        l3_output: List[Union[str, Path]] = None,
        duration: DurationType = DurationType.P1M,
        doi: str = None,
        processing_level: ProductProcessingLevels = ProductProcessingLevels.LEVEL3_COLLATED,
        data_record: DataRecordType = None,
        exclude_months: List[int] = None,
) -> None:
    """

    :param processing_period:
    :param l2_product_directory:
    :param l3_product_directory:
    :param l3_settings:
    :param l3_grid_id:
    :param l3_output:
    :param duration:
    :param doi:
    :param processing_level:
    :param data_record:
    :param exclude_months:
        List of months to exclude from processing (1-12).
        If None, no months are excluded.
        If provided, the months will be filtered out from the processing period.
    :return:
    """

    # --- Get the period segments for the Level-3 processor ---
    # NOTE: These depend on the chosen total time range and the duration period for the grid.
    if duration == "custom":
        period_segments = [processing_period]
        n_periods = 1
    else:
        period_segments = processing_period.get_segments(duration)
        if exclude_months is not None:
            period_segments.filter_month(exclude_months)
        n_periods = period_segments.n_periods

    # Get the output grid
    grid = Level3GridDefinition.from_grid_id(l3_grid_id)

    # Initialize the interface to the l2i products
    l2i_handler = L2iDataHandler(l2_product_directory)

    # Initialize the output handler
    # Currently,  overwrite protection is disabled per default
    output = []
    for l3_output_file in l3_output:
        output_handler = Level3OutputHandler(
            output_def=l3_output_file,
            base_directory=l3_product_directory,
            period=duration,
            doi=doi,
            data_record_type=data_record,
            overwrite_protection=False
        )
        output.append(output_handler)

    # Compile the product def
    product_def = Level3ProductDefinition(l3_settings, grid, output, processing_period)

    # Initialize the Processor
    l3_processor = Level3Processor(product_def)

    # Loop over all iterations
    for i, time_range in enumerate(period_segments):

        # Report processing period
        msg = "# Processing %s period (%g of %g): %s"
        msg %= (processing_period.date_label, i+1, n_periods, time_range.date_label)
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

        # Parse the command line arguments
        args = self.parser.parse_args() if args_list is None else self.parser.parse_args(args_list)

        # Add variable processing level (L3C for single L2(i) directory, L3S for multiple)
        args.processing_level = self.autocomplete_product_processing_level(args.l2_product_directory)

        # Autocomplete the L3 product directory if not provided
        args.l3_product_directory = self.autocomplete_l3_product_directory(
            args.l3_product_directory,
            args.l2_product_directory,
            args.processing_level,
        )

        # Autocomplete the L3 product directory if not provided
        args.duration = self.autocomplete_l3_product_duration(args.duration, args.processing_period)

        return args

    @staticmethod
    def autocomplete_product_processing_level(l2_product_directory: Union[List[Path]]) -> str:
        """
        Autocomplete the product processing level based on the provided arguments.

        :param l2_product_directory: The L2 product directory path or multiple thereof.

        :return: product processing level as a string.
        """

        # Check if the l2_product_directory is a list of Path objects
        assert (
            isinstance(l2_product_directory, (list, Path)),
            "l2_product_directory must be a list of Path objects"
        )

        # Set the processing level based on the number of L2(i) directories
        processing_level = (
            ProductProcessingLevels.LEVEL3_COLLATED if len(l2_product_directory) == 1
            else ProductProcessingLevels.LEVEL3_SUPERCOLLATED
        )
        logger.info(f"Set output product processing level to `{processing_level}`")
        return processing_level

    @staticmethod
    def autocomplete_l3_product_directory(
            l3_product_directory: Union[str, Path, None],
            l2_product_directory: Union[str, Path, List[Path]],
            processing_level: str,
    ) -> Path:
        """
        Autocomplete the L3 product directory based on the provided arguments.

        :param l3_product_directory: The L3 product directory path or None if not provided.
        :param l2_product_directory: The L2 product directory path or multiple thereof.
        :param processing_level: The target processing level.

        :return: Updated Level-3 product directory path.
        """
        # Set the l3 product directory if not provided
        if l3_product_directory is not None:
            return l3_product_directory

        # Check if processing level is L3C. Only in this case the L3 product directory
        # can be derived from the L2(i) product directory. Otherwise, it does not fit in the
        # single platform product directory structure. Hence, the expected setting of the
        # L3S product directory is requested.
        if processing_level == ProductProcessingLevels.LEVEL3_SUPERCOLLATED:
            raise argparse.ArgumentTypeError(
                """
                The L3 product directory (`-O/--l3-output`) has not been specified, but is mandatory if 
                multiple L2(i) directories are provided."
                """
            )
        dirs = list(l2_product_directory[0].parts[:-1]) + [processing_level]
        l3_product_directory = Path().joinpath(*dirs)
        logger.info(f"Set L3 product directory to: `{l3_product_directory}`")
        return l3_product_directory

    @staticmethod
    def autocomplete_l3_product_duration(
            duration: Union[str, None],
            processing_period: DatePeriod,
    ) -> str:
        """
        Autocomplete the L3 product directory based on the provided arguments.

        :param duration:
        :param processing_period:

        :return: Updated arguments namespace with the L3 product directory set (if required)
        """
        if duration is not None:
            return duration
        duration_dict = {"Y": DurationType.P1M, "M": DurationType.P1M, "D": DurationType.P1D}
        least_significant_period = processing_period.duration.isoformat[-1]  # Remove the trailing 'Z'
        duration = duration_dict.get(least_significant_period)
        logger.info(f"Set L3 product duration to: `{duration}`")
        return duration

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
            L3Grid(),
            ProcessingPeriod(),
            L2iDirectory(required=True),
            # Mandatory arguments
            L3Outputs(required=True),
            # Optional arguments
            L3Directory(),
            Duration(),
            DataRecord(),
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
