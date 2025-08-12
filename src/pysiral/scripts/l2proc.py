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
from pysiral.core.legacy_classes import DefaultLoggingClass, ErrorStatus
from pysiral.l2proc import Level2Processor, Level2ProductDefinition


def l2proc(
    start_date: List[int] = None,
    stop_date: List[int] = None,
    l2_settings: Union[str, Path] = None,
    l2_outputs: List[Union[str, Path]] = None,
    exclude_month: List[int] = None,
    input_version: str = None,
    l1p_version: str = None,
    mp_cpu_count: int = None,
    force_l2def_record_type: bool = False,
    output_directory: Union[str, Path] = None,
    **kwargs: dict
) -> None:
    """
    Main entry point for the Level-2 Processor job.

    This function initializes the Level-2 processor with the provided settings
    and processes the Level-1P data files for the specified time range.

    :param start_date: Start date for processing in [year, month, [day]] format.
    :param stop_date: Stop date for processing in [year, month, [day]] format
    :param l2_settings: Path to the Level-2 settings file or its identifier.
    :param l2_outputs: List of output definitions for Level-2 products.
    :param exclude_month: List of months to exclude from processing (1-12).
    :param input_version: Version of the input data (optional).
    :param l1p_version: Version of the Level-1P data (optional).
    :param mp_cpu_count: Number of CPUs to use for multiprocessing (optional).
    :param force_l2def_record_type: If True, forces the use of a specificLevel-2 definition record type.
    :param output_directory: Directory where the output files will be saved.
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

    # get the period for the Level-2 Processor
    period = DatePeriod(start_date, stop_date)

    # Clip the time range to the valid time range of the target platform
    period = period.intersect(psrlcfg.get_platform_period(platform))
    if period is None:
        msg = f"Invalid period definition ({start_date}-{stop_date}) for platform {platform}"
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


class Level2ProcArgParser(DefaultLoggingClass):

    def __init__(self):
        super(Level2ProcArgParser, self).__init__(self.__class__.__name__)
        self.error = ErrorStatus()
        self._args = self.get_parser_args()

    def get_parser_args(self) -> argparse.ArgumentParser:
        # XXX: Move back to caller

        # Take the command line options from default settings
        # -> see config module for data types, destination variables, etc.
        clargs = DefaultCommandLineArguments()

        # List of command line option required for pre-processor
        # (argname, argtype (see config module), destination, required flag)
        options = [
            ("-l2-settings", "l2-settings", "l2_settings", True),
            ("-start", "date", "start_date", False),
            ("-stop", "date", "stop_date", False),
            ("-exclude-month", "exclude-month", "exclude_month", False),
            ("-input-version", "input-version", "input_version", False),
            ("-l1p-version", "l1p-version", "l1p_version", False),
            ("-l2-output", "l2-output", "l2_output", False),
            ("-mp-cpu-count", "mp-cpu-count", "mp_cpu_count", False),
            ("--force-l2def-record-type", "force-l2def-record-type", "force_l2def_record_type", False),
        ]

        # create the parser
        parser = argparse.ArgumentParser()
        for option in options:
            argname, argtype, destination, required = option
            argparse_dict = clargs.get_argparse_dict(argtype, destination, required)
            parser.add_argument(argname, **argparse_dict)
        parser.set_defaults(overwrite_protection=False)

        return parser

    @property
    def l2_settings_file(self):
        l2_settings = self._args.l2_settings
        filename = self.pysiral_config.get_settings_file("proc", "l2", l2_settings)
        if filename is not None:
            return filename
        msg = "Invalid l2 settings filename or id: %s\n" % l2_settings
        msg = msg + " \nRecognized Level-2 processor setting ids:\n"
        for l2_settings_id in psrlcfg.get_setting_ids("proc", "l2"):
            msg = f'{msg}  {l2_settings_id}' + "\n"
        self.error.add_error("invalid-l2-settings", msg)
        self.error.raise_on_error()

    @property
    def l2_output(self):
        filenames = []
        for l2_output in self._args.l2_output.split(";"):
            filename = psrlcfg.get_settings_file("output", "l2i", l2_output)

            if filename is None:
                msg = "Invalid l2 outputdef filename or id: %s\n" % l2_output
                msg = msg + " \nRecognized Level-2 output definitions ids:\n"
                l2_output_ids = psrlcfg.get_setting_ids("output", "l2i")
                for l2_output_id in l2_output_ids:
                    msg = f'{msg}    - {l2_output_id}' + "\n"
                self.error.add_error("invalid-l2-outputdef", msg)
                self.error.raise_on_error()
            else:
                filenames.append(filename)

        if not filenames:
            msg = "No valid output definition file found for argument: %s"
            msg %= (str(self._args.l3_output))
            self.error.add_error("invalid-outputdef", msg)
            self.error.raise_on_error()

        return filenames
