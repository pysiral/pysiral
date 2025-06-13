#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys
import time
from datetime import timedelta
from pathlib import Path

from dateperiods import DatePeriod
from loguru import logger

from pysiral import psrlcfg
from pysiral.core import DefaultLoggingClass
from pysiral.core.config import DefaultCommandLineArguments
from pysiral.core.datahandler import L2iDataHandler
from pysiral.core.errorhandler import ErrorStatus
from pysiral.l2preproc import (Level2PreProcessor,
                               Level2PreProcProductDefinition)


def pysiral_l2preproc():
    """ Caller for converting Level-2 Intermediate (l2i) into
    Level-2 Pre-Processed (l2p) data products.
    NOTE: At the moment that only means summary of valid freeboard/thickness
          data points into daily summary files. """

    # Collect job settings from pysiral configuration data and
    # command line arguments
    args = Level2PreProcArgParser()

    # Parse and validate the command line arguments
    args.parse_command_line_arguments()

    # Get confirmation for critical choices (if necessary)
    args.critical_prompt_confirmation()

    # Start the level-2 pre-processor
    # Get start time of processor run
    t0 = time.process_time()

    # Get the product definition
    product_def = Level2PreProcProductDefinition()

    # Specifically add an output handler
    product_def.add_output_definition(
            args.l2i_product_dir,
            args.l2p_output,
            period="daily",
            doi=args.doi,
            overwrite_protection=args.overwrite_protection)

    # Prepare DataHandler
    # The l2 pre-processor requires l2i input files
    l2i_handler = L2iDataHandler(args.l2i_product_dir)

    # Get list of days for processing
    # start and/or stop can be omitted. In this case fall back to the
    # start and/or stop of l2i product availability
    start = args.start if args.start is not None else l2i_handler.start_month
    stop = args.stop if args.stop is not None else l2i_handler.stop_month
    period = DatePeriod(start, stop)
    days = period.get_segments("day")
    if args.exclude_month is not None:
        days.filter_month(args.exclude_month)

    # Processor Initialization
    # NOTE: This is only for later cases. Not much is done here at this
    #       point
    l2preproc = Level2PreProcessor(product_def)

#    # Loop over iterations (one per day)
    for day in days:

        # Do some extra logging
        logger.info("Processing Day [%s]" % day.label)

#        XXX: This needs a bit more thought
#        # Product Data Management
#        if args.remove_old:
#            for output_handler in product_def.output_handler:
#                output_handler.remove_old(day)

        # Get input files
        l2i_daily_files = l2i_handler.get_files_for_day(day.tcs.dt)
        if len(l2i_daily_files) == 0:
            logger.info("- no l2i products, skip day")
            continue
        logger.info("- Found %g l2i product files" % len(l2i_daily_files))

        # Process the orbits
        l2preproc.process_l2i_files(l2i_daily_files, day)

    # All done, log processor time
    t1 = time.process_time()
    seconds = int(t1-t0)
    logger.info("Run completed in %s" % str(timedelta(seconds=seconds)))


class Level2PreProcArgParser(DefaultLoggingClass):

    def __init__(self):
        super(Level2PreProcArgParser, self).__init__(self.__class__.__name__)
        self.error = ErrorStatus()
        self._args = None

    def parse_command_line_arguments(self):
        # use python module argparse to parse the command line arguments
        # (first validation of required options and data types)
        self._args = self.parser.parse_args()

    def critical_prompt_confirmation(self):

        # Any confirmation prompts can be overriden by --no-critical-prompt
        no_prompt = self._args.no_critical_prompt

        # if --remove_old is set, all previous l1bdata files will be
        # erased for all month
        if self._args.remove_old and not no_prompt:
            message = "You have selected to remove all previous " + \
                "l2p files for the requested period\n" + \
                "(Note: use --no-critical-prompt to skip confirmation)\n" + \
                "Enter \"YES\" to confirm and continue: "
            result = input(message)

            if result != "YES":
                sys.exit(1)

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

    @property
    def arg_dict(self):
        """ Return the arguments as dictionary """
        return self._args.__dict__

    @property
    def start(self):
        return self._args.start_date

    @property
    def stop(self):
        return self._args.stop_date

    @property
    def exclude_month(self):
        return self._args.exclude_month

    @property
    def doi(self):
        return self._args.doi

    @property
    def overwrite_protection(self):
        return self._args.overwrite_protection

    @property
    def l2i_product_dir(self):
        l2i_product_directories = self._args.l2i_product_dir
        if len(l2i_product_directories) != 1:
            raise ValueError("Must specify only 1 l2i directory")
        directory = Path(l2i_product_directories[0]).resolve()
        if not directory.is_dir():
            raise IOError(f"Not a valid l2i product directory: {directory}")
        return Path(l2i_product_directories[0])

    @property
    def l2p_output(self):
        l2p_output = self._args.l2p_output
        filename = psrlcfg.get_settings_file("output", "l2p", l2p_output)
        if filename is not None:
            return filename
        msg = "Invalid l2p outputdef filename or id: %s\n" % l2p_output
        msg = msg + " \nRecognized Level-2 output definitions ids:\n"
        l2p_output_ids = psrlcfg.get_setting_ids("output", "l2p")
        for l2p_output_id in l2p_output_ids:
            msg = f'{msg}    - {l2p_output_id}' + "\n"
        self.error.add_error("invalid-l2p-outputdef", msg)
        self.error.raise_on_error()

    @property
    def remove_old(self):
        return self._args.remove_old and not self._args.overwrite_protection


if __name__ == "__main__":
    pysiral_l2preproc()
