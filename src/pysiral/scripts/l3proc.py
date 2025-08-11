#!/usr/bin/env python
# -*- coding: utf-8 -*-


import argparse
import sys
import time
from datetime import timedelta
from pathlib import Path
from typing import List

from dateperiods import DatePeriod
from loguru import logger

from pysiral import psrlcfg
from pysiral.core import DefaultLoggingClass
from pysiral.core.config import DefaultCommandLineArguments
from pysiral.core.datahandler import L2iDataHandler
from pysiral.core.legacy_classes import ErrorStatus
from pysiral.l3proc import (Level3GridDefinition, Level3OutputHandler,
                            Level3Processor, Level3ProductDefinition)


def pysiral_l3proc():
    # parse command line arguments
    args = Level3ProcArgParser()
    args.parse_command_line_arguments()

    # Get start time of processor run
    t0 = time.process_time()

    # --- Get the period segments for the Level-3 processor ---
    # NOTE: These depend on the chosen total time range and the duration period for the grid.
    period = DatePeriod(args.start, args.stop)
    if args.period == "custom":
        period_segments = [period]
        n_periods = 1
    else:
        period_segments = period.get_segments(args.period)
        n_periods = period_segments.n_periods

    # Get the output grid
    grid = Level3GridDefinition(args.l3_griddef)

    # Initialize the interface to the l2i products
    l2i_handler = L2iDataHandler(args.l2i_product_directories, search_str="l2")

    # Initialize the output handler
    # Currently,  overwrite protection is disabled per default
    output = []
    for l3_output_file in args.l3_output_file:
        output_handler = Level3OutputHandler(output_def=l3_output_file,
                                             base_directory=args.l3_product_basedir,
                                             period=args.period,
                                             doi=args.doi,
                                             data_record_type=args.data_record_type,
                                             overwrite_protection=False)
        output.append(output_handler)

    # Compile the product def
    product_def = Level3ProductDefinition(args.l3_settings_file, grid, output, period)

    # Initialize the Processor
    l3proc = Level3Processor(product_def)

    # Loop over all iterations
    for i, time_range in enumerate(period_segments):

        # Report processing period
        msg = "# Processing %s period (%g of %g): %s"
        msg %= (args.period, i+1, n_periods, time_range.date_label)
        logger.info(msg)

        # Retrieve files
        l2i_files = l2i_handler.get_files_from_time_range(time_range)
        logger.info("Num l2i files: %g" % len(l2i_files))
        if len(l2i_files) == 0:
            logger.info("Skip data period")
            continue

        # Start the Level-3 processing
        l3proc.process_l2i_files(l2i_files, time_range)

    # Final reporting
    t1 = time.process_time()
    seconds = int(t1 - t0)
    logger.info("Run completed in %s" % str(timedelta(seconds=seconds)))


class Level3ProcArgParser(DefaultLoggingClass):

    def __init__(self):
        super(Level3ProcArgParser, self).__init__(self.__class__.__name__)
        self.error = ErrorStatus()
        self._args = None

    def parse_command_line_arguments(self):
        # use python module argparse to parse the command line arguments
        # (first validation of required options and data types)
        self._args = self.parser.parse_args()

        # Add addtional check to make sure either `l1b-files` or
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

    def critical_prompt_confirmation(self):

        # Any confirmation prompts can be overriden by --no-critical-prompt
        no_prompt = self._args.no_critical_prompt

        # if --remove_old is set, all previous l1bdata files will be
        # erased for all month
        if self._args.remove_old and not no_prompt:
            message = "You have selected to remove all previous " + \
                      "l3 files for the requested period\n" + \
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
    def period(self):
        return self._args.period

    @property
    def doi(self):
        return self._args.doi

    @property
    def data_record_type(self):
        return self._args.data_record_type

    @property
    def l2i_product_directories(self) -> List[Path]:
        return [Path(base_dir) for base_dir in self._args.l2i_basedir]

    @property
    def l3_settings_file(self):
        l3_settings = self._args.l3_settings
        filename = psrlcfg.get_settings_file("proc", "l3", l3_settings)
        if filename is not None:
            return filename
        msg = "Invalid l3 settings filename or id: %s\n" % l3_settings
        msg = msg + " \nRecognized Level-3 processor setting ids:\n"
        for l3_settings_id in psrlcfg.get_setting_ids("proc", "l3"):
            msg = f"{msg}  {l3_settings_id}" + "\n"
        self.error.add_error("invalid-l3-settings", msg)
        self.error.raise_on_error()

    @property
    def l3_griddef(self):
        l3_griddef = self._args.l3_griddef
        filename = psrlcfg.get_settings_file("grid", None, l3_griddef)
        if filename is not None:
            return filename
        msg = "Invalid griddef filename or id: %s\n" % l3_griddef
        msg = msg + "    Recognized grid definition ids:\n"
        for griddef_id in psrlcfg.get_setting_ids("grid"):
            msg = f"{msg}    - {griddef_id}" + "\n"
        self.error.add_error("invalid-griddef", msg)
        self.error.raise_on_error()

    @property
    def l3_output_file(self):
        """
        Get the full output definition file path. Multiple output filenames are possible if
        the command line argument is semicolon-separated list.
        :return: A list of output definition filenames
        """

        filenames = []
        for l3_output in self._args.l3_output.split(";"):
            filename = psrlcfg.get_settings_file("output", "l3", l3_output)
            if filename is None:
                msg = "Invalid output definition filename or id: %s\n" % l3_output
                msg = msg + "    Recognized output definition ids:\n"
                for output_id in psrlcfg.get_setting_ids("output", "l3"):
                    msg = f"{msg}    - {output_id}" + "\n"
                self.error.add_error("invalid-outputdef", msg)
            else:
                filenames.append(filename)

        if not filenames:
            msg = "No valid output definition file found for argument: %s"
            msg %= (str(self._args.l3_output))
            self.error.add_error("invalid-outputdef", msg)
            self.error.raise_on_error()

        return filenames

    @property
    def l3_product_basedir(self) -> Path:
        """
        Returns the base directory for the L3 output. There are several cases depending on
        whether multiple l2i product dirs are specified and most importantly if the L3 output
        directory has been specified (`--l3-product-dir`)

        L3C (single L2 data source):

        - L3C target dir: ../l2i/../l3c

        L3S (multiple L2 data sources):

        - Multiple l2i dirs, no dedicated l3 target dir -> ../l2i/../l3s of first l2i directory [Warning displayed]
        - dedicated l3 target dir -> used if speficied

        :return:
        """

        # The specification of `--l3-product-dir` trumps all other conditions
        if self._args.l3_product_dir:
            output_dir = Path(self._args.l3_product_dir)
            if not output_dir.is_dir():
                output_dir.mkdir(parents=True, exist_ok=True)
            return output_dir

        # 1. Clean up the path
        if not self.l2i_product_directories[0].resolve().is_dir():
            raise IOError(f"Not a valid l2i product directory: {self.l2i_product_directories[0]}")
        product_basedir = self.l2i_product_directories[0]
        dirs = product_basedir.parts
        l3_product_basedir = Path(*dirs[:-1]) if dirs[-1] in ["l2i", "l2"] else product_basedir
        if len(self.l2i_product_directories) > 1:
            logger.warning(
                f"Multiple l2i directories, but no dedicated l3 product directory. Defaulting to {l3_product_basedir}"
            )
        return l3_product_basedir

    @property
    def remove_old(self):
        return self._args.remove_old and not self._args.overwrite_protection


if __name__ == "__main__":
    pysiral_l3proc()
