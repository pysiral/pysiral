# -*- coding: utf-8 -*-

from pysiral.config import (ConfigInfo, DefaultCommandLineArguments,
                            TimeRangeRequest)
from pysiral.errorhandler import ErrorStatus
from pysiral.datahandler import DefaultL1bDataHandler
from pysiral.l2proc import Level2Processor, Level2ProductDefinition
from pysiral.logging import DefaultLoggingClass
from pysiral.path import file_basename

from datetime import timedelta
import argparse
import glob
import time
import sys
import os
import re


def pysiral_l2proc():

    # Collect job settings from pysiral configuration data and
    # command line arguments
    args = Level2ProcArgParser()

    # Parse and validate the command line arguments
    args.parse_command_line_arguments()

    # Get confirmation for critical choices (if necessary)
    args.critical_prompt_confirmation()

    # From here on there are two options
    # a. Time range given -> Get l1bdata input with datahandler
    # b. Predefined set of l1b input files
    # Splitting into functions for clarity
    if args.is_time_range_request:
        pysiral_l2proc_time_range_job(args)
    else:
        pysiral_l2proc_l1b_predef_job(args)


def pysiral_l2proc_time_range_job(args):
    """ This is a Level-2 Processor job for a given time range """

    # Get start time of processor run
    t0 = time.clock()

    # Get the product definition
    product_def = Level2ProductDefinition(args.run_tag, args.l2_settings_file)
    mission_id = product_def.l2def.mission.id
    hemisphere = product_def.l2def.hemisphere

    # Specifically add an output handler
    product_def.add_output_definition(args.l2_output, overwrite_protection=args.overwrite_protection)

    # Break down the time range in individual month
    start, stop = args.start, args.stop
    job = TimeRangeRequest(start, stop, exclude_month=args.exclude_month)
    job.clip_to_mission(mission_id)
    job.raise_if_empty()

    # Prepare DataHandler
    l1b_data_handler = DefaultL1bDataHandler(mission_id, hemisphere, version=args.l1b_version)

    # Processor Initialization
    l2proc = Level2Processor(product_def)

#    # Loop over iterations (one per month)
    for time_range in job.iterations:

        # Do some extra logging
        l2proc.log.info("Processing period: %s" % time_range.label)

        # Product Data Management
        if args.remove_old:
            for output_handler in product_def.output_handler:
                output_handler.remove_old(time_range)

        # Get input files
        l1b_files = l1b_data_handler.get_files_from_time_range(time_range)
        l2proc.log.info("Found %g files in %s" % (len(l1b_files), l1b_data_handler.last_directory))

        # Process the orbits
        l2proc.process_l1b_files(l1b_files)

    # All done
    t1 = time.clock()
    seconds = int(t1-t0)
    l2proc.log.info("Run completed in %s" % str(timedelta(seconds=seconds)))


def pysiral_l2proc_l1b_predef_job(args):
    """ A more simple Level-2 job with a predefined list of l1b data files """

    # Get start time of processor run
    t0 = time.clock()

    # Get the product definition
    product_def = Level2ProductDefinition(args.run_tag, args.l2_settings_file)

    # Specifically add an output handler
    product_def.add_output_definition(args.l2_output, overwrite_protection=args.overwrite_protection)

    # Processor Initialization
    l2proc = Level2Processor(product_def)
    l2proc.process_l1b_files(args.l1b_predef_files)

    # All done
    t1 = time.clock()
    seconds = int(t1-t0)
    l2proc.log.info("Run completed in %s" % str(timedelta(seconds=seconds)))


class Level2ProcArgParser(DefaultLoggingClass):

    def __init__(self):
        super(Level2ProcArgParser, self).__init__(self.__class__.__name__)
        self.error = ErrorStatus()
        self.pysiral_config = ConfigInfo()
        self._args = None

    def parse_command_line_arguments(self):
        # use python module argparse to parse the command line arguments
        # (first validation of required options and data types)
        self._args = self.parser.parse_args()

        # Add additional check to make sure either `l1b-files` or
        # `start ` and `stop` are set
        l1b_file_preset_is_set = self._args.l1b_files_preset is not None
        start_and_stop_is_set = self._args.start_date is not None and \
            self._args.stop_date is not None

        if l1b_file_preset_is_set and start_and_stop_is_set:
            self.parser.error("-start & -stop and -l1b-files are exclusive")

        if not l1b_file_preset_is_set and not start_and_stop_is_set:
            self.parser.error("either -start & -stop or -l1b-files required")

    def critical_prompt_confirmation(self):

        # Any confirmation prompts can be overridden by --no-critical-prompt
        no_prompt = self._args.no_critical_prompt

        # if --remove_old is set, all previous l1bdata files will be
        # erased for all month
        if self._args.remove_old and not no_prompt:
            message = "You have selected to remove all previous " + \
                "l2 files for the requested period\n" + \
                "(Note: use --no-critical-prompt to skip confirmation)\n" + \
                "Enter \"YES\" to confirm and continue: "
            result = raw_input(message)

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
            ("-l2-settings", "l2-settings", "l2_settings", True),
            ("-run-tag", "run-tag", "run_tag", False),
            ("-start", "date", "start_date", False),
            ("-stop", "date", "stop_date", False),
            ("-l1b-files", "l1b_files", "l1b_files_preset", False),
            ("-exclude-month", "exclude-month", "exclude_month", False),
            ("-input-version", "input-version", "input_version", False),
            ("-l2-output", "l2-output", "l2_output", False),
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

        parser.set_defaults(overwrite_protection=True)

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
    def run_tag(self):
        """ run_tag is a str or relative path that determines the output directory for
        the Level-2 processor. If the -run-tag option is not specified, the output
        directory will be the `product_repository` specification in `local_machine_def`
        with the l2 settings file basename as subfolder.

        One can however specify a custom string, or a relative path, with subfolders
        defined by `\` or `/`, e.g.

        Examples:
            -run-tag cs2awi_v2p0_nrt
            -run-tag c3s/cdr/cryosat2/v1p0/nh
        """

        # Get from command line arguments (default: None)
        run_tag = self._args.run_tag

        # If argument is empty use the basename of the l2 settings file
        if run_tag is None:
            run_tag = self._args.l2_settings
            # Settings file may be specified as full path and not just the id
            if os.path.isfile(run_tag):
                run_tag = file_basename(run_tag)

        # split the run-tag on potential path separators
        run_tag = re.split(r'[\\|/]', run_tag)
        return run_tag

    @property
    def exclude_month(self):
        return self._args.exclude_month

    @property
    def overwrite_protection(self):
        return self._args.overwrite_protection

    @property
    def l2_settings_file(self):
        l2_settings = self._args.l2_settings
        filename = self.pysiral_config.get_settings_file("proc","l2", l2_settings)
        if filename is None:
            msg = "Invalid l2 settings filename or id: %s\n" % l2_settings
            msg = msg + " \nRecognized Level-2 processor setting ids:\n"
            for l2_settings_id in self.pysiral_config.get_setting_ids("proc", "l2"):
                msg = msg + "  " + l2_settings_id+"\n"
            self.error.add_error("invalid-l2-settings", msg)
            self.error.raise_on_error()
        else:
            return filename

    @property
    def l1b_version(self):
        return self._args.input_version

    @property
    def l1b_predef_files(self):
        l1b_files = glob.glob(self._args.l1b_files_preset)
        return l1b_files

    @property
    def l2_output(self):
        l2_output = self._args.l2_output
        filename = self.pysiral_config.get_settings_file("output", "l2i", l2_output)
        if filename is None:
            msg = "Invalid l2 outputdef filename or id: %s\n" % l2_output
            msg = msg + " \nRecognized Level-2 output definitions ids:\n"
            l2_output_ids = self.pysiral_config.get_setting_ids("output", "l2i")
            for l2_output_id in l2_output_ids:
                msg = msg + "    - " + l2_output_id+"\n"
            self.error.add_error("invalid-l2-outputdef", msg)
            self.error.raise_on_error()
        else:
            return filename

    @property
    def is_time_range_request(self):
        return self._args.l1b_files_preset is None

    @property
    def remove_old(self):
        return self._args.remove_old and not self._args.overwrite_protection


if __name__ == "__main__":
    pysiral_l2proc()
