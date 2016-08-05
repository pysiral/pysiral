# -*- coding: utf-8 -*-

# from pysiral.config import (TimeRangeRequest, DefaultCommandLineArguments)
from pysiral.config import DefaultCommandLineArguments
# from pysiral.errorhandler import ErrorStatus
from pysiral.l1bpreproc import L1bPreProcJob
# from pysiral.logging import DefaultLoggingClass

import argparse
import sys


def pysiral_l1bpreproc():
    """ Main call for l1b preprocessing """

    # Collect job settings from pysiral configuration data and
    # command line arguments
    args = L1bPreProcArgParser()

    # Parse and validate the command line arguments
    args.parse_command_line_arguments()

    # Get confirmation for critical choices (if necessary)
    args.critical_prompt_confirmation()

    # Set up the pre-processing job
    jobdef = L1bPreProcJob()

    # configure the job definition
    jobdef.options.from_dict(args.arg_dict)

    # There are several options of how the time range can be requested
    # a) start and stop month (no days specified)
    # b) full start and stop date
    # In addition, full month can be excluded from the pre-processing
    # (e.g. summer Arctic)
    jobdef.process_requested_time_range()

    # Check if everything is ok
    jobdef.error.raise_on_error()

    # L1b data is grouped by month
    # -> if requested time range exceeds one month, the preprocessor
    #    will run in several iterations
    jobdef.generate_preprocessor_iterations()

    # Get the mission-specific preprocessor class
    preprocessor = jobdef.get_mission_preprocessor()

    # Loop over iterations (one per month)
    for time_range in jobdef.iterations:

        # Start the pre-processing
        job = preprocessor()
        job.log.info("Processing period: %s" % time_range.label)

        # Set the pre-processing parameter
        job.set_job_definition(jobdef)

        # Get the list of input files from local machine def
        version = jobdef.input_version
        job.get_input_files_local_machine_def(time_range, version=version)

        # Check error, e.g. problems with local_machine_def, ...
        job.error.raise_on_error()

        # List might be empty
        if job.has_empty_file_list:
            job.log.info(" Skipping period: %s" % time_range.label)
            continue

        # Empty output folder (if --remove_old is set)
        if jobdef.remove_old:
            job.remove_old_l1bdata()

        # Pre-process data for one month
        job.execute()


class L1bPreProcArgParser(object):

    def __init__(self):
        self.args = None

    def parse_command_line_arguments(self):
        # use python module argparse to parse the command line arguments
        # (first validation of required options and data types)
        self.args = self.parser.parse_args()

    def critical_prompt_confirmation(self):

        # Any confirmation prompts can be overriden by --no-critical-prompt
        no_prompt = self.args.no_critical_prompt

        # if --remove_old is set, all previous l1bdata files will be
        # erased for all month
        if self.args.remove_old and not no_prompt:
            message = "You have selected to remove all previous " + \
                "l1bdata files for the requested period\n" + \
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
            ("-mission", "mission", "mission_id", True),
            ("-start", "date", "start_date", True),
            ("-stop", "date", "stop_date", True),
            ("-hemisphere", "hemisphere", "hemisphere", False),
            ("-exclude-month", "exclude-month", "exclude_month", False),
            ("-input-version", "input-version", "input_version", False),
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
        return self.args.__dict__


if __name__ == "__main__":
    pysiral_l1bpreproc()
