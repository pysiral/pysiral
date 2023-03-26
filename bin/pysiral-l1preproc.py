#!/usr/bin/env python
# NOTE: pysiral-l1preproc is a complete re-design of pysiral-l1bpreproc.py and will successively replace the
#       older version

import argparse
import sys

from loguru import logger

from pysiral import get_cls
from pysiral.config import DefaultCommandLineArguments
from pysiral.l1preproc import (Level1POutputHandler, Level1PreProcJobDef,
                               get_preproc)


def pysiral_l1preproc(job):
    """
    Workflow script of the pysiral l1b preprocessor.

    :param job: A pysiral.l1preproc.Level1PreProcJobDef instance
    :return: None
    """

    # Take the time
    job.stopwatch.start()

    # 1. Get the input handler
    input_handler_def = job.l1pprocdef.input_handler
    input_handler_cls = get_cls(input_handler_def.module_name, input_handler_def.class_name, relaxed=False)
    input_handler = input_handler_cls(input_handler_def.options)

    # 2. Get the adapter class that transfers
    adapter_def = job.l1pprocdef.input_adapter
    input_adapter_cls = get_cls(adapter_def.module_name, adapter_def.class_name, relaxed=False)
    input_adapter = input_adapter_cls(adapter_def.options)

    # 3. Get the output handler
    output_handler_def = job.l1pprocdef.output_handler
    output_handler = Level1POutputHandler(output_handler_def.options)
    output_handler.cfg.update(**job.output_handler_cfg)

    # 4. Get the pre-processor
    preproc_def = job.l1pprocdef.level1_preprocessor
    l1preproc = get_preproc(preproc_def.type, input_adapter, output_handler, preproc_def.options)

    # 5. Loop over monthly periods
    for period in job.period_segments:

        # 5.1 Get input files
        file_list = input_handler.get_file_for_period(period)
        if len(file_list) == 0:
            logger.warning(f"No input files found for period: {period.date_label}, skipping")

        # 5.2 Output management
        # Note: This is only relevant, if the --remove-old keyword is set
        output_handler.remove_old_if_applicable(period)

        # 5.3 Run the pre-processor
        l1preproc.process_input_files(file_list)

    # Report processing time
    job.stopwatch.stop()
    logger.info(f"Level-1 PreProcessor finished in {job.stopwatch.get_duration()}")


class Level1PreProcArgParser(object):

    def __init__(self):
        """
        Command line parser class for the pysiral Level-1 Pre-Processor
        """
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
                "l1p files for the requested period\n" + \
                "(Note: use --no-critical-prompt to skip confirmation)\n" + \
                "Enter \"YES\" to confirm and continue: "
            result = input(message)

            if result != "YES":
                sys.exit(1)

    @property
    def parser(self):

        # Take the command line options from default settings
        # -> see config module for data types, destination variables, etc.
        clargs = DefaultCommandLineArguments()

        # List of command line option required for pre-processor
        # (argname, argtype (see config module), destination, required flag)
        options = [
            ("-l1p-settings", "l1p-settings", "l1p_settings", True),
            ("-platform", "platform", "platform", False),
            ("-source-repo-id", "source-repo-id", "source_repo_id", False),
            ("-start", "date", "start_date", True),
            ("-stop", "date", "stop_date", True),
            ("-exclude-month", "exclude-month", "exclude_month", False),
            ("-hemisphere", "hemisphere", "hemisphere", False),
            ("--remove-old", "remove-old", "remove_old", False),
            ("--no-critical-prompt", "no-critical-prompt", "no_critical_prompt", False),
            ("--no-overwrite-protection", "no-overwrite-protection", "overwrite_protection", False),
            ("--overwrite-protection", "overwrite-protection", "overwrite_protection", False)]

        # create the parser
        parser = argparse.ArgumentParser()
        for option in options:
            argname, argtype, destination, required = option
            argparse_dict = clargs.get_argparse_dict(argtype, destination, required)
            parser.add_argument(argname, **argparse_dict)
        parser.set_defaults(overwrite_protection=False)

        return parser

    @property
    def arg_dict(self):
        """ Return the arguments as dictionary """
        return self._args.__dict__

    @property
    def args(self):
        return self._args


if __name__ == "__main__":

    # Get the command line arguments
    cmd_args = Level1PreProcArgParser()
    cmd_args.parse_command_line_arguments()

    # Create the job definitions
    job = Level1PreProcJobDef.from_args(cmd_args.args)

    # Execute Level-1 Pre-Processor Workflow
    pysiral_l1preproc(job)
