#!/usr/bin/env python
# NOTE: pysiral-l1preproc is a complete re-design of pysiral-l1bpreproc.py and will successively replace the
#       older version

import argparse
import sys

from loguru import logger
from pathlib import Path
from typing import List, Union, Optional

from pysiral import get_cls, set_psrl_cpu_count
from pysiral.core.config import DefaultCommandLineArguments
from pysiral.l1preproc import (Level1POutputHandler, Level1PreProcJobDef, get_preproc)


def l1preproc(
    l1p_settings: Union[str, Path],
    start_date: List[int],
    stop_date: List[int],
    exclude_month: Optional[List[int]] = None,
    hemisphere: str = "global",
    platform: str = None,
    output_handler_cfg: dict = None,
    source_repo_id: str = None,
    **_: Optional[dict]
) -> None:
    """
    Workflow script of the pysiral l1b preprocessor.

    :param l1p_settings: Level-1 preprocessor settings file or ID.
    :param start_date: Start date of the time coverage as a list [year, month, day].
    :param stop_date: End date of the time coverage as a list [year, month, day].
    :param exclude_month: List of months to exclude from processing (optional).
    :param hemisphere: Hemisphere to process data for, default is "global".
    :param platform: Platform identifier (optional).
    :param output_handler_cfg: Configuration for the output handler (optional).
    :param source_repo_id: Identifier for the source repository (optional).

    :return: None
    """

    job = Level1PreProcJobDef(
        l1p_settings,
        start_date,
        stop_date,
        exclude_month,
        hemisphere,
        platform,
        output_handler_cfg,
        source_repo_id
    )

    # Take the time
    job.stopwatch.start()

    # 1. Get the input handler
    input_handler_def = job.l1pprocdef.input_handler
    input_handler_cls, err = get_cls(input_handler_def.module_name, input_handler_def.class_name, relaxed=False)
    input_handler = input_handler_cls(input_handler_def.options)

    # 2. Get the adapter class that transfers
    adapter_def = job.l1pprocdef.input_adapter
    input_adapter_cls, err = get_cls(adapter_def.module_name, adapter_def.class_name, relaxed=False)
    input_adapter = input_adapter_cls(adapter_def.options)

    # 3. Get the output handler
    output_handler_def = job.l1pprocdef.output_handler
    output_handler = Level1POutputHandler(output_handler_def.options)
    output_handler.cfg.update(**job.output_handler_cfg)

    # 4. Get the pre-processor
    preproc_def = job.l1pprocdef.level1_preprocessor
    pre_processor = get_preproc(preproc_def.type, input_adapter, output_handler, preproc_def.options)

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
        pre_processor.process_input_files(file_list)

    # Report processing time
    job.stopwatch.stop()
    logger.info(f"Level-1 PreProcessor finished in {job.stopwatch.get_duration()}")


class L1PreProcScriptArguments(object):

    def __init__(self):
        """
        Command line parser class for the pysiral Level-1 Pre-Processor
        """
        self.parser = self.define_argument_parser()

    def get(self, args_list: List[str] = None) -> "argparse.Namespace":
        args = self.parser.parse_args() if args_list is None else self.parser.parse_args(args_list)
        if args.mp_cpu_count is not None:
            set_psrl_cpu_count(args.mp_cpu_count)
        return args

    @staticmethod
    def define_argument_parser() -> argparse.ArgumentParser:
        """
        Set up the command line argument parser for the Level-1 Pre-Processor.

        :return: The argument parser object.
        """

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
            ("-mp-cpu-count", "mp-cpu-count", "mp_cpu_count", False),
        ]

        # create the parser
        parser = argparse.ArgumentParser()
        for option in options:
            argname, argtype, destination, required = option
            argparse_dict = clargs.get_argparse_dict(argtype, destination, required)
            parser.add_argument(argname, **argparse_dict)
        parser.set_defaults(overwrite_protection=False)

        return parser
