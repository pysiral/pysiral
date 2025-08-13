#!/usr/bin/env python
# NOTE: pysiral-l1preproc is a complete re-design of pysiral-l1bpreproc.py and will successively replace the
#       older version

import argparse

from loguru import logger
from pathlib import Path
from typing import List, Union, Optional

from dateperiods import DatePeriod


from pysiral import get_cls, set_psrl_cpu_count
from pysiral.scripts.parser_items import (
    ProcessingPeriod, ExcludeMonths, Hemisphere, PlatformID,
    L1PSettings, SourceDatasetID, MultiProcesssingNumCores, USeMultiProcesssing
)
from pysiral.l1preproc import (Level1POutputHandler, Level1PreProcJobDef, get_preproc)


def l1preproc(
    l1p_settings: Union[str, Path] = None,
    processing_period: DatePeriod = None,
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
    :param processing_period: The processing period to run the preprocessor for.
    :param exclude_month: List of months to exclude from processing (optional).
    :param hemisphere: Hemisphere to process data for, default is "global".
    :param platform: Platform identifier (optional).
    :param output_handler_cfg: Configuration for the output handler (optional).
    :param source_repo_id: Identifier for the source repository (optional).

    :return: None
    """

    job = Level1PreProcJobDef(
        l1p_settings,
        processing_period,
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
        self.parser = self.get_argument_parser()

    def get(self, args_list: List[str] = None) -> "argparse.Namespace":
        args = self.parser.parse_args() if args_list is None else self.parser.parse_args(args_list)
        if args.multiprocesssing_num_cores is not None:
            set_psrl_cpu_count(args.multiprocesssing_num_cores)
        return args

    @staticmethod
    def get_argument_parser() -> argparse.ArgumentParser:
        """
        Set up the command line argument parser for the Level-1 Pre-Processor.

        :return: The argument parser object.
        """

        # Level-1 pre-processor specific help text for the hemisphere argument
        hemisphere_l1p_help = """
        Target hemisphere for processing. Options are 'global', 'nh', or 'sh'. The 
        latitude limit of the hemisphere is defined in the Level-1 pre-processor settings file.
        If 'global' is selected, the processor will run for both hemispheres, but still within the
        latitude limits.
        """

        # List of command line option required for the Level-1 pre-processor
        arg_item_list = [
            # Positional arguments
            L1PSettings(),
            ProcessingPeriod(),
            # Optional arguments
            ExcludeMonths(),
            Hemisphere(help=hemisphere_l1p_help),
            PlatformID(),
            SourceDatasetID(),
            USeMultiProcesssing(),
            MultiProcesssingNumCores()
        ]

        # create the parser
        parser = argparse.ArgumentParser(
            prog="pysiral l1preproc",
            description="""
            This script is used to generate Level-1 files (l1p) from individual source files.
            the Level-1 pre-processor is a tool that harmonizes the input data from various sources
            and prepares it for further processing in the Level-2 processor. Processing steps include:
            Generating continuous trajectories over the polar oceans, pre-computing of waveform shape
            parameters and harmonization of the data format. 
            """,
            epilog="For more information, see: https://pysiral.readthedocs.io",
            formatter_class=lambda prog: argparse.HelpFormatter(prog, width=96, indent_increment=4)  # noqa: E501
        )
        for arg_item in arg_item_list:
            arg_flags, args_dict = arg_item.get()
            parser.add_argument(*arg_flags, **args_dict)

        return parser
