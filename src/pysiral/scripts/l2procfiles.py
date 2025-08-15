# -*- coding: utf-8 -*-

"""
"""

import argparse
from typing import List, Union
from pathlib import Path
from loguru import logger

from pysiral import set_psrl_cpu_count
from pysiral.l2proc import Level2Processor, Level2ProductDefinition

from pysiral.scripts.parser_items import (
    L2Settings, L2Outputs, MultiProcesssingNumCores,
    UseMultiProcesssing, ForceL2DefRecordType, L1PFiles,
    DOI
)


__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"


def l2procfiles(
        l1p_files: list[Path] = None,
        l2_settings: Union[Path, str] = None,
        l2_outputs: list[Path] = None,
        mp_cpu_count: int = None,
        force_l2def_record_type: bool = False,
        **_: dict
):
    """
    A more simple Level-2 job with a predefined list of l1b data files

    :param l1p_files: A list of Level-1P files to process.
    :param l2_settings: Path to the Level-2 settings file or its identifier.
    :param l2_outputs: List of output definitions for Level-2 products.
    :param mp_cpu_count: Number of CPUs to use for multiprocessing (optional).
    :param force_l2def_record_type: If True, forces the use of a specificLevel-2 definition record type.
    :param _:

    :return:
    """

    # Update pysiral multiprocessing settings
    if mp_cpu_count is not None:
        logger.info(f"Using {mp_cpu_count} CPU cores.")
        set_psrl_cpu_count(mp_cpu_count)

    # Get the product definition
    product_def = Level2ProductDefinition(
        l2_settings,
        force_l2def_record_type=force_l2def_record_type
    )

    # Specifically add an output handler
    for l2_output in l2_outputs:
        product_def.add_output_definition(l2_output)

    # Processor Initialization
    l2proc = Level2Processor(product_def)
    l2proc.process_l1b_files(l1p_files)


class L2ProcFilesScriptArguments(object):

    def __init__(self):
        self.parser = self.get_argument_parser()

    def get(self, args_list: List[str] = None) -> "argparse.Namespace":
        args = self.parser.parse_args() if args_list is None else self.parser.parse_args(args_list)
        if args.multiprocessing_num_cores is not None:
            set_psrl_cpu_count(args.multiprocessing_num_cores)
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
            L2Settings(),
            L1PFiles(),
            L2Outputs(required=True),
            # Optional arguments
            DOI(),
            UseMultiProcesssing(),
            MultiProcesssingNumCores(),
            ForceL2DefRecordType()
        ]

        # create the parser
        parser = argparse.ArgumentParser(
            prog="pysiral l2procfiles",
            description="""
                    The Level-2 Processor (l2proc) generates Level-2 files (l2/l2i) from l1p input files.
                    Level-2 files contain geophysical information and auxiliary data along the 
                    orbit at full sensor resolution. The processor uses a Level-2 product definition
                    file to define the product metadata, the list of auxiliary data files and the
                    algorithm steps to be applied to the Level-1P data. The output can be written
                    into multiple files.
                    Input Level-1 files are automatically selected based on the the source dataset ID 
                    and l1p version. ((see also: `pysiral l2proc --help` for running the 
                    Level-2 processor for a time range).
                    """,
            epilog="""
            For more information, see: https://pysiral.readthedocs.io
            """,
            formatter_class=lambda prog: argparse.HelpFormatter(prog, width=96, indent_increment=4)  # noqa: E501
        )
        for arg_item in arg_item_list:
            arg_flags, args_dict = arg_item.get()
            parser.add_argument(*arg_flags, **args_dict)

        return parser
