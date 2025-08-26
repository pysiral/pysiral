# -*- coding: utf-8 -*-

"""

"""

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"

import argparse
from loguru import logger
from typing import Literal, Union, List

from pysiral.core.flags import PysiralProcessingLevels, ProductProcessingLevels


def info(
    target: str = None,
    search_string: str = None,
) -> None:
    """
    Main entry point for the info script.

    This function provides information about the specified target.

    :param target: The target to get information about.
    :param search_string: Optional search string to filter information.
    """

    target = "basic" if target is None else target.lower()

    match target:

        case "basic":
            pysiral_basic_info()

        case "l1p_settings":
            pysiral_settings_info(
                PysiralProcessingLevels.LEVEL1,
                "proc",
                search_string
            )

        case "l2_settings":
            pysiral_settings_info(
                PysiralProcessingLevels.LEVEL2,
                "proc",
                search_string
            )
        case "l2_output":
            pysiral_settings_info(
                PysiralProcessingLevels.LEVEL2,
                "output",
                search_string
            )

        case "l2p_output":
            pysiral_settings_info(
                ProductProcessingLevels.LEVEL2_PREPROCESSED,
                "output",
                search_string
            )

        case "l3_settings":
            pysiral_settings_info(
                PysiralProcessingLevels.LEVEL3,
                "proc",
                search_string
            )
        case "l3_output":
            pysiral_settings_info(
                PysiralProcessingLevels.LEVEL3,
                "output",
                search_string
            )

def pysiral_basic_info() -> None:
    """
    Prints basic information about pysiral.
    """
    print("pysiral - Python SAR Interferometry and Radar Altimetry Library")
    print("Version: 0.9.0")
    print("Author: Stefan Hendricks")


def pysiral_settings_info(
    processing_level: Union[PysiralProcessingLevels, ProductProcessingLevels],
    settings_type: Literal["proc", "output"],
    search_string: str = None
) -> None:
    """
    Prints information about pysiral settings.

    :param processing_level: The processing level (e.g., LEVEL1, LEVEL2).
    :param settings_type: The type of settings (e.g., "proc", "output").
    :param search_string: Optional search string to filter information.
    """
    breakpoint()


class InfoScriptArguments(object):

    def __init__(self) -> None:
        self.parser = self.get_argument_parser()

    def get(self, args_list: List[str] = None) -> "argparse.Namespace":
        return self.parser.parse_args() if args_list is None else self.parser.parse_args(args_list)

    @staticmethod
    def get_argument_parser() -> argparse.ArgumentParser:
        """
            Set up the command line argument parser for the Level-2 Processor.

            :return: The argument parser object.
            """

        # create the parser
        parser = argparse.ArgumentParser(
            prog="pysiral info",
            description="""
                        List pysiral package information and allows to 
                        query settings information.
                        """,
            epilog="For more information, see: https://pysiral.readthedocs.io",
            formatter_class=lambda prog: argparse.HelpFormatter(prog, width=96, indent_increment=4)  # noqa: E501
        )

        target_choices = [
            "basic",
            "l1p_settings",
            "l2_settings",
            "l2_output",
            "l2p_output",
            "l3_settings",
            "l3_output"
        ]
        parser.add_argument(
            "-t", "--target",
            type=str,
            choices=target_choices,
            metavar="<target>",
            default="basic",
            help="""
            Name of configuration settings to get information about.
            """+ f" Options are: {target_choices}."
        )

        return parser
