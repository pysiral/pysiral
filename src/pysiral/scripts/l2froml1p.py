# -*- coding: utf-8 -*-

"""
"""


import argparse
import glob
import re
import sys
import time
from datetime import timedelta

from dateperiods import DatePeriod
from loguru import logger

from pysiral import psrlcfg, set_psrl_cpu_count
from pysiral.scripts.parser_items import DefaultCommandLineArguments
from pysiral.core.datahandler import L1PDataHandler
from pysiral.core.legacy_classes import DefaultLoggingClass, ErrorStatus
from pysiral.l2proc import Level2Processor, Level2ProductDefinition


__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"


def l2froml1p(
        run_tag: str,
        l2_settings_file: str,
        l1b_predef_files: list[str],
        l2_output: list[str],
        force_l2def_record_type: str = None,
        **_: dict
):
    """ A more simple Level-2 job with a predefined list of l1b data files """

    # Get start time of processor run
    t0 = time.time()

    # Get the product definition
    product_def = Level2ProductDefinition(
        run_tag,
        l2_settings_file,
        force_l2def_record_type=force_l2def_record_type
    )

    # Specifically add an output handler
    for l2_output in l2_output:
        product_def.add_output_definition(l2_output)

    # Processor Initialization
    l2proc = Level2Processor(product_def)
    l2proc.process_l1b_files(l1b_predef_files)

    # All done
    t1 = time.time()
    seconds = int(t1 - t0)
    logger.info(f"Run completed in {str(timedelta(seconds=seconds))}")
