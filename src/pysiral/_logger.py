# -*- coding: utf-8 -*-

"""
A small module to handle logging in pysiral. This module is intended to be the first module
called in pysiral.__init__.py, so that the logger is set up before any other modules are imported.
"""

import logging
import os
import sys
import warnings
from datetime import datetime
from typing import Dict

from dateutil.relativedelta import relativedelta
from loguru import logger


class InterceptHandler(logging.Handler):

    def emit(self, record):
        # Get corresponding Loguru level if it exists
        try:
            level = logger.level(record.levelname).name
        except ValueError:
            level = record.levelno

        # Find caller from where originated the logged message
        frame, depth = logging.currentframe(), 2
        while frame.f_code.co_filename == logging.__file__:
            frame = frame.f_back
            depth += 1

        logger.opt(depth=depth, exception=record.exc_info).log(level, record.getMessage())


def get_duration(elapsed_seconds: int, fmt_str: str = "%H:%M:%S") -> str:
    """
    Get a formatted string representing the duration from a given number of elapsed seconds.
    """
    datum = datetime(1900, 1, 1)
    duration = datum + relativedelta(seconds=elapsed_seconds)
    return duration.strftime(fmt_str)


def logger_format(record: Dict) -> str:
    """
    Formatter function for loguru (Normal mode

    :param record:

    :return: logger string
    """
    elapsed_timestr = get_duration(record["elapsed"].total_seconds())
    fmt_str = (
        '<green>{time:YYYY-MM-DD HH:mm:ss.SS} - {elapsed_timestr}</green> | '
        '<level>{level:<8}</level> | '
        '<level>{message} </level>\n'
    )
    return fmt_str.format(elapsed_timestr=elapsed_timestr, **record)


def logger_format_debug(record: Dict) -> str:
    elapsed_timestr = get_duration(record["elapsed"].total_seconds())
    if record["function"] == "<module>":
        record["function"] = ""
    else:
        record["function"] = "." + record["function"]
    code_line = "{name}{function}:L{line}".format(**record)
    fmt_str = (
        '<green>{time:YYYY-MM-DD HH:mm:ss.SS} - {elapsed_timestr}</green> | '
        '<level>{level: <8}</level> | '
        '<cyan>{code_line}</cyan> : '
        '<level>{message}</level>\n'
    )
    return fmt_str.format(elapsed_timestr=elapsed_timestr, code_line=code_line, **record)


def set_logger() -> None:
    """
    Set up the logger for pysiral.
    """
    logger.remove()
    debug_mode = "PYSIRAL_DEBUG_LOGGER" in os.environ
    log_level = "DEBUG" if debug_mode else "INFO"
    logger_format_func = logger_format_debug if debug_mode else logger_format
    logger.add(sys.stderr, format=logger_format_func, enqueue=True, level=log_level)


# Suppress warnings
warnings.filterwarnings("ignore")


set_logger()


