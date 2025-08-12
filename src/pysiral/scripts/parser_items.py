# -*- coding: utf-8 -*-

"""
This module defines default command line arguments for pysiral scripts.
"""

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"

import argparse
from typing import List, Union, Callable, Any, Tuple, ClassVar
from dataclasses import dataclass, asdict, field


from pysiral import psrlcfg
from pysiral.core.flags import Hemispheres
from pysiral.scripts._argparse_types import file_type


@dataclass(frozen=True, kw_only=True)
class ArgparseArgumentsArgs:
    name_or_flags: ClassVar[list[str]] = []
    action: Union[str, Callable] = field(default=None)
    nargs: str = field(default=None)
    const: Union[str, int] = field(default=None)
    default: Any = field(default=None)
    type: Callable = field(default=None)
    choices: ClassVar[list[Any] | None] = None
    required: bool = field(default=False)
    help: str = field(default=None)
    metavar: str = field(default=None)
    dest: str = field(default=None)

    def get(self) -> Tuple[List[str], dict[str, Any]]:
        """
        Returns the input arguments and keyword arguments to `parser.add_argument()`

        :return: A tuple containing a list of argument names or flags and a
            ictionary of keyword arguments
        """
        # Convert the `args` property to a dictionary and remove empty values`
        args_dict = asdict(self)
        args_dict = {k: v for k, v in args_dict.items() if v is not None}

        # Remove the `name_or_flags` key from the dictionary and return it separately
        arg_flags = args_dict.pop("name_or_flags", [])
        return arg_flags, args_dict


@dataclass(frozen=True, kw_only=True)
class PlatformID(ArgparseArgumentsArgs):
    name_or_flags: ClassVar[list[str]] = ["-p", "--platform-id"]
    action: str = "store"
    dest: str = "platform_id"
    type: Callable = str
    required: bool = True
    choices: ClassVar[list[Any]] = psrlcfg.platform_ids
    help: str = "radar altimeter platform id (see choices)"


@dataclass(frozen=True, kw_only=True)
class StartDate(ArgparseArgumentsArgs):
    name_or_flags: ClassVar[list[str]] = ["start_date"]
    nargs: str = "+"
    type: Callable = int
    metavar: str = "YYYY MM [DD]"
    help: str = "Start date for processing period, given as year, month, and optionally day"


@dataclass(frozen=True, kw_only=True)
class EndDate(ArgparseArgumentsArgs):
    name_or_flags: ClassVar[list[str]] = ["end_date"]
    nargs: str = "+"
    type: Callable = int
    metavar: str = "YYYY MM [DD]"
    help: str = "End date for processing period, given as year, month, and optionally day"


@dataclass(frozen=True, kw_only=True)
class ExcludeMonths(ArgparseArgumentsArgs):
    name_or_flags: ClassVar[list[str]] = ["-E", "--exclude-months"]
    action: str = "store"
    dest: str = "exclude_months"
    metavar: str = "[month_num, month_num, ...]"
    type: Callable = int
    help: str = "List of months to be excluded from processing, given as integers (1-12)."


@dataclass(frozen=True, kw_only=True)
class InputVersion(ArgparseArgumentsArgs):
    name_or_flags: ClassVar[list[str]] = ["-v", "--input-version"]
    action: str = "store"
    dest: str = "input_version"
    default: Any = "default"
    metavar: str = "v{minor}p{major}|default"
    type: Callable = str
    help: str = """
    Input version name (e.g., v1p0, v2p1, or 'default'). 
    This is used to identify the version of the input data.
    """

@dataclass(frozen=True, kw_only=True)
class InputDataset(ArgparseArgumentsArgs):
    name_or_flags: ClassVar[list[str]] = ["-s", "--source-dataset-id"]
    action: str = "store"
    dest: str = "source_dataset_id"
    type: Callable = str
    help: str = """
    Identifier of the source dataset to be used for processing, summarizing the platform, version, 
    and timeliness information. The source dataset ID must be specified in the local machine definition file
    ({pysiral-cfg-location}/local_machine_def.yaml, e.g. `root.l1b_repository.<platform>.<source_dataset_id>`)
    """

@dataclass(frozen=True, kw_only=True)
class L1PFile(ArgparseArgumentsArgs):
    name_or_flags: ClassVar[list[str]] = ["l1p_file"]
    type: Callable = file_type(suffix=".nc")
    help: str = """
    Path to Level-1P (l1p) input file for the Level-2 processor.
    """



class DefaultCommandLineArguments(object):

    def __init__(self):

        self._args = {

            # Mission id
            "mission": {
                "action": 'store',
                "dest": 'mission_id',
                "choices": psrlcfg.platform_ids,
                "required": True,
                "help": "pysiral recognized mission id"
            },

            # platform (same as mission, but proper name)
            "platform": {
                "action": 'store',
                "dest": 'platform',
                "choices": psrlcfg.platform_ids,
                "required": True,
                "default": None,
                "help": "pysiral recognized platform id"
            },

            # Default date parameter
            "date": {
                "action": "store",
                "dest": "stop_date",
                "nargs": "+",
                "type": int,
                "required": False,
                "help": 'list as year and month and day (optional)'
            },

            # Default date parameter
            "hemisphere": {
                "action": "store",
                "dest": "hemisphere",
                "choices": ["global", "north", "south"],
                "default": "global",
                "required": False,
                "help": 'hemisphere flag for processing)'
            },

            # List of month to exclude from monthly iterations
            "exclude-month": {
                "action": "store",
                "dest": "exclude_month",
                "nargs": "+",
                "type": int,
                "required": False,
                "default": None,
                "help": 'list of months to be excluded from processing'
            },

            # Flag that indicates if previous versions shall be removed
            # before processing / plotting etc.
            "remove-old": {
                "action": "store_true",
                "dest": "remove_old",
                "default": False,
                "required": False,
                "help": 'remove all existing product in target directory'
            },

            # version tag of input data
            "input-version": {
                "action": "store",
                "dest": "input_version",
                "default": "default",
                "required": False,
                "help": 'input version name (see documentation)'
            },

            # same as input-version, but better worded
            "source-repo-id": {
                "action": "store",
                "dest": "source_repo_id",
                "default": None,
                "required": False,
                "help": 'specific tag in local_machine_def.yaml (root.l1b_repository.<platform>.<source_repo_od>'
            },

            # override any critical prompts for cronjobs etc
            "no-critical-prompt": {
                "action": "store_true",
                "dest": "no_critical_prompt",
                "default": False,
                "required": False,
                "help": 'set to skip any required command line inputs'
            },

            # preset for level-1b (l1bdata) fiels
            "l1b_files": {
                "action": "store",
                "dest": "l1b_files_preset",
                "default": None,
                "required": False,
                "help": 'Path to one or many l1bdata files (e.g.: path/*.nc)'
            },

            # fetch the level-1p file version
            "l1p-version": {
                "action": "store",
                "dest": "l1p-version",
                "default": None,
                "required": False,
                "help": 'file version of the l1p file'
            },

            # fetch the level-2 settings file
            "l1p-settings": {
                "action": "store",
                "dest": "l1p_settings",
                "default": None,
                "required": True,
                "help": 'id or path to Level-1P processor definition file file'
            },

            # fetch the level-2 settings file
            "l2-settings": {
                "action": "store",
                "dest": "l2_settings",
                "default": None,
                "required": True,
                "help": 'id or path to Level-2 settings file'
            },

            # fetch the level-2 settings file
            "l2-output": {
                "action": "store",
                "dest": "l2_output",
                "default": "l2i_default",
                "required": False,
                "help": 'l2 outputdef id'
            },

            # fetch the level-2 settings file
            "l2p-output": {
                "action": "store",
                "dest": "l2p_output",
                "default": "l2p_default",
                "required": False,
                "help": 'l2p outputdef id'
            },

            # set the run tag for the Level-2 Processor
            "run-tag": {
                "action": "store",
                "dest": "run_tag",
                "default": None,
                "required": False,
                "help": 'tag for the Level-2 output'
            },

            # no overwrite protection for level-2 outputs
            "no-overwrite-protection": {
                "action": "store_false",
                "dest": "overwrite_protection",
                "default": False,
                "required": False,
                "help": 'disable writing Level-2 output to unique directory'
            },

            # no overwrite protection for level-2 outputs
            "overwrite-protection": {
                "action": "store_true",
                "dest": "overwrite_protection",
                "default": False,
                "required": False,
                "help": 'enable writing Level-2 output to unique directory (default)'
            },

            "period": {
                "action": "store",
                "dest": "period",
                "default": "month",
                "required": False,
                "help": 'data period tag (default: month)'
            },

            "l2i-product-dir": {
                "action": "store",
                "dest": "l2i_product_dir",
                "nargs": "+",
                "default": None,
                "required": True,
                "help": "l2i input directory"
            },

            "l3-product-dir": {
                "action": "store",
                "dest": "l3_product_dir",
                "default": None,
                "required": False,
                "help": "l3 output directory"
            },

            "l3-settings": {
                "action": "store",
                "dest": "l3_settings",
                "default": "l3_default",
                "required": False,
                "help": "l3 settings definition id or filename"
            },

            "l3-griddef": {
                "action": "store",
                "dest": "l3_griddef",
                "default": None,
                "required": True,
                "help": "l3 grid definition id or filename"
            },

            "l3-output": {
                "action": "store",
                "dest": "l3_output",
                "default": "default",
                "required": True,
                "help": "l3 output id"
            },

            "doi": {
                "action": "store",
                "dest": "doi",
                "default": "None",
                "required": False,
                "type": str,
                "help": "doi number to be written in global attributes"
            },

            "data_record_type": {
                "action": "store",
                "dest": "data_record_type",
                "default": "None",
                "required": False,
                "type": str,
                "help": "type of data record [cdr, icdr]"
            },

            "force-l2def-record-type": {
                "action": "store_true",
                "dest": "force_l2def_record_type",
                "default": False,
                "required": False,
                "help": "overwrite l1p record type [cdr, icdr, nrt, ..] with metadata.record_type tag in l2def"
            },

            "mp-cpu-count": {
                "action": "store",
                "dest": "mp_cpu_count",
                "default": None,
                "type": int,
                "required": False,
                "help": "Number of CPU's to be used for multi-processing"
            },
        }

    def get_argparse_dict(self, name, destination, required):
        options = self._args[name]
        options["dest"] = destination
        options["required"] = required
        return options