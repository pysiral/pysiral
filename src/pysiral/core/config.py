# -*- coding: utf-8 -*-
#
# Copyright © 2015 Stefan Hendricks
#
# Licensed under the terms of the GNU GENERAL PUBLIC LICENSE
#
# (see LICENSE for details)

"""
Purpose:
    Returns content of configuration and definition files

Created on Mon Jul 06 10:38:41 2015

@author: Stefan
"""

import yaml

from pathlib import Path
from typing import Dict, Union

from core.legacy_classes import AttrDict
from src.pysiral import psrlcfg


# TODO: Marked as obsolete -> flag_dict now in mission_def yaml.
class RadarModes(object):

    flag_dict = {"lrm": 0, "sar": 1, "sin": 2}

    def __init__(self):
        pass

    @classmethod
    def get_flag(cls, mode_name: str) -> Union[int, None]:
        try:
            return cls.flag_dict[mode_name]
        except KeyError:
            return None

    @classmethod
    def get_name(cls, flag: int) -> Union[str, None]:
        return next(
            (
                mode_name
                for mode_name, mode_flag in cls.flag_dict.items()
                if flag == mode_flag
            ),
            None,
        )

    def name(self, index: int) -> str:
        i = list(self.flag_dict.values()).index(index)
        return list(self.flag_dict.keys())[i]

    @property
    def num(self) -> int:
        return len(self.flag_dict.keys())


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


def get_yaml_config(filename: Union[str, Path], output: str = "attrdict") -> Union[Dict, AttrDict]:
    """
    Parses the contents of a configuration file in .yaml format
    and returns the content in various formats

    :param filename: The full file path of the config file
    :param output: Dict or AttrDict (depends on `output` keyword)
    :return:
    """
    with open(str(filename), 'r') as fileobj:
        content_dict = yaml.safe_load(fileobj)

    return AttrDict(content_dict) if output == "attrdict" else content_dict
