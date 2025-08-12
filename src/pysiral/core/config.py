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

from pathlib import Path
from typing import Dict, Union

import yaml

from pysiral.core.legacy_classes import AttrDict


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
