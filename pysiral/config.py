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

import os
import yaml
from treedict import TreeDict


class ConfigInfo(object):
    """
    Container for the content of the pysiral definition files
    (in pysiral/configration) and the local machine definition file
    (local_machine_definition.yaml)
    """

    # Global variables
    _DEFINITION_FILES = {
        "mission": os.path.join("config", "mission_def.yaml"),
        "area": os.path.join("config", "area_def.yaml"),
        "auxdata": os.path.join("config", "auxdata_def.yaml"),
        "products": os.path.join("config", "product_def.yaml"),
        "parameter": os.path.join("config", "parameter_def.yaml"),
        "local_machine": "local_machine_def.yaml"
    }

    def __init__(self):
        """ Read all definition files """
        self.pysiral_local_path = get_pysiral_local_path()
        for key in self._DEFINITION_FILES.keys():
            content = get_yaml_config(
                os.path.join(
                    self.pysiral_local_path,
                    self._DEFINITION_FILES[key]))
            setattr(self, key, content)

    @property
    def doc_path(self):
        return self._return_path("doc")

    def _return_path(self, subfolder):
        return os.path.join(self.pysiral_local_path, subfolder)


def get_yaml_config(filename, output="treedict"):
    """
    Parses the contents of a configuration file in .yaml format
    and returns the content in various formats

    Arguments:
        filename (str)
            path the configuration file

    Keywords:
        output (str)
            "treedict" (default): Returns a treedict object
            "dict": Returns a python dictionary
    """
    with open(filename, 'r') as f:
        content_dict = yaml.load(f)

    if output == "treedict":
        return TreeDict.fromdict(content_dict, expand_nested=True)
    else:
        return content_dict
