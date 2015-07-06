# -*- coding: utf-8 -*-
"""
Created on Mon Jul 06 17:57:33 2015

@author: Stefan
"""

import unittest

from treedict import TreeDict
import os
import numpy as np

from pysiral.config import get_pysiral_local_path, get_yaml_config


class TestDefinitionfiles(unittest.TestCase):

    def setUp(self):
        self.main_path = get_pysiral_local_path()

    def testYamlSyntaxOfDefinitionFiles(self):
        def_files = ["mission_def.yaml",
                     "parameter_def.yaml",
                     "area_def.yaml",
                     "auxdata_def.yaml",
                     "product_def.yaml"]
        for def_file in def_files:
            filename = os.path.join(self.main_path, "config", def_file)
            self.assertIsInstance(
                get_yaml_config(filename), TreeDict,
                msg=def_file)

    def testParameterIDsAreUnique(self):
        pars = get_parameters()
        ids = [par.short_name for par in pars.parameter.iterbranches()]
        ids_unique = list(np.unique(ids))
        self.assertListEqual(sorted(ids), sorted(ids_unique))

    def testParameterDefinitionsAreComplete(self):
        pars = get_parameters()
        tags = ["short_name", "level", "unit",
                "docstr", "dtype", "ascii_format"]
        for par in pars.parameter.iterbranches():
            for tag in tags:
                msg = tag+" in "+par.branchName()
                self.assertTrue(type(par[tag]) is str, msg=msg)


def get_parameters():
    """ Returns a Treedict of Parameter Definitions """
    config = get_yaml_config(
        os.path.join(get_pysiral_local_path(), "config", "parameter_def.yaml"))
    return config


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDefinitionfiles)
    unittest.TextTestRunner(verbosity=2).run(suite)
