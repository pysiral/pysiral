# -*- coding: utf-8 -*-
"""
Created on Mon Jul 06 17:57:33 2015

@author: Stefan
"""

import unittest

from treedict import TreeDict
import os

from pysiral import USER_CONFIG_PATH
from pysiral.config import get_yaml_config


class TestDefinitionfiles(unittest.TestCase):

    def setUp(self):
        pass

    def testYamlSyntaxOfDefinitionFiles(self):
        def_files = ["mission_def.yaml", "auxdata_def.yaml"]
        for def_file in def_files:
            filename = os.path.join(USER_CONFIG_PATH, def_file)
            self.assertIsInstance(get_yaml_config(filename), TreeDict, msg=def_file)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDefinitionfiles)
    unittest.TextTestRunner(verbosity=2).run(suite)
