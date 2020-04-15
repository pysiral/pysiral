# -*- coding: utf-8 -*-
"""
Created on Mon Jul 06 17:57:33 2015

@author: Stefan
"""

import unittest

from attrdict import AttrDict

from pysiral import USER_CONFIG_PATH
from pysiral.config import get_yaml_config


class TestDefinitionfiles(unittest.TestCase):

    def setUp(self):
        pass

    def testYamlSyntaxOfDefinitionFiles(self):
        def_files = ["mission_def.yaml", "auxdata_def.yaml"]
        for def_file in def_files:
            filename = USER_CONFIG_PATH / def_file
            content = get_yaml_config(filename)
            self.assertIsInstance(content, AttrDict, msg=def_file)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDefinitionfiles)
    unittest.TextTestRunner(verbosity=2).run(suite)
