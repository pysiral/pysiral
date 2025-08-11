# -*- coding: utf-8 -*-
"""
Testing all Level-2 output definition files for compabilility with
Level2Processor conventions

#TODO: Add CF/ACDD compliance checks

@author: Stefan Hendricks
"""

import unittest

from core.legacy_classes import AttrDict
from loguru import logger

from src.pysiral import psrlcfg
from core.config import get_yaml_config

logger.disable("pysiral")


class TestL2OutputDef(unittest.TestCase):

    def setUp(self):

        # Get a list of processor definition files in the code repository
        l2i_files = psrlcfg.get_settings_files("output", "l2i")
        l2p_files = psrlcfg.get_settings_files("output", "l2p")
        self.l2_output_files = l2i_files + l2p_files

    def testYamlSyntaxOfDefinitionFiles(self):
        for filename in self.l2_output_files:
            content = get_yaml_config(filename)
            self.assertIsInstance(content, AttrDict, msg=filename)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestL2OutputDef)
    unittest.TextTestRunner(verbosity=2).run(suite)
