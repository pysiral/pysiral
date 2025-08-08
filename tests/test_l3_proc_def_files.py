# -*- coding: utf-8 -*-
"""
Testing all Level-2 processor definition files for compabilility with
Level2Processor conventions

@author: Stefan Hendricks
"""

import unittest

import yaml
from pysiral.core.legacy_packages import AttrDict
from loguru import logger

from pysiral import psrlcfg
from pysiral.core.config import get_yaml_config
from pysiral.l2proc.procsteps import Level2ProcessorStepOrder

logger.disable("pysiral")


class TestL3ProcDef(unittest.TestCase):

    def setUp(self):

        # Get a list of processor definition files in the code repository
        self.l2procdef_files = psrlcfg.get_settings_files("proc", "l3")

    def testYamlSyntaxOfDefinitionFiles(self):
        for filename in self.l2procdef_files:
            content = get_yaml_config(filename)
            self.assertIsInstance(content, AttrDict, msg=filename)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestL3ProcDef)
    unittest.TextTestRunner(verbosity=2).run(suite)
