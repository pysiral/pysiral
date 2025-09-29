# -*- coding: utf-8 -*-
"""
Testing the pysiral configuration management

@author: Stefan
"""

import unittest
from pathlib import Path

from pysiral.core.legacy_classes import AttrDict
from loguru import logger

from pysiral import psrlcfg

logger.disable("pysiral")


class TestConfig(unittest.TestCase):

    def setUp(self):
        pass

    def testVersionFiles(self):
        from pysiral import __version__, __git_version__, __git_branch__, __git_origin__
        self.assertIsInstance(__version__, str)
        self.assertIsInstance(__git_version__, str)
        self.assertIsInstance(__git_branch__, str)
        self.assertIsInstance(__git_origin__, str)

    def testMissionConfig(self):
        self.assertIsInstance(psrlcfg.platforms.content, AttrDict)
        self.assertIsInstance(psrlcfg.platforms.ids, list)

    def testConfigPath(self):
        self.assertTrue(Path(psrlcfg.package_config_path).is_dir())
        self.assertTrue(Path(psrlcfg.package_path).is_dir())
        self.assertTrue(Path(psrlcfg.config_path).is_dir())

    def testLocalMachineDefinition(self):
        """
        Test the local machine definition
        :return:
        """
        if psrlcfg.local_machine is not None:
            self.assertTrue(hasattr(psrlcfg, "local_machine_def_filepath"))
            lmd_filepath = psrlcfg.local_machine_def_filepath
            self.assertIsInstance(lmd_filepath, Path)
            self.assertTrue(lmd_filepath.is_file())

    def testL1ProcessorDefinitions(self):
        """
        Test the processor definitions for L1 Pre-processor
        :return:
        """

        # Loop over all processor levels
        for processor_level in psrlcfg.processor_levels:

            # Get a list of all ids
            proc_defs = psrlcfg.get_processor_definition_ids(processor_level)
            label = "procdef:{}".format(processor_level)

            # lists of ids must be a list and longer than 0
            self.assertIsInstance(proc_defs, list, msg="Type is not list: {} [{}]".format(type(proc_defs), label))
            self.assertGreater(len(proc_defs), 0, msg="No definitions found for {}".format(label))

            # Load all processor definitions (must return a valid AttrDict)
            for proc_def_id in proc_defs:
                filepath = psrlcfg.get_settings_file("proc", processor_level, proc_def_id)
                self.assertIsInstance(filepath, Path)
                self.assertTrue(filepath.is_file())


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestConfig)
    unittest.TextTestRunner(verbosity=2, descriptions=True).run(suite)
