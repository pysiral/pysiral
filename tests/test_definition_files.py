# -*- coding: utf-8 -*-
"""
Created on Mon Jul 06 17:57:33 2015

@author: Stefan
"""

import datetime
import unittest

from attrdict import AttrDict

from pysiral import psrlcfg
from pysiral.config import get_yaml_config


class TestDefinitionfiles(unittest.TestCase):

    def setUp(self):
        pass

    def testYamlSyntaxOfDefinitionFiles(self):
        def_files = ["mission_def.yaml", "auxdata_def.yaml"]
        for def_file in def_files:
            filename = psrlcfg.config_path / def_file
            content = get_yaml_config(filename)
            self.assertIsInstance(content, AttrDict, msg=def_file)

    def testPlatformDefinitions(self):
        self.assertIsInstance(psrlcfg.platforms.ids, list)
        for platform_id in psrlcfg.platforms.ids:
            self.assertIsNotNone(psrlcfg.platforms.get_name(platform_id))
            self.assertIsNotNone(psrlcfg.platforms.get_sensor(platform_id))
            self.assertIsNotNone(psrlcfg.platforms.get_orbit_inclination(platform_id))
            tcs, tce = psrlcfg.platforms.get_time_coverage(platform_id)
            self.assertIsInstance(tcs, datetime.datetime)
            self.assertIsInstance(tce, datetime.datetime)

    def testAuxdataDefinitionBasic(self):
        self.assertTrue(hasattr(psrlcfg, "auxdef"))
        keys = psrlcfg.auxdef.iter_keys
        self.assertGreater(len(keys), 0)

    def testAuxdataDefinitionContent(self):
        """
        Test the content of the auxdata definitions (both in the userhome as well as the package,
        since they might be different)
        :return:
        """
        current_config_target = psrlcfg.current_config_target
        for config_target in psrlcfg.VALID_CONFIG_TARGETS:
            psrlcfg.set_config_target(config_target)
            psrlcfg.reload()
            self._testAuxdataDefinitionContent(config_target)
        psrlcfg.set_config_target(current_config_target)
        psrlcfg.reload()

    def _testAuxdataDefinitionContent(self, target):
        for category, id, item in psrlcfg.auxdef.items:
            aux_id = "{}:{}:{}".format(target, category, id)
            config_dict_keys = item.keys
            for required_key in ["options", "long_name", "pyclass", "local_repository"]:
                self.assertTrue(required_key in config_dict_keys, msg="{} has {}".format(aux_id, required_key))

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDefinitionfiles)
    unittest.TextTestRunner(verbosity=2).run(suite)
