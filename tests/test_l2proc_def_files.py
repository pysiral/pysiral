# -*- coding: utf-8 -*-
"""
Testing all Level-2 processor definition files for compabilility with
Level2Processor conventions

@author: Stefan Hendricks
"""

import yaml
import unittest
from attrdict import AttrDict
from pysiral import psrlcfg
from loguru import logger
logger.disable("pysiral")


class TestL2ProcDef(unittest.TestCase):

    def setUp(self):

        # Get a list of processor definition files in the code repository
        l2proc_ids = psrlcfg.get_setting_ids("proc", "l2")
        self.l2procdef_files = [psrlcfg.get_settings_file("proc", "l2", l2proc_id) for l2proc_id in l2proc_ids]

    def testYAMLSyntax(self):
        required_tags = ["id", "version_tag", "hemisphere", "mission", "auxdata", "procsteps"]
        for l2procdef_file in self.l2procdef_files:
            with open(str(l2procdef_file)) as fh:
                content = AttrDict(yaml.safe_load(fh))
                for required_tag in required_tags:
                    msg = "Search for tag {} in {}".format(required_tag, l2procdef_file)
                    self.assertTrue(required_tag in content, msg=msg)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestL2ProcDef)
    unittest.TextTestRunner(verbosity=2).run(suite)
