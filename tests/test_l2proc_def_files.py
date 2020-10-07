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
from pysiral.l2proc.procsteps import Level2ProcessorStepOrder
from loguru import logger
logger.disable("pysiral")


class TestL2ProcDef(unittest.TestCase):

    def setUp(self):

        # Get a list of processor definition files in the code repository
        l2proc_ids = psrlcfg.get_setting_ids("proc", "l2")
        self.l2procdef_files = [psrlcfg.get_settings_file("proc", "l2", l2proc_id) for l2proc_id in l2proc_ids]

    def testConfigFileRootTags(self):
        required_tags = ["metadata", "auxdata", "procsteps"]
        for l2procdef_file in self.l2procdef_files:
            with open(str(l2procdef_file)) as fh:
                content = AttrDict(yaml.safe_load(fh))
                for required_tag in required_tags:
                    msg = "Search for root tag `{}` in {}".format(required_tag, l2procdef_file)
                    self.assertTrue(required_tag in content, msg=msg)

    def testConfigMetadataTags(self):
        required_tags = ["label", "product_line", "record_type", "platform", "file_version_tag", "hemisphere"]
        for l2procdef_file in self.l2procdef_files:
            with open(str(l2procdef_file)) as fh:
                content = AttrDict(yaml.safe_load(fh))
                for required_tag in required_tags:
                    msg = "Search for root tag `{}` in {}".format(required_tag, l2procdef_file)
                    self.assertTrue(required_tag in content.metadata, msg=msg)

    def testConfigFileL2ProcStepContent(self):
        for l2procdef_file in self.l2procdef_files:
            with open(str(l2procdef_file)) as fh:
                content = AttrDict(yaml.safe_load(fh))
                # This configuration error will be caught in method testConfigFileRootTags
                # No need to repeat it here
                if "procsteps" not in content:
                    continue
                is_valid = True
                try:
                    procsteps = Level2ProcessorStepOrder(content.procsteps)
                    procsteps.validate()
                except SystemExit:
                    is_valid = False
                self.assertTrue(is_valid, "Validating procsteps definition in {}".format(l2procdef_file))


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestL2ProcDef)
    unittest.TextTestRunner(verbosity=2).run(suite)
