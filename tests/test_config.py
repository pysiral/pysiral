# -*- coding: utf-8 -*-
"""
Testing the pysiral configuration management

@author: Stefan
"""

import unittest
from attrdict import AttrDict

from pysiral import psrlcfg


class TestPackage(unittest.TestCase):

    def testMissionConfig(self):
        self.assertIsInstance(psrlcfg.mission_def.content, AttrDict, msg="")
        self.assertIsInstance(psrlcfg.mission_def.platform_ids, list)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestPackage)
    unittest.TextTestRunner(verbosity=2, descriptions=True).run(suite)
