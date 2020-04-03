# -*- coding: utf-8 -*-
"""
Testing the pysiral configuration management

@author: Stefan
"""

import unittest
from attrdict import AttrDict

from pysiral import psrlcfg


class TestPackage(unittest.TestCase):

    def testConfigLoadOk(self):
        self.assertIsInstance(psrlcfg.mission, AttrDict)
        self.assertIsInstance(psrlcfg.auxdata, AttrDict)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestPackage)
    unittest.TextTestRunner(verbosity=2).run(suite)
