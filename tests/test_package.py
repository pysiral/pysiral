# -*- coding: utf-8 -*-
"""
@author: Stefan
"""

import unittest


class TestPackage(unittest.TestCase):

    def testAllImports(self):
        import pysiral


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestPackage)
    unittest.TextTestRunner(verbosity=2).run(suite)
