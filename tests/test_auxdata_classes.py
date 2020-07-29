# -*- coding: utf-8 -*-
"""
Created on Mon Jul 06 17:57:33 2015

@author: Stefan
"""

import unittest
from loguru import logger
logger.disable("pysiral")

from pysiral.auxdata import get_all_auxdata_classes


class TestAuxdataClasses(unittest.TestCase):

    def setUp(self):
        # Get all auxiliary data classes
        self.auxdata_classes = get_all_auxdata_classes()

    def testAuxdataClassHasL2Method(self):
        """
        Test if auxdata class has the mandatory get_l2_track_vars method
        :return:
        """
        for module, name, class_instance in self.auxdata_classes:
            self.assertTrue(hasattr(class_instance, "get_l2_track_vars"))


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestAuxdataClasses)
    unittest.TextTestRunner(verbosity=2).run(suite)
