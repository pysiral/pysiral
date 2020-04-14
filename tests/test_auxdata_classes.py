# -*- coding: utf-8 -*-
"""
Created on Mon Jul 06 17:57:33 2015

@author: Stefan
"""

import unittest

from pysiral.auxdata import AuxdataBaseClass

class TestAuxdataClassesfiles(unittest.TestCase):

    def setUp(self):

        # Get all auxiliary data classes
        all_classes = AuxdataBaseClass.__inheritors__
        print(all_classes)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestAuxdataClassesfiles)
    unittest.TextTestRunner(verbosity=2).run(suite)
