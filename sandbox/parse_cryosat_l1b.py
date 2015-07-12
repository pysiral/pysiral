# -*- coding: utf-8 -*-
"""
Created on Tue Jul 07 13:53:17 2015

@author: Stefan
"""

from pysiral.config import ConfigInfo
from pysiral.cryosat2.cryosat2_l1b import CryoSatL1B

import glob
import os


def parse_cryosat_l1b():
    """ Sandbox script to parse a cryosat-l1b file """

    # Get Path information
    info = ConfigInfo()

    # Get an L1B SAR file
    l1b_directory = info.local_machine.local_l1b_repository.cryosat2.sar
    l1b_directory = os.path.join(l1b_directory, "2015", "04")
    l1b_files = glob.glob(os.path.join(l1b_directory, "*.DBL"))

    # Just select one
    l1b_file = l1b_files[0]

    # Parse the file into an L1B object
    l1b = CryoSatL1B()
    l1b.filename = l1b_file
    l1b.parse()

    # Tests
#    print l1b.mph
#    print l1b.sph
#    print l1b.dsd


if __name__ == "__main__":
    parse_cryosat_l1b()
