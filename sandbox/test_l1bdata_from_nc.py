# -*- coding: utf-8 -*-

import os
import glob

from pysiral.config import ConfigInfo
from pysiral.l1bdata import L1bdataNCFile


def test_l1bdata_from_nc():
    """ create an L1BData object from a CryoSat-2 l1bdata nc file """

    # Get Configuration Information
    config = ConfigInfo()

    # Get an L1B SAR file
    l1b_directory = config.local_machine.l1b_repository.cryosat2.l1bdata
    l1b_directory = os.path.join(l1b_directory, "north", "2015", "03")
    l1b_files = glob.glob(os.path.join(l1b_directory, "*.nc"))

    # Read the file
    l1b = L1bdataNCFile(l1b_files[0])
    l1b.parse()


if __name__ == "__main__":
    test_l1bdata_from_nc()