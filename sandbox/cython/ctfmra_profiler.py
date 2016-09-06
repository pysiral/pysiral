# -*- coding: utf-8 -*-
"""
Created on Tue Sep 06 21:28:11 2016

@author: Stefan
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 14:06:28 2016

@author: shendric

Purpose: Test a faster implementation of the TFMRA retracker

"""

from pysiral.config import ConfigInfo
from pysiral.l1bdata import L1bConstructor
from pysiral.retracker import cTFMRA

import numpy as np

import cPickle as pickle
import os

TEST_FILE = r"CS_OFFL_SIR_SAR_1B_20160401T013146_20160401T013933_C001.DBL"


def ctfrma_profiler():

#    pysiral_config = ConfigInfo()
#
#    # Get the full path to a test data file
#    cs2_filename = os.path.join(pysiral_config.sandbox_path, "testdata",
#                                "cryosat2", TEST_FILE)
#
#    # Create l1b data object
#    l1b = L1bConstructor(pysiral_config)
#    l1b.mission = "cryosat2"
#    l1b.filename = cs2_filename
#    l1b.construct()
#
#
#    # Get reference data set
#    wfm = l1b.waveform.power
#    rng = l1b.waveform.range
#    rmode = l1b.waveform.radar_mode
#    is_valid = l1b.waveform.is_valid
#    indices = l1b.surface_type.get_by_name("ocean").indices
#
#
#    data = [wfm, rng, rmode, is_valid, indices, l1b.n_records]
#    with open("profile_test_data.pkl", "wb") as f:
#        pickle.dump(data, f)

    with open("profile_test_data.pkl", "rb") as f:
        wfm, rng, rmode, is_valid, indices, n_records =  pickle.load(f)

    ctfmra = cTFMRA()
    ctfmra.set_default_options()
    ctfmra.init(n_records)
    ctfmra.l2_retrack(rng, wfm, indices, rmode, is_valid)


if __name__ == "__main__":
    ctfrma_profiler()

