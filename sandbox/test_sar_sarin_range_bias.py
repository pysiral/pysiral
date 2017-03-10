# -*- coding: utf-8 -*-
"""
Created on Wed Feb 08 13:30:43 2017

@author: shendric
"""

from pysiral.cryosat2.functions import (get_cryosat2_wfm_range_userhandbook,
                                        get_cryosat2_wfm_range)

import numpy as np
import matplotlib.pyplot as plt

window_delay = 0.005 # sec
n_range_bins_sar = 512
n_range_bins_sin = 1024


# test sar
sar_range = get_cryosat2_wfm_range(window_delay, n_range_bins_sar)
sar_range_uh = get_cryosat2_wfm_range_userhandbook(window_delay, n_range_bins_sar)

print sar_range_uh[0]-sar_range[0]

sin_range = get_cryosat2_wfm_range(window_delay, n_range_bins_sin)
sin_range_uh = get_cryosat2_wfm_range_userhandbook(window_delay, n_range_bins_sin)

print sin_range_uh[0]-sin_range[0]


range_bin_index = np.arange(n_range_bins_sar)
plt.figure("sar")
plt.plot(range_bin_index, sar_range, lw=4, label="l1bdata")
plt.plot(range_bin_index, sar_range_uh, label="user handbook")



range_bin_index = np.arange(n_range_bins_sin)
plt.figure("sarin")
plt.plot(range_bin_index, sin_range, lw=4, label="l1bdata")
plt.plot(range_bin_index, sin_range_uh, label="user handbook")
plt.show()

