# -*- coding: utf-8 -*-
"""
Created on Wed Sep 07 19:10:08 2016

@author: Stefan
"""

import timeit

setup = """
import numpy as np
a = np.arange(1000)
"""

print timeit.timeit('b=np.amax(a)', setup=setup)

setup2 = """
import numpy as np
import bottleneck as bn
a = np.arange(1000)
"""

print timeit.timeit('b=bn.nanmax(a)', setup=setup2)
