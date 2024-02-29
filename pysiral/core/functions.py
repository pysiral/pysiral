# -*- coding: utf-8 -*-

"""

"""

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"

import numpy as np


def inverse_power(
        x: np.ndarray,
        scale: float,
        decay: float,
        offset: float
) -> np.ndarray:
    return scale * np.power(x, -1.0*decay) + offset
