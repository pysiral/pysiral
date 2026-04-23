# -*- coding: utf-8 -*-

"""

"""

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"

import numpy as np


def inverse_power(
        x: np.ndarray,
        decay: float,
        offset: float,
        scale: float = 1.0
) -> np.ndarray:
    return scale * np.power(x, -1.0*decay) + offset
