# -*- coding: utf-8 -*-

import numpy as np

__all__ = ["SnowParameterContainer"]


class SnowParameterContainer(object):

    def __init__(self):
        self.depth = None
        self.density = None
        self.depth_uncertainty = None
        self.density_uncertainty = None

    def set_invalid(self, indices):
        self.depth[indices] = np.nan
        self.density[indices] = np.nan
        self.depth_uncertainty[indices] = np.nan
        self.density_uncertainty[indices] = np.nan

    def set_dummy(self, n_records):
        self.depth = np.full(n_records, np.nan)
        self.density = np.full(n_records, np.nan)
        self.depth_uncertainty = np.full(n_records, np.nan)
        self.density_uncertainty = np.full(n_records, np.nan)