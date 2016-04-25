# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 17:15:39 2016

@author: shendric
"""

from pysiral.config import options_from_dictionary
import numpy as np


class L2ThicknessAlgorithmBaseClass(object):

    def __init__(self):
        pass

    def set_options(self, **opt_dict):
        self._options = options_from_dictionary(**opt_dict)

    def get_thickness(self, l2):
        thickness, msg = self._get_thickness(l2)
        return thickness, msg


class SeaIceFreeboardDefault(L2ThicknessAlgorithmBaseClass):
    """ Classical sea ice freeboard to thickness conversion """

    def __init__(self):
        super(SeaIceFreeboardDefault, self).__init__()

    def _get_thickness(self, l2):
        # Initial thickness value with nans
        sit = np.ndarray(shape=(l2.n_records))*np.nan
        # Calculate thickness only for sea ice waveforms
        sea_ice = l2.surface_type.sea_ice.indices
        # Short cuts
        rho_w = self._options.water_density
        rho_s = l2.snow_dens
        rho_i = self._get_rho_i(l2.sitype)
        # Classical sea ice freeboard to thickness c
        sit[sea_ice] = l2.frb[sea_ice] * (rho_w)/(rho_w - rho_i[sea_ice]) + \
            l2.snow_depth[sea_ice] * (rho_s[sea_ice])/(rho_w - rho_i[sea_ice])
        return sit, ""

    def _get_rho_i(self, myi_fraction):
        """ sea ice density scaled between fyi and myi by myi_fraction """
        rho_i_fyi = self._options.fyi_density
        rho_i_myi = self._options.myi_density
        return rho_i_fyi - myi_fraction * (rho_i_fyi - rho_i_myi)


def get_sit_algorithm(name):
    return globals()[name]()
