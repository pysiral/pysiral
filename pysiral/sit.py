# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 17:15:39 2016

@author: shendric
"""

from pysiral.errorhandler import ErrorStatus
from pysiral.config import options_from_dictionary
import numpy as np


class L2ThicknessAlgorithmBaseClass(object):

    def __init__(self):
        self.error = ErrorStatus()

    def set_options(self, **opt_dict):
        self._options = options_from_dictionary(**opt_dict)

    def get_thickness(self, l2):
        sit, sit_unc, ice_dens, ice_dens_unc = self._get_thickness(l2)
        return sit, sit_unc, ice_dens, ice_dens_unc


class SeaIceFreeboardDefault(L2ThicknessAlgorithmBaseClass):
    """
    Classical sea ice freeboard to thickness conversion with
    Ricker et al. 2014 style uncertainty estimation
    """

    def __init__(self):
        super(SeaIceFreeboardDefault, self).__init__()

    def _get_thickness(self, l2):
        """
        Compute sea ice thickness based on a Level-2 data object
        - Assuming freeboard is sea ice freeboard
        - Sea ice density scaled based on FYI/MYI densities and MYI fraction
        """

        # Initial thickness value with nans
        sit = np.full((l2.n_records), np.nan, dtype=np.float32)
        sit_unc = np.full((l2.n_records), np.nan, dtype=np.float32)

        # Calculate thickness only for sea ice waveforms
        sea_ice = l2.surface_type.sea_ice.indices

        # Short cuts
        rho_w = self._options.water_density
        rho_s = l2.snow_dens

        # Get sea ice density and uncertainty
        rho_i, rho_i_unc = self._get_sea_ice_density(l2.sitype)

        # Classical sea ice freeboard to thickness concersion
        sit[sea_ice] = l2.frb[sea_ice] * (rho_w)/(rho_w - rho_i[sea_ice]) + \
            l2.snow_depth[sea_ice] * (rho_s[sea_ice])/(rho_w - rho_i[sea_ice])

        # Compute uncertainty
        sit_unc = self._get_thickness_uncertainty(l2, rho_i, rho_i_unc)

        return sit, sit_unc, rho_i, rho_i_unc

    def _get_sea_ice_density(self, myi_fraction):
        """ sea ice density scaled between fyi and myi by myi_fraction """

        # Boundary conditions (all fyi, all myi)
        rho_i_fyi = self._options.fyi_density
        rho_i_myi = self._options.myi_density

        # Scales with MYI fraction
        rho_i = rho_i_fyi - myi_fraction * (rho_i_fyi - rho_i_myi)
        rho_i_unc = np.full(rho_i.shape, 0.0, dtype=np.float32)

        return rho_i, rho_i_unc

    def _get_thickness_uncertainty(self, l2, rho_i, rho_i_unc):
        """
        Compute the sea ice thickness uncertainty based on the uncertainties
        from the input parameters.

        XXX: For now it is assumed that the uncertainty of the sea water
             density is negligible
        """

        rho_w = self._options.water_density

        # Compute the partial derivates for the error propagation
        deriv_fb = rho_w/(rho_w-rho_i)
        deriv_rho_s = l2.snow_depth/(rho_w-rho_i)
        deriv_rho_i = (l2.frb*rho_w + l2.snow_depth * l2.snow_dens) / \
            (rho_w-rho_i)**2.
        deriv_sd = (l2.snow_dens)/(rho_w-rho_i)

        # Error propagation for statistical and systematic errors
        sit_unc = np.sqrt((deriv_fb**2. * l2.frb.uncertainty**2) +
                          (deriv_rho_i**2 * rho_i_unc**2) +
                          (deriv_sd**2. * l2.snow_depth.uncertainty**2.) +
                          (deriv_rho_s**2 * l2.snow_dens.uncertainty**2.))


#        # Error propagation for statistical and systematic errors
#        sit_stat_unc = np.sqrt((deriv_fb**2. * l2.frb.uncertainty**2) +
#                               (deriv_rho_i**2 * rho_i_unc**2))
#
#        sit_syst_unc = np.sqrt((deriv_sd**2. * l2.snow_depth.bias**2.) +
#                               (deriv_rho_i**2. * rho_i.bias**2.) +
#                               (deriv_fb**2. * l2.frb.bias**2.) +
#                               (deriv_rho_s**2 * l2.snow_dens.bias**2.))
#
#        # Total uncertainty sum of both uncertainties
#        sit_unc = sit_stat_unc + sit_syst_unc

        return sit_unc


def snowfreeboard2thickness(frb, sd, rho_w, rho_i, rho_s):
    return rho_w/(rho_w-rho_i)*frb + (rho_s-rho_w)/(rho_w-rho_i)*sd


def get_sit_algorithm(name):
    pyclass = globals().get(name, None)
    if pyclass is not None:
        return pyclass()
    else:
        return pyclass
