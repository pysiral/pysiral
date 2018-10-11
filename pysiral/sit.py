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

        # Calculate thickness only for sea ice waveforms
        sea_ice = l2.surface_type.sea_ice.indices

        # Short cuts
        rho_w = self._options.water_density
        rho_s = l2.sdens

        # Get sea ice density and uncertainty
        rho_i, rho_i_unc = self._get_sea_ice_density(l2.sitype)

        # Classical sea ice freeboard to thickness conversion
        sit[sea_ice] = l2.frb[sea_ice] * (rho_w)/(rho_w - rho_i[sea_ice]) + \
            l2.sd[sea_ice] * (rho_s[sea_ice])/(rho_w - rho_i[sea_ice])

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

        # Compute uncertainty of sea ice density
        # Note: it has been decided to use a simplified computation of
        # sea ice density uncertainty, since a simple Gaussian error
        # propagation does not cut it, as the the two (fyi & myi) reference
        # uncertainties are not independent (linked via scaling factor).
        # One could add the covariance matrix, but this is hardly justified
        # given the already very basic assumptions for fyi and myi density
        # uncertainty.
        # Here, the uncertainty is based on the fyi and myi density
        # uncertainties scaled with the myi fraction. To reflect the
        # myi_fraction uncertainty, the sea ice density uncertainty is
        # slightly increased in regions of higher myi_fraction uncertainty
        if "uncertainty" in self._options:
            unc = self._options.uncertainty
            rho_i_unc = unc.fyi_density
            rho_i_unc += myi_fraction * (unc.myi_density - unc.fyi_density)
            rho_i_unc += (unc.fyi_density - unc.myi_density) * (myi_fraction.uncertainty)

        # Use fixed value for uncertainty (or 0 if unspecified)
        else:
            unc = 0.0
            if "uncertainty_fixed_value" in self._options:
                unc = self._options.uncertainty_fixed_value
            rho_i_unc = np.full(rho_i.shape, unc, dtype=np.float32)

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
        deriv_rho_s = l2.sd/(rho_w-rho_i)
        deriv_rho_i = (l2.frb*rho_w + l2.sd * l2.sdens) / (rho_w-rho_i)**2.
        deriv_sd = (l2.sdens)/(rho_w-rho_i)

        # Error propagation for statistical and systematic errors
        sit_unc = np.sqrt((deriv_fb**2. * l2.frb.uncertainty**2) +
                          (deriv_rho_i**2 * rho_i_unc**2) +
                          (deriv_sd**2. * l2.sd.uncertainty**2.) +
                          (deriv_rho_s**2 * l2.sdens.uncertainty**2.))

        return sit_unc


class SnowFreeboardDefault(L2ThicknessAlgorithmBaseClass):
    """
    Classical sea ice freeboard to thickness conversion with
    Ricker et al. 2014 style uncertainty estimation
    """

    def __init__(self):
        super(SnowFreeboardDefault, self).__init__()

    def _get_thickness(self, l2):
        """
        Compute sea ice thickness based on a Level-2 data object
        - Assuming freeboard is sea ice freeboard
        - Sea ice density scaled based on FYI/MYI densities and MYI fraction
        """

        # Initial thickness value with nans
        sit = np.full((l2.n_records), np.nan, dtype=np.float32)

        # Calculate thickness only for sea ice waveforms
        sea_ice = l2.surface_type.sea_ice.indices

        # Short cuts
        rho_w = self._options.water_density
        rho_s = l2.sdens

        # Get sea ice density and uncertainty
        rho_i, rho_i_unc = self._get_sea_ice_density(l2.sitype)

        # Classical sea ice freeboard to thickness conversion
        sit[sea_ice] = l2.frb[sea_ice] * (rho_w)/(rho_w - rho_i[sea_ice]) + \
            l2.sd[sea_ice] * (rho_s[sea_ice])/(rho_w - rho_i[sea_ice])

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

        # Compute uncertainty of sea ice density
        # Note: it has been decided to use a simplified computation of
        # sea ice density uncertainty, since a simple Gaussian error
        # propagation does not cut it, as the the two (fyi & myi) reference
        # uncertainties are not independent (linked via scaling factor).
        # One could add the covariance matrix, but this is hardly justified
        # given the already very basic assumptions for fyi and myi density
        # uncertainty.
        # Here, the uncertainty is based on the fyi and myi density
        # uncertainties scaled with the myi fraction. To reflect the
        # myi_fraction uncertainty, the sea ice density uncertainty is
        # slightly increased in regions of higher myi_fraction uncertainty
        if "uncertainty" in self._options:
            unc = self._options.uncertainty
            rho_i_unc = unc.fyi_density
            rho_i_unc += myi_fraction * (unc.myi_density - unc.fyi_density)
            rho_i_unc += (unc.fyi_density - unc.myi_density) * (myi_fraction.uncertainty)

        # Use fixed value for uncertainty (or 0 if unspecified)
        else:
            unc = 0.0
            if "uncertainty_fixed_value" in self._options:
                unc = self._options.uncertainty_fixed_value
            rho_i_unc = np.full(rho_i.shape, unc, dtype=np.float32)

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
        deriv_fb = rho_w / (rho_w - rho_i)
        deriv_rho_s = l2.sd / (rho_w - rho_i)
        deriv_rho_i = (l2.frb * rho_w + l2.sd * l2.sdens) / (rho_w - rho_i)**2.
        deriv_sd = (l2.sdens) / (rho_w - rho_i)

        # Error propagation for statistical and systematic errors
        sit_unc = np.sqrt((deriv_fb**2. * l2.frb.uncertainty**2) +
                          (deriv_rho_i**2 * rho_i_unc**2) +
                          (deriv_sd**2. * l2.sd.uncertainty**2.) +
                          (deriv_rho_s**2 * l2.sdens.uncertainty**2.))

        return sit_unc


def snowfreeboard2thickness(frb, sd, rho_w, rho_i, rho_s):
    return rho_w / (rho_w - rho_i)*frb + (rho_s - rho_w)/(rho_w - rho_i)*sd


def frb2sit_errprop(frb, sd, rho_w, rho_i, rho_s, frb_unc, sd_unc, rho_i_unc, rho_s_unc):

    # Compute the partial derivates for the error propagation
    deriv_fb = rho_w/(rho_w-rho_i)
    deriv_rho_s = sd/(rho_w-rho_i)
    deriv_rho_i = (frb*rho_w + sd*rho_s) / (rho_w-rho_i)**2.
    deriv_sd = (rho_s)/(rho_w-rho_i)

    # Error propagation
    sit_unc = np.sqrt((deriv_fb**2. * frb_unc**2) +
                      (deriv_rho_i**2 * rho_i_unc**2) +
                      (deriv_sd**2. * sd_unc**2.) +
                      (deriv_rho_s**2 * rho_s_unc**2.))

    return sit_unc

def get_sit_algorithm(name):
    pyclass = globals().get(name, None)
    if pyclass is not None:
        return pyclass()
    else:
        return pyclass
