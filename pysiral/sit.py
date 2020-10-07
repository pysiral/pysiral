# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 17:15:39 2016

@author: shendric
"""

import numpy as np

from pysiral.l2proc.procsteps import Level2ProcessorStep


class AlexandrovSeaIceDensity(Level2ProcessorStep):
    """
    A Level-2 processor step that adds sea ice density and sea ice density uncertainty
    to the L2 data object based on sea ice type classification
    """

    def __init__(self, *args, **kwargs):
        super(AlexandrovSeaIceDensity, self).__init__(*args, **kwargs)

    def execute_procstep(self, l1b, l2):
        """
        Mandatory method for the Level2ProcessorStep class. Will add sea ice density
        and uncertainty computed from the sea ice type (interpreted as multi-year ice fraction)
        and add output parameters to the L2 data object
        :param l1b:
        :param l2:
        :return: Error Status
        """

        # Get the mandatory error status output
        error_status = self.get_clean_error_status(l2.n_records)

        # Boundary conditions (all fyi, all myi)
        rho_i_fyi = self.cfg.options.fyi_density
        rho_i_myi = self.cfg.options.myi_density

        # Scales with MYI fraction
        myi_fraction = l2.sitype
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
        if "uncertainty" in self.cfg.options:
            unc = self.cfg.options.uncertainty
            rho_i_unc = unc.fyi_density
            rho_i_unc += myi_fraction * (unc.myi_density - unc.fyi_density)
            rho_i_unc += (unc.fyi_density - unc.myi_density) * myi_fraction.uncertainty

        # Use fixed value for uncertainty (or 0 if unspecified)
        else:
            unc = 0.0
            if "uncertainty_fixed_value" in self.cfg.options:
                unc = self.cfg.options.uncertainty_fixed_value
            rho_i_unc = np.full(rho_i.shape, unc, dtype=np.float32)

        # Add data to the Level-2 data object
        l2.set_auxiliary_parameter("idens", "sea_ice_density", rho_i, rho_i_unc)

        # Done
        return error_status

    @property
    def l2_input_vars(self):
        return ["sitype"]

    @property
    def l2_output_vars(self):
        return ["idens"]

    @property
    def error_bit(self):
        return self.error_flag_bit_dict["sit"]


class SeaIceFreeboard2SIT(Level2ProcessorStep):
    """
    Classical freeboard to thickness conversion assuming that the variable freeboard
    represent sea ice freeboard.
    """

    def __init__(self, *args, **kwargs):
        super(SeaIceFreeboard2SIT, self).__init__(*args, **kwargs)

    def execute_procstep(self, l1b, l2):
        """
        Compute sea ice thickness based on a Level-2 data object
        - Assuming freeboard is sea ice freeboard
        :param l1b:
        :param l2:
        :return:
        """

        # Short cuts
        rho_w = self.cfg.options.water_density

        # Compute sea ice thickness ...
        sit = self.func(l2.frb, l2.sd, rho_w, l2.idens, l2.sdens)

        # ... and sea ice thickness uncertainty
        sit_unc = self.uncfunc(l2.frb, l2.sd, rho_w, l2.idens, l2.sdens,
                               l2.frb.uncertainty, l2.sd.uncertainty,
                               l2.idens.uncertainty, l2.sdens.uncertainty)

        # Add parameter to Level-2 data object
        l2.sit.set_value(sit)
        l2.sit.set_uncertainty(sit_unc)

        # Remove any values for non-sea ice waveforms
        is_not_ice = np.logical_not(l2.surface_type.sea_ice.flag)
        l2.sit.set_nan_indices(np.where(is_not_ice))

        # Return the not sea-ice as error status flag
        return is_not_ice

    @property
    def l2_input_vars(self):
        return ["frb", "sd", "sdens", "idens"]

    @property
    def l2_output_vars(self):
        return ["sit"]

    @property
    def error_bit(self):
        return self.error_flag_bit_dict["sit"]

    @property
    def func(self):
        """
        Return the function for transforming sea-ice freeboard to thickness
        :return:
        """
        return icefreeboard2thickness

    @property
    def uncfunc(self):
        """
        Return the function for computing sea-ice thickness uncertainty
        :return:
        """
        return frb2sit_errprop


class SnowFreeboard2SIT(SeaIceFreeboard2SIT):
    """
    Variant of the SeaIceFreeboard2SIT, but for snow freeboard to thickness conversion.
    This is implemented as the children to SeaIceFreeboard2SIT, with only different
    pointer to the correct functions
    """

    def __init__(self, *args, **kwargs):
        super(SnowFreeboard2SIT, self).__init__(*args, **kwargs)

    @property
    def func(self):
        """
        Return the function for transforming sea-ice freeboard to thickness
        :return:
        """
        return snowfreeboard2thickness

    @property
    def uncfunc(self):
        """
        Return the function for computing sea-ice thickness uncertainty
        :return:
        """
        # TODO: Is this correct?
        return frb2sit_errprop


class L2SeaIceDraft(Level2ProcessorStep):
    """
    A Level-2 processor step class for computing sea ice draft and its uncertainty
    """

    def __init__(self, *args, **kwargs):
        super(L2SeaIceDraft, self).__init__(*args, **kwargs)

    def execute_procstep(self, l1b, l2):
        """
        API class for the Level-2 processor. Functionality is to compute sea ice draft and its uncertainty
        :param l1b: Level-1 data instance
        :param l2: A Level-2 data instance
        :return: None, Level-2 object is change in place
        """

        # Get error status
        error_status = self.get_clean_error_status(l2.n_records)

        # Compute sea-ice draft
        sea_ice_draft = l2.sit - l2.frb
        sea_ice_draft_uncertainty = np.sqrt(l2.sit.uncertainty**2. + l2.frb.uncertainty**2.)

        # Add to l2
        l2.set_auxiliary_parameter("sid", "sea_ice_draft", sea_ice_draft, uncertainty=sea_ice_draft_uncertainty)

        return error_status

    @property
    def l2_input_vars(self):
        return ["frb", "sit"]

    @property
    def l2_output_vars(self):
        return ["sid"]

    @property
    def error_bit(self):
        return self.error_flag_bit_dict["sit"]


def icefreeboard2thickness(frb, sd, rho_w, rho_i, rho_s):
    return rho_w/(rho_w - rho_i)*frb + rho_s/(rho_w - rho_i)*sd


def snowfreeboard2thickness(frb, sd, rho_w, rho_i, rho_s):
    return rho_w / (rho_w - rho_i)*frb + (rho_s - rho_w)/(rho_w - rho_i)*sd


def frb2sit_errprop(frb, sd, rho_w, rho_i, rho_s, frb_unc, sd_unc, rho_i_unc, rho_s_unc):

    # Compute the partial derivates for the error propagation
    deriv_fb = rho_w/(rho_w-rho_i)
    deriv_rho_s = sd/(rho_w-rho_i)
    deriv_rho_i = (frb*rho_w + sd*rho_s) / (rho_w-rho_i)**2.
    deriv_sd = rho_s/(rho_w-rho_i)

    # Error propagation
    sit_unc = np.sqrt((deriv_fb**2. * frb_unc**2) +
                      (deriv_rho_i**2 * rho_i_unc**2) +
                      (deriv_sd**2. * sd_unc**2.) +
                      (deriv_rho_s**2 * rho_s_unc**2.))

    return sit_unc

