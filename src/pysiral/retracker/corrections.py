# -*- coding: utf-8 -*-

"""
This module contains range corrections associated with the retracking step
"""

import numpy as np
from loguru import logger

try:
    from samosa.help_functions import calc_ssb_j2Table
    SAMOSA_OK = True
except ImportError:
    SAMOSA_OK = False

from pysiral.l2proc.procsteps import Level2ProcessorStep


class ERSPulseDeblurring(Level2ProcessorStep):
    """
    A processing step applying the pulse deblurring correction of
    Peacock 1998 to retracked ranges.

    NOTE: Should only be used for ERS-1/2 data
    """

    def __init__(self, *args, **kwargs):
        super(ERSPulseDeblurring, self).__init__(*args, **kwargs)

    def execute_procstep(self, l1b, l2):
        """
        Compute the pulse deblurring correction

            Hcor = h + eps / 5. (eps < 0)

        based on the classifier data transferred from the l1p (epss: epsilon in seconds).
        :param l1b:
        :param l2:
        :return: error_status_flag
        """

        # Get a clean error status
        error_status = self.get_clean_error_status(l2.n_records)

        # Compute epsilon in meter (eps_m = eps_sec * c / 2.)
        eps = l2.epss * 0.5 * 299792458.

        # Compute and apply pulse deblurring correction
        slope = self.constant_slope
        pulse_deblurring_correction = np.array(eps < 0.).astype(float) * eps / slope
        for target_variable in self.target_variables:
            var = l2.get_parameter_by_name(target_variable)
            var[:] = var[:] + pulse_deblurring_correction
            l2.set_parameter(target_variable, var[:], var.uncertainty[:])

        # Add pulse deblurring correction to level-2 auxiliary data
        l2.set_auxiliary_parameter("pdbc", "pulse_deblurring_correction", pulse_deblurring_correction)

        # Return clean error status (for now)
        return error_status

    @property
    def l2_input_vars(self):
        return ["elev", "epss"]

    @property
    def l2_output_vars(self):
        return ["pdbc"]

    @property
    def target_variables(self):
        return self.cfg.options.get("target_variables", ["elevation"])

    @property
    def error_bit(self):
        return self.error_flag_bit_dict["retracker"]

    @property
    def constant_slope(self):
        # return slope term depending on ERS mission
        # defaults to 5.0, the ERS-2 value
        return self.cfg.options.get("slope", 5.0)


class SSBCorrectionJason2(Level2ProcessorStep):
    """
    A processing step applying a sea state bias correction.

    NOTE: Designed for use with SAMOSA retracker and taken from the Jason 2 mission

    Sea state bias model by Tran et al. (2012), CLS-CNES [1]

    [1]Tran N., S. Philipps, J.-C. Poisson, S. Urien, E. Bronner, N. Picot, "Impact of GDR_D standards on SSB
    corrections", Presentation OSTST2012 in Venice, http://www.aviso.altimetry.fr/fileadmin/
    documents/OSTST/2012/oral/02_friday_28/01_instr_processing_I/01_IP1_Tran.pdf

    """

    def __init__(self, *args, **kwargs):
        super(SSBCorrectionJason2, self).__init__(*args, **kwargs)
        if not SAMOSA_OK:
            logger.error("Unable to import the SAMOSA retracker helper functions. Has it been installed?")

    def execute_procstep(self, l1b, l2):
        """
        Compute the ssb correction.
        :param l1b:
        :param l2:
        :return: error_status_flag
        """

        if not SAMOSA_OK:
            raise ImportError("Unable to import the SAMOSA retracker helper functions. Has it been installed?")

        # Get a clean error status
        error_status = self.get_clean_error_status(l2.n_records)

        surface_lead = l2.surface_type.get_by_name('lead')
        surface_ocean = l2.surface_type.get_by_name('ocean')
        ssb_mask = np.bitwise_or(surface_lead.flag, surface_ocean.flag)

        if np.count_nonzero(ssb_mask) > 0:
            # Compute ssb
            try:
                # Need nan to num because SSB is -0.297281 for nan inputs!
                ssb = calc_ssb_j2Table(np.nan_to_num(l2.samswh), np.nan_to_num(l2.samwsp))
                logger.info('Computing sea state bias correction')
            except AttributeError:
                # This happens if no SAMOSA retracking because no ocean or lead records
                ssb = np.zeros(l2.n_records)

            # Change NaN values for correction to zero, and mask out non-ocean/lead
            ssb = np.nan_to_num(ssb, copy=False) * np.array(ssb_mask).astype(float)

            for target_variable in self.target_variables:
                var = l2.get_parameter_by_name(target_variable)
                var[:] = var[:] - ssb
                l2.set_parameter(target_variable, var[:], var.uncertainty[:])

            # Add sea state bias correction to level-2 auxiliary data
            l2.set_auxiliary_parameter("ssb", "sea_state_bias", ssb)

        # Return clean error status (for now)
        return error_status

    @property
    def l2_input_vars(self):
        return ["elev", "samswh", "samwsp"]

    @property
    def l2_output_vars(self):
        return ["ssb"]

    @property
    def target_variables(self):
        return self.cfg.options.get("target_variables", ["elevation"])

    @property
    def error_bit(self):
        return self.error_flag_bit_dict["retracker"]
