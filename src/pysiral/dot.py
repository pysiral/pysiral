# -*- coding: utf-8 -*-

"""
@author: Stefan Hendricks

pysiral module for estimating dynamic ocean topography from existing sea-level anomaly data
child classes of pysiral.l2proc.procsteps.Level2ProcessorStep.

NOTES:

    1. The convention of the pysiral Level-2 processor is to compute dynamic ocean topography as

            dot = sla + mdt

       Thus, the auxiliary data set mean dynamic topography (mdt) is a mandatory auxiliary data
       set for the functionality of all classes in this module and the sla anomaly must be
       computed first.

"""

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"

import numpy as np

from pysiral.l2proc.procsteps import Level2ProcessorStep


class DynamicOceanTopography(Level2ProcessorStep):
    """
    A Level-2 processor step that adds sea ice density and sea ice density uncertainty
    to the L2 data object based on sea ice type classification
    """

    def __init__(self, *args, **kwargs):
        super(DynamicOceanTopography, self).__init__(*args, **kwargs)

    def execute_procstep(self, l1b, l2):
        """
        Mandatory Level-2 processor method that will execute the processing step and modify the
        L2 data object in-place.
        This method only sets the value of the dynamic ocean topography fields (dot) since
        DOT is a standard Level-2 data item.
        :param l1b:
        :param l2:
        :return:
        """

        # The dynamic ocean topography is either defined as
        #    1. DOT = SSH - geoid
        #    2. DOT = SLA + MDT
        # Here we use the second variant
        dot = l2.sla[:] + l2.mdt[:]
        dot_unc = np.sqrt((1./np.sqrt(3.)*l2.mdt.uncertainty[:])**2. + (l2.sla.uncertainty[:])**2.)
        l2.dot.set_value(dot)
        l2.dot.set_uncertainty(dot_unc)

        return np.isnan(l2.dot[:])

    @property
    def l2_input_vars(self):
        """
        Mandatory property for Level2ProcessorStep children
        :return: list (str)
        """
        return ["mdt", "sla"]

    @property
    def l2_output_vars(self):
        """
        Mandatory property for Level2ProcessorStep children
        :return: list (str)
        """
        return ["dot"]

    @property
    def error_bit(self):
        return self.error_flag_bit_dict["sla"]
