# -*- coding: utf-8 -*-

"""

"""

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"

import numpy as np
from typing import Tuple

from pysiral.l3proc import Level3ProcessorItem
from pysiral.sit import frb2sit_errprop


class Level3GridUncertainties(Level3ProcessorItem):
    """
    A Level-3 processor item to compute uncertainties of key geophysical variables on a grid.
    NOTE: As a concession to backward compability: sea ice draft uncertainty will be computed, but
          the sea ice draft is not a required input parameter
    """

    # Mandatory properties
    required_options = ["water_density", "snow_depth_correction_factor", "max_l3_uncertainty"]
    l2_variable_dependencies = ["radar_freeboard_uncertainty"]
    l3_variable_dependencies = ["sea_ice_freeboard", "snow_depth",
                                "snow_density", "snow_depth_uncertainty", "snow_density_uncertainty"]
    l3_output_variables = dict(radar_freeboard_l3_uncertainty=dict(dtype="f4", fill_value=np.nan),
                               freeboard_l3_uncertainty=dict(dtype="f4", fill_value=np.nan),
                               sea_ice_thickness_l3_uncertainty=dict(dtype="f4", fill_value=np.nan),
                               sea_ice_draft_l3_uncertainty=dict(dtype="f4", fill_value=np.nan))

    def __init__(self, *args, **kwargs):
        """
        Compute surface type statistics
        :param args:
        :param kwargs:
        """
        super(Level3GridUncertainties, self).__init__(*args, **kwargs)

    def apply(self):
        """ Compute a level 3 uncertainty. The general idea is to compute the error propagation of average
        error components, where for components for random error the error of the l2 average
        is used and for systematic error components the average of the l2 error """

        # Loop over grid items
        map(self.compute_gridded_uncertainty, self.l3grid.grid_indices)

    def compute_gridded_uncertainty(self, grid_index: Tuple[int, int]) -> None:
        """
        Compute the level-3 uncertainty for a given grid cell

        :param grid_index: (xi, yj) grid cell index tuple

        :return: None: the result is written to the l3grid object
        """

        xi, yj = grid_index

        # Options
        rho_w = self.water_density
        sd_corr_fact = self.snow_depth_correction_factor

        # Get random uncertainty
        # Note: this applies only to the radar freeboard uncertainty.
        #       Thus, we need to recalculate the sea ice freeboard uncertainty

        # Get the stack of radar freeboard uncertainty values and remove NaN's
        # rfrb_unc = self.l3["radar_freeboard_uncertainty"][yj, xi]
        rfrb_uncs = np.array(self.l3grid.l2.stack["radar_freeboard_uncertainty"][yj][xi])
        rfrb_uncs = rfrb_uncs[~np.isnan(rfrb_uncs)]

        # Compute radar freeboard uncertainty as error or the mean from values with individual
        # error components (error of a weighted mean)
        weight = np.nansum(1. / rfrb_uncs ** 2.)
        rfrb_unc = 1. / np.sqrt(weight)
        self.l3grid.vars["radar_freeboard_l3_uncertainty"][yj, xi] = rfrb_unc

        # Get parameters
        frb = self.l3grid.vars["sea_ice_freeboard"][yj, xi]
        sd = self.l3grid.vars["snow_depth"][yj, xi]
        rho_s = self.l3grid.vars["snow_density"][yj, xi]

        # Get systematic error components
        sd_unc = self.l3grid.vars["snow_depth_uncertainty"][yj, xi]
        rho_s_unc = self.l3grid.vars["snow_density_uncertainty"][yj, xi]

        # Calculate the level-3 freeboard uncertainty with updated radar freeboard uncertainty
        deriv_snow = sd_corr_fact
        frb_unc = np.sqrt((deriv_snow * sd_unc) ** 2. + rfrb_unc ** 2.)
        self.l3grid.vars["freeboard_l3_uncertainty"][yj, xi] = frb_unc

        # Calculate the level-3 thickness uncertainty

        # # Check of sea ice thickness exists
        if "sea_ice_thickness" not in self.l3grid.vars:
            return

        if np.isnan(self.l3grid.vars["sea_ice_thickness"][yj, xi]):
            return

        rho_i = self.l3grid.vars["sea_ice_density"][yj, xi]
        rho_i_unc = self.l3grid.vars["sea_ice_density_uncertainty"][yj, xi]
        errprop_args = [frb, sd, rho_w, rho_i, rho_s, frb_unc, sd_unc, rho_i_unc, rho_s_unc]
        sit_l3_unc = frb2sit_errprop(*errprop_args)

        # Cap the uncertainty
        # (very large values may appear in extreme cases)
        if sit_l3_unc > self.max_l3_uncertainty:
            sit_l3_unc = self.max_l3_uncertainty

        # Assign Level-3 uncertainty
        self.l3grid.vars["sea_ice_thickness_l3_uncertainty"][yj, xi] = sit_l3_unc

        # Compute sea ice draft uncertainty
        if "sea_ice_draft" not in self.l3grid.vars:
            return

        sid_l3_unc = np.sqrt(sit_l3_unc ** 2. + frb_unc ** 2.)
        self.l3grid.vars["sea_ice_draft_l3_uncertainty"][yj, xi] = sid_l3_unc


class Level3GridUncertaintiesV2(Level3ProcessorItem):
    """
    A Level-3 processor item to compute uncertainties of key geophysical variables on a grid.
    NOTE: As a concession to backward compability: sea ice draft uncertainty will be computed, but
          the sea ice draft is not a required input parameter
    """

    # Mandatory properties
    required_options = ["water_density", "snow_depth_correction_factor", "max_l3_uncertainty"]
    l2_variable_dependencies = ["radar_freeboard_uncertainty"]
    l3_variable_dependencies = ["sea_ice_freeboard", "snow_depth",
                                "snow_density", "snow_depth_uncertainty", "snow_density_uncertainty"]
    l3_output_variables = dict(radar_freeboard_l3_uncertainty=dict(dtype="f4", fill_value=np.nan),
                               freeboard_l3_uncertainty=dict(dtype="f4", fill_value=np.nan),
                               sea_ice_thickness_l3_uncertainty=dict(dtype="f4", fill_value=np.nan),
                               sea_ice_draft_l3_uncertainty=dict(dtype="f4", fill_value=np.nan))

    def __init__(self, *args, **kwargs):
        """
        Compute surface type statistics
        :param args:
        :param kwargs:
        """
        super(Level3GridUncertaintiesV2, self).__init__(*args, **kwargs)

    def apply(self):
        """ Compute a level 3 uncertainty. The general idea is to compute the error propagation of average
        error components, where for components for random error the error of the l2 average
        is used and for systematic error components the average of the l2 error """

        # Loop over grid items
        map(self.compute_gridded_uncertainty, self.l3grid.grid_indices)

    def compute_gridded_uncertainty(self, grid_index: Tuple[int, int]) -> None:
        """
        Compute the level-3 uncertainty for a given grid cell

        :param grid_index: (xi, yj) grid cell index tuple

        :return: None: the result is written to the l3grid object
        """

        xi, yj = grid_index

        # Options
        rho_w = self.water_density
        sd_corr_fact = self.snow_depth_correction_factor

        # Get random uncertainty
        # Note: this applies only to the radar freeboard uncertainty.
        #       Thus, we need to recalculate the sea ice freeboard uncertainty

        # Get the stack of radar freeboard uncertainty values and remove NaN's
        # rfrb_unc = self.l3["radar_freeboard_uncertainty"][yj, xi]
        rfrb_uncs = np.array(self.l3grid.l2.stack["radar_freeboard_uncertainty"][yj][xi])
        rfrb_uncs = rfrb_uncs[~np.isnan(rfrb_uncs)]

        # Compute radar freeboard uncertainty as error or the mean from values with individual
        # error components (error of a weighted mean)
        weight = np.nansum(1. / rfrb_uncs ** 2.)
        rfrb_unc = 1. / np.sqrt(weight)
        self.l3grid.vars["radar_freeboard_l3_uncertainty"][yj, xi] = rfrb_unc

        # Get parameters
        frb = self.l3grid.vars["sea_ice_freeboard"][yj, xi]
        sd = self.l3grid.vars["snow_depth"][yj, xi]
        rho_s = self.l3grid.vars["snow_density"][yj, xi]

        # Get systematic error components
        sd_unc = self.l3grid.vars["snow_depth_uncertainty"][yj, xi]
        rho_s_unc = self.l3grid.vars["snow_density_uncertainty"][yj, xi]

        # Calculate the level-3 freeboard uncertainty with updated radar freeboard uncertainty
        deriv_snow = sd_corr_fact
        frb_unc = np.sqrt((deriv_snow * sd_unc) ** 2. + rfrb_unc ** 2.)
        self.l3grid.vars["freeboard_l3_uncertainty"][yj, xi] = frb_unc

        # Calculate the level-3 thickness uncertainty

        # # Check of sea ice thickness exists
        if "sea_ice_thickness" not in self.l3grid.vars:
            return

        if np.isnan(self.l3grid.vars["sea_ice_thickness"][yj, xi]):
            return

        rho_i = self.l3grid.vars["sea_ice_density"][yj, xi]
        rho_i_unc = self.l3grid.vars["sea_ice_density_uncertainty"][yj, xi]
        errprop_args = [frb, sd, rho_w, rho_i, rho_s, frb_unc, sd_unc, rho_i_unc, rho_s_unc]
        sit_l3_unc = frb2sit_errprop(*errprop_args)

        # Cap the uncertainty
        # (very large values may appear in extreme cases)
        if sit_l3_unc > self.max_l3_uncertainty:
            sit_l3_unc = self.max_l3_uncertainty

        # Assign Level-3 uncertainty
        self.l3grid.vars["sea_ice_thickness_l3_uncertainty"][yj, xi] = sit_l3_unc

        # Compute sea ice draft uncertainty
        if "sea_ice_draft" not in self.l3grid.vars:
            return

        sid_l3_unc = np.sqrt(sit_l3_unc ** 2. + frb_unc ** 2.)
        self.l3grid.vars["sea_ice_draft_l3_uncertainty"][yj, xi] = sid_l3_unc
