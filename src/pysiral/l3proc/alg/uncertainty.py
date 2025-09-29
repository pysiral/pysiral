# -*- coding: utf-8 -*-

"""

"""

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"

import numpy as np
from loguru import logger
from typing import Dict, List, Tuple

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
    Computes the Level-3 uncertainties for radar freeboard, sea ice freeboard, sea ice thickness and sea ice draft.
    NOTE: This is a re-implementation of Level3GridUncertainties with more flexibility in terms of variable names
          and which uncertainties to compute.
    """

    def __init__(self, *args, **kwargs):
        """
        Compute surface type statistics
        :param args:
        :param kwargs:
        """
        self.var_name_dict = {}
        super(Level3GridUncertaintiesV2, self).__init__(*args, **kwargs)

    def _process_dynamic_dependencies(self) -> None:
        """
        Process dynamic dependencies based on the configuration.
        This method updates the required options, l2 and l3 variable dependencies,
        and l3 output variables based on the configuration provided during initialization.
        """

        # Update variable name dictionary
        self.var_name_dict = self.get_variable_name_dict()

    def apply(self):
        """
        Compute a level 3 uncertainty. The general idea is to compute the error propagation of average
        error components, where for components for random error the error of the l2 average
        is used and for systematic error components the average of the l2 error
        """

        if "radar_freeboard" in self.cfg:
            self.batch_func("radar_freeboard_uncertainty", **self.cfg["radar_freeboard"])

        if "sea_ice_freeboard" in self.cfg:
            self.batch_func("sea_ice_freeboard_uncertainty", **self.cfg["sea_ice_freeboard"])

        if "sea_ice_thickness" in self.cfg:
            self.batch_func("sea_ice_thickness_uncertainty", **self.cfg["sea_ice_thickness"])

        if "sea_ice_draft" in self.cfg:
            self.batch_func("sea_ice_draft_uncertainty", **self.cfg["sea_ice_draft"])

    def batch_func(self, func_name: str, **kwargs) -> None:
        """
        Batch processing function to apply a specific uncertainty computation method across all grid indices.

        :param func_name: Name of the uncertainty computation method to apply.
        :param kwargs: Additional keyword arguments to pass to the computation method.

        :return: None
        """
        func_map = {
            "radar_freeboard_uncertainty": self.compute_radar_freeboard_uncertainty,
            "sea_ice_freeboard_uncertainty": self.compute_sea_ice_freeboard_uncertainty,
            "sea_ice_thickness_uncertainty": self.compute_gridded_sea_ice_thickness_uncertainty,
            "sea_ice_draft_uncertainty": self.compute_sea_ice_draft_uncertainty
        }

        if func_name not in func_map:
            logger.error(f"Function '{func_name}' is not recognized.")
            return

        compute_func = func_map[func_name]
        for ix, iy in self.l3grid.grid_indices:
            compute_func(ix, iy, **kwargs)

    def get_variable_name_dict(self) -> Dict:
        """
        Get a dictionary with variable names for source and target variables.
        This allows to compute uncertainties for custom variable names and multiple solutions.

        If not specified, default variable names are used.

        :return: Variable name dictionary
        """
        required_vars = [
            "radar_freeboard_uncertainty",
            "radar_freeboard_l3_uncertainty",
            "sea_ice_freeboard",
            "sea_ice_freeboard_l3_uncertainty",
            "snow_depth",
            "snow_depth_uncertainty",
            "snow_density",
            "snow_density_uncertainty",
            "sea_ice_density",
            "sea_ice_density_uncertainty",
            "sea_ice_thickness",
            "sea_ice_thickness_l3_uncertainty",
            "sea_ice_draft",
            "sea_ice_draft_l3_uncertainty"
        ]
        var_name_dict = {var: var for var in required_vars}
        custom_source_variable_dict = self.cfg.get("source_variable_dict", {})
        var_name_dict.update(custom_source_variable_dict)
        return var_name_dict

    def compute_radar_freeboard_uncertainty(
            self,
            xi: int,
            yj: int,
            valid_uncertainty_range: Tuple[float, float] = None,
    ) -> None:
        """
        Compute the level-3 radar freeboard uncertainty for a given grid cell

        :param xi, x grid cell index
        :param yj, y grid cell index
        :param valid_uncertainty_range: tuple with (min, max) valid range for radar freeboard uncertainty

        :return: None: the result is written to the l3grid object
        """

        radar_freeboard_uncertainty_var = self.var_name_dict.get("radar_freeboard_uncertainty")
        target_variable = self.var_name_dict.get("radar_freeboard_l3_uncertainty")

        # Get the stack of radar freeboard uncertainty values and remove NaN's
        rfrb_uncs = np.array(self.l3grid.l2.stack[radar_freeboard_uncertainty_var][yj][xi])
        rfrb_uncs = rfrb_uncs[~np.isnan(rfrb_uncs)]

        if len(rfrb_uncs) == 0:
            self.l3grid.vars[target_variable][yj, xi] = np.nan
            return

        # Compute radar freeboard uncertainty as error or the mean from values with individual
        # error components (error of a weighted mean)
        weight = np.nansum(1. / rfrb_uncs ** 2.)
        rfrb_unc = 1. / np.sqrt(weight)

        if valid_uncertainty_range is not None:
            rfrb_unc = np.clip(rfrb_unc, valid_uncertainty_range[0], valid_uncertainty_range[1])

        self.l3grid.vars[target_variable][yj, xi] = rfrb_unc

    def compute_sea_ice_freeboard_uncertainty(
            self,
            xi: int,
            yj: int,
            snow_depth_correction_factor: float = 0.22,
            valid_uncertainty_range: Tuple[float, float] = None,
    ) -> None:
        """
        Compute the level-3 sea ice freeboard uncertainty for a given grid cell

        :param xi, x grid cell index
        :param yj, y grid cell index
        :param snow_depth_correction_factor: correction factor to account for snow loading on freeboard
        :param valid_uncertainty_range: tuple with (min, max) valid range for radar freeboard uncertainty

        :return: None: the result is written to the l3grid object
        """

        rfrc_unc_var_name = self.var_name_dict.get("radar_freeboard_l3_uncertainty")
        frb_var_name = self.var_name_dict.get("sea_ice_freeboard")
        sd_var_name = self.var_name_dict.get("snow_depth")
        sd_unc_var_name = self.var_name_dict.get("snow_depth_uncertainty")
        target_variable = self.var_name_dict.get("sea_ice_freeboard_l3_uncertainty")

        # Get parameters
        frb = self.l3grid.vars[frb_var_name][yj, xi]
        sd = self.l3grid.vars[sd_var_name][yj, xi]

        # Get uncertainties
        rfrb_unc = self.l3grid.vars[rfrc_unc_var_name][yj, xi]
        sd_unc = self.l3grid.vars[sd_unc_var_name][yj, xi]

        if np.isnan(frb) or np.isnan(sd) or np.isnan(rfrb_unc) or np.isnan(sd_unc):
            return

        # Calculate the level-3 freeboard uncertainty with updated radar freeboard uncertainty
        deriv_snow = snow_depth_correction_factor
        frb_unc = np.sqrt((deriv_snow * sd_unc) ** 2. + rfrb_unc ** 2.)

        if valid_uncertainty_range is not None:
            frb_unc = np.clip(frb_unc, valid_uncertainty_range[0], valid_uncertainty_range[1])

        self.l3grid.vars[target_variable][yj, xi] = frb_unc

    def compute_gridded_sea_ice_thickness_uncertainty(
            self,
            xi: int,
            yj: int,
            sea_water_density: float = 1024.0,
            valid_uncertainty_range: Tuple[float, float] = None,
    ) -> None:
        """
        Compute the level-3 uncertainty for a given grid cell

        :param xi, x grid cell index
        :param yj, y grid cell index
        :param sea_water_density: Seawater density [kg/m^3]
        :param valid_uncertainty_range: tuple with (min, max) valid range for sea ice thickness uncertainty

        :return: None: the result is written to the l3grid object
        """
        required_vars = [
            "radar_freeboard_uncertainty",
            "sea_ice_freeboard",
            "sea_ice_freeboard_l3_uncertainty",
            "snow_depth",
            "snow_density",
            "snow_depth_uncertainty",
            "snow_density_uncertainty",
            "sea_ice_density",
            "sea_ice_density_uncertainty",
            "sea_ice_thickness"
        ]
        values = {
            var_name: self.l3grid.vars[self.var_name_dict.get(var_name)][yj, xi]
            for var_name in required_vars
        }

        if np.isnan(values["sea_ice_thickness"]):
            return

        errprop_args = [
            values["sea_ice_freeboard"],
            values["snow_depth"],
            sea_water_density,
            values["sea_ice_density"],
            values["snow_density"],
            values["sea_ice_freeboard_l3_uncertainty"],
            values["snow_depth_uncertainty"],
            values["sea_ice_density_uncertainty"],
            values["snow_density_uncertainty"]
        ]
        sit_l3_unc = frb2sit_errprop(*errprop_args)

        # Cap the uncertainty
        # (very large values may appear in extreme cases)
        if valid_uncertainty_range is not None:
            sit_l3_unc = np.clip(sit_l3_unc, valid_uncertainty_range[0], valid_uncertainty_range[1])

        # Assign Level-3 uncertainty
        target_variable = self.var_name_dict.get("sea_ice_thickness_l3_uncertainty")
        self.l3grid.vars[target_variable][yj, xi] = sit_l3_unc

    def compute_sea_ice_draft_uncertainty(
            self,
            xi: int,
            yj: int,
            valid_uncertainty_range: Tuple[float, float] = None,
    ) -> None:
        """
        Compute the sea ice draft uncertainty for a given grid cell

        :param xi, x grid cell index
        :param yj, y grid cell index
        :param valid_uncertainty_range: tuple with (min, max) valid range for sea ice draft uncertainty

        :return: None: the result is written to the l3grid object
        """
        target_variable = self.var_name_dict.get("sea_ice_draft_l3_uncertainty")
        frb_unc = self.l3grid.vars[self.var_name_dict.get("sea_ice_freeboard_l3_uncertainty")][yj, xi]
        sit_l3_unc = self.l3grid.vars[self.var_name_dict.get("sea_ice_thickness_l3_uncertainty")][yj, xi]

        sid_l3_unc = np.sqrt(sit_l3_unc ** 2. + frb_unc ** 2.)

        if valid_uncertainty_range is not None:
            sid_l3_unc = np.clip(sid_l3_unc, valid_uncertainty_range[0], valid_uncertainty_range[1])

        self.l3grid.vars[target_variable][yj, xi] = sid_l3_unc

    @property
    def required_options(self) -> List[str]:
        return []

    @property
    def l2_variable_dependencies(self) -> List[str]:
        return [self.var_name_dict.get("radar_freeboard_uncertainty")]

    @property
    def l3_variable_dependencies(self) -> List[str]:
        dependencies = []
        if "sea_ice_freeboard" in self.cfg:
            dependencies.append(self.var_name_dict.get("sea_ice_freeboard"))
            dependencies.append(self.var_name_dict.get("snow_depth"))
            dependencies.append(self.var_name_dict.get("snow_depth_uncertainty"))
        if "sea_ice_thickness" in self.cfg:
            dependencies.append(self.var_name_dict.get("sea_ice_density"))
            dependencies.append(self.var_name_dict.get("sea_ice_density_uncertainty"))
            dependencies.append(self.var_name_dict.get("snow_density"))
            dependencies.append(self.var_name_dict.get("snow_density_uncertainty"))
            dependencies.append(self.var_name_dict.get("sea_ice_thickness"))
        if "sea_ice_draft" in self.cfg:
            dependencies.append(self.var_name_dict.get("sea_ice_draft"))
        return dependencies

    @property
    def l3_output_variables(self) -> Dict[str, Dict]:
        default_kwargs = {"dtype": "f4", "fill_value": np.nan}
        output_variables = {self.var_name_dict["radar_freeboard_l3_uncertainty"]: default_kwargs.copy()}
        if "sea_ice_freeboard" in self.cfg:
            output_variables[self.var_name_dict["sea_ice_freeboard_l3_uncertainty"]] = default_kwargs.copy()
        if "sea_ice_thickness" in self.cfg:
            output_variables[self.var_name_dict["sea_ice_thickness_l3_uncertainty"]] = default_kwargs.copy()
        if "sea_ice_draft" in self.cfg:
            output_variables[self.var_name_dict["sea_ice_draft_l3_uncertainty"]] = default_kwargs.copy()
        return output_variables

    @property
    def allow_overwrite(self) -> bool:
        """
        Repeated use may overwrite existing data
        :return:
        """
        return True
