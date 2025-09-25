# -*- coding: utf-8 -*-

"""

"""

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"

import numpy as np
from scipy.ndimage import maximum_filter

from pysiral import psrlcfg
from pysiral.l3proc import Level3ProcessorItem


class Level3StatusFlag(Level3ProcessorItem):
    """
    A Level-3 processor item to compute the status flag
    """

    # Mandatory properties
    required_options = ["retrieval_status_target", "sic_thrs", "flag_values"]
    l2_variable_dependencies = []
    l3_variable_dependencies = ["sea_ice_concentration", "n_valid_waveforms", "landsea"]
    l3_output_variables = dict(status_flag=dict(dtype="i1", fill_value=1))

    def __init__(self, *args, **kwargs):
        """
        Compute surface type statistics
        :param args:
        :param kwargs:
        """
        super(Level3StatusFlag, self).__init__(*args, **kwargs)

    def apply(self):
        """
        Computes the status flag
        :return:
        """

        # Get the flag values from the l3 settings file
        flag_values = self.flag_values

        # Get status flag (fill value should be set to zero)
        sf = np.copy(self.l3grid.vars["status_flag"])

        # Init the flag with not data flag value
        sf[:] = flag_values["no_data"]

        # get input parameters
        par = np.copy(self.l3grid.vars[self.retrieval_status_target])
        sic = self.l3grid.vars["sea_ice_concentration"]
        nvw = self.l3grid.vars["n_valid_waveforms"]
        lnd = self.l3grid.vars["landsea"]

        # --- Compute conditions for flags ---

        # Get sea ice mask
        is_below_sic_thrs = np.logical_and(sic >= 0., sic < self.sic_thrs)

        # Get the pole hole information
        mission_ids = self.l3grid.metadata.mission_ids.split(",")
        orbit_inclinations = [psrlcfg.platforms.get_orbit_inclination(mission_id) for mission_id in mission_ids]

        # NOTE: due to varying grid cell size, it is no sufficient to just check which grid cell coordinate
        #       is outside the orbit coverage
        is_pole_hole = np.logical_and(
            np.abs(self.l3grid.vars["latitude"]) > (np.amin(orbit_inclinations)-1.0),
            flag_values["no_data"])

        # Check where the retrieval has failed
        is_land = lnd > 0
        has_data = nvw > 0
        has_retrieval = np.isfinite(par)
        retrieval_failed = np.logical_and(
            np.logical_and(has_data, np.logical_not(is_below_sic_thrs)),
            np.logical_not(has_retrieval))

        # Set sic threshold
        sf[np.where(is_below_sic_thrs)] = flag_values["is_below_sic_thrs"]

        # Set pole hole (Antarctica: Will be overwritten below)
        sf[np.where(is_pole_hole)] = flag_values["is_pole_hole"]

        # Set land mask
        sf[np.where(is_land)] = flag_values["is_land"]

        # Set failed retrieval
        sf[np.where(retrieval_failed)] = flag_values["retrieval_failed"]

        # Set retrieval successful
        sf[np.where(has_retrieval)] = flag_values["has_retrieval"]

        # Write Status flag
        self.l3grid.vars["status_flag"] = sf


class Level3QualityFlag(Level3ProcessorItem):
    """
    A Level-3 processor item to compute the status flag
    """

    # Mandatory properties
    required_options = ["add_rule_flags", "rules"]
    l2_variable_dependencies = []
    l3_variable_dependencies = ["sea_ice_thickness", "n_valid_waveforms", "negative_thickness_fraction",
                                "lead_fraction"]
    l3_output_variables = dict(quality_flag=dict(dtype="i1", fill_value=3))

    def __init__(self, *args, **kwargs):
        """
        Compute surface type statistics
        :param args:
        :param kwargs:
        """
        super(Level3QualityFlag, self).__init__(*args, **kwargs)

    def apply(self):
        """ Computation of quality flag indicator based on several rules defined in the l3 settings file """

        # Get the quality flag indicator array
        # This array will be continously updated by the quality check rules
        qif = np.copy(self.l3grid.vars["quality_flag"])
        sit = np.copy(self.l3grid.vars["sea_ice_thickness"])
        nvw = np.copy(self.l3grid.vars["n_valid_waveforms"])
        ntf = np.copy(self.l3grid.vars["negative_thickness_fraction"])
        lfr = np.copy(self.l3grid.vars["lead_fraction"])

        # As first step set qif to 1 where data is availabe
        qif[np.where(np.isfinite(sit))] = 0

        # Get a list of all the rules
        quality_flag_rules = self.rules.keys()

        # Simple way of handling rules (for now)

        # Use the Warren99 validity maslk
        # XXX: Not implemented yet
        if "qif_warren99_valid_flag" in quality_flag_rules:
            w99 = self.l3grid.vars["warren99_is_valid"]
            # mask = 0 means warren99 is invalid
            rule_options = self.rules["qif_warren99_valid_flag"]
            flag = np.full(qif.shape, 0, dtype=qif.dtype)
            flag[np.where(w99 == 0)] = rule_options["target_flag"]
            qif = np.maximum(qif, flag)

        # Elevate the quality flag for SARin or mixed SAR/SARin regions
        # (only sensible for CryoSat-2)
        if "qif_cs2_radar_mode_is_sin" in quality_flag_rules:
            radar_modes = self.l3grid.vars["radar_mode"]
            rule_options = self.rules["qif_cs2_radar_mode_is_sin"]
            flag = np.full(qif.shape, 0, dtype=qif.dtype)
            flag[np.where(radar_modes >= 2.)] = rule_options["target_flag"]
            qif = np.maximum(qif, flag)

        # Check the number of waveforms (less valid waveforms -> higher warning flag)
        if "qif_n_waveforms" in quality_flag_rules:
            flag = np.full(qif.shape, 0, dtype=qif.dtype)
            rule_options = self.rules["qif_n_waveforms"]
            for threshold, target_flag in zip(rule_options["thresholds"], rule_options["target_flags"]):
                flag[np.where(nvw < threshold)] = target_flag
            qif = np.maximum(qif, flag)

        # Check the availiability of leads in an area adjacent to the grid cell
        if "qif_lead_availability" in quality_flag_rules:
            flag = np.full(qif.shape, 0, dtype=qif.dtype)
            rule_options = self.rules["qif_lead_availability"]
            # get the window size
            grid_res = self.l3grid.griddef.resolution
            window_size = np.ceil(rule_options["search_radius_m"] / grid_res)
            window_size = int(2 * window_size + 1)
            # Use a maximum filter to get best lead fraction in area
            area_lfr = maximum_filter(lfr, size=window_size)
            thrs = rule_options["area_lead_fraction_minimum"]
            flag[np.where(area_lfr <= thrs)] = rule_options["target_flag"]
            qif = np.maximum(qif, flag)

        if "qif_miz_flag" in quality_flag_rules:
            flag = np.full(qif.shape, 0, dtype=qif.dtype)
            rule_options = self.rules["qif_miz_flag"]
            for source_flag, target_flag in zip(rule_options["source_flags"], rule_options["source_flags"]):
                flag[np.where(self.l3grid.vars["flag_miz"] == source_flag)] = target_flag
            qif = np.maximum(qif, flag)

        # Check the negative thickness fraction (higher value -> higher warnung flag)
        if "qif_high_negative_thickness_fraction" in quality_flag_rules:
            flag = np.full(qif.shape, 0, dtype=qif.dtype)
            rule_options = self.rules["qif_high_negative_thickness_fraction"]
            for threshold, target_flag in zip(rule_options["thresholds"], rule_options["target_flags"]):
                flag[np.where(ntf > threshold)] = target_flag
            qif = np.maximum(qif, flag)

        # Set all flags with no data to last flag value again
        qif[np.where(np.isnan(sit))] = 3

        # Set flag again
        self.l3grid.vars["quality_flag"] = qif
