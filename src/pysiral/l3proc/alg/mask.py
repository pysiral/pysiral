# -*- coding: utf-8 -*-

"""

"""

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"

import numpy as np
import xarray as xr

from loguru import logger
from pathlib import Path

from pysiral import psrlcfg
from pysiral.l3proc import Level3ProcessorItem
from pysiral.mask import L3Mask


class Level3LoadMasks(Level3ProcessorItem):
    """
    A Level-3 processor item to load external masks
    """

    # Mandatory properties
    required_options = ["mask_names"]
    l2_variable_dependencies = []
    l3_variable_dependencies = []
    # Note: the output names depend on mask name, thus these will be
    #       created in apply (works as well)
    l3_output_variables = dict()

    def __init__(self, *args, **kwargs):
        """
        Compute surface type statistics
        :param args:
        :param kwargs:
        """
        super(Level3LoadMasks, self).__init__(*args, **kwargs)

    def apply(self):
        """
        Load masks and add them as grid variable (variable name -> mask name)
        :return:
        """

        # Get the mask names and load each
        for mask_name in self.mask_names:

            # The masks are stored in external files that can be automatically
            #  found with the grid id
            mask = L3Mask(mask_name, self.l3grid.griddef.grid_id)

            # Add the mask to the l3grid as variable
            if not mask.error.status:
                self.l3grid.add_grid_variable(mask_name, np.nan, mask.mask.dtype)
                self.l3grid.vars[mask_name] = mask.mask

            # If fails, only add an empty variable
            else:
                self.l3grid.add_grid_variable(mask_name, np.nan, "f4")
                error_msgs = mask.error.get_all_messages()
                for error_msg in error_msgs:
                    logger.error(error_msg)


class Level3LoadCCILandMask(Level3ProcessorItem):
    """
    A Level-3 processor item to load the CCI land mask
    """

    # Mandatory properties
    required_options = ["local_machine_def_mask_tag", "mask_name_dict"]
    l2_variable_dependencies = []
    l3_variable_dependencies = []
    # Note: the output names depend on mask name, thus these will be
    #       created in apply (works as well)
    l3_output_variables = dict()

    def __init__(self, *args, **kwargs):
        """
        Initiate the class
        :param args:
        :param kwargs:
        """
        super(Level3LoadCCILandMask, self).__init__(*args, **kwargs)

    def apply(self):
        """
        Load masks and add them as grid variable (variable name -> mask name)
        :return:
        """

        # Short cut
        grid_id = self.l3grid.griddef.grid_id

        # Get mask target path:
        mask_tag = self.cfg["local_machine_def_mask_tag"]
        lookup_directory = psrlcfg.local_machine.auxdata_repository.mask.get(mask_tag, None)
        if lookup_directory is None:
            msg = "Missing local machine def tag: auxdata_repository.mask.{}".format(mask_tag)
            self.error.add_error("invalid-local-machine-def", msg)
            logger.error(msg)
            return
        lookup_directory = Path(lookup_directory)

        # Get the mask target filename
        try:
            filename = self.cfg["mask_name_dict"][grid_id.replace("_", "").lower()]
        except KeyError:
            logger.error(f"Could not find mask for grid id {grid_id} -> aborting")
            return
        mask_filepath = lookup_directory / filename
        if not mask_filepath.is_file():
            msg = "Missing input file: {}".format(mask_filepath)
            self.error.add_error("invalid-local-machine-def", msg)
            logger.error(msg)
            return

        # Load the data and extract the flag
        nc = xr.load_dataset(str(mask_filepath), decode_times=False)

        # The target land sea flag should be 1 for land and 0 for sea,
        # but the CCI landsea mask provides also fractional values for
        # mixed surfaces types. Thus, we add two arrays to the
        # L3grid object
        #   1. a classical land/sea mask with land:1 and sea: 0. In
        #      this notation the mixed pixels are attributed to sea
        #      because there might be some valid retrieval there
        #   2. the ocean density value as is
        density_of_ocean = np.flipud(nc.density_of_ocean.values)
        landsea_mask = density_of_ocean < 1e-5

        # Add mask to l3 grid
        mask_variable_name = self.cfg["mask_variable_name"]
        self.l3grid.add_grid_variable(mask_variable_name, np.nan, landsea_mask.dtype)
        self.l3grid.vars[mask_variable_name] = landsea_mask.astype(int)

        density_variable_name = self.cfg["density_variable_name"]
        self.l3grid.add_grid_variable(density_variable_name, np.nan, nc.density_of_ocean.values.dtype)
        self.l3grid.vars[density_variable_name] = density_of_ocean


class Level3ParameterMask(Level3ProcessorItem):
    """
    A Level-3 processor item to load external masks
    """

    # Mandatory properties
    required_options = ["source", "condition", "targets"]
    l2_variable_dependencies = []
    l3_variable_dependencies = []
    # Note: the output names depend on mask name, thus these will be
    #       created in apply (works as well)
    l3_output_variables = dict()

    def __init__(self, *args, **kwargs):
        """
        Compute surface type statistics
        :param args:
        :param kwargs:
        """
        super(Level3ParameterMask, self).__init__(*args, **kwargs)

    def apply(self):
        """
        Mask certain parameters based on condition of one other parameter
        :return:
        """

        # Get the source parameter
        source = self.l3grid.vars[self.source]

        # Compute the masking condition
        conditions = self.condition.split(";")
        n_conditions = len(conditions)

        if n_conditions == 0:
            msg = "Missing condition in %s" % self.__class__.__name__
            self.error.add_error("invalid-l3mask-def", msg)
            return

        # Start with the first (and maybe only condition)
        filter_mask = self._get_l3_mask(source, conditions[0], self.cfg)

        # Add conditions
        if n_conditions >= 2:
            for i in range(1, n_conditions):
                new_filter = self._get_l3_mask(source, conditions[i], self.cfg)
                if self.cfg["connect_conditions"] == "or":
                    filter_mask = np.logical_or(filter_mask, new_filter)
                elif self.cfg["connect_conditions"] == "and":
                    filter_mask = np.logical_and(filter_mask, new_filter)
                else:
                    msg = "Invalid l3 mask operation: %s"
                    msg %= self.cfg["connect_conditions"]
                    self.error.add_error("invalid-l3mask-def", msg)
                    self.error.raise_on_error()

        # Apply mask
        masked_indices = np.where(filter_mask)
        for target in self.targets:
            try:
                self.l3grid.vars[target][masked_indices] = np.nan
            except ValueError:
                if self.l3grid.vars[target].dtype.kind == "i":
                    self.l3grid.vars[target][masked_indices] = -1
                else:
                    msg = "Cannot set nan (or -1) as mask value to parameter: %s " % target
                    logger.warning(msg)

    def _get_l3_mask(self, source_param, condition, options):
        """ Return bool array based on a parameter and a predefined
        masking operation """
        if condition.strip() == "is_nan":
            return np.isnan(source_param)
        elif condition.strip() == "is_zero":
            return np.array(source_param <= 1.0e-9)
        elif condition.strip() == "is_smaller":
            return np.array(source_param < options["is_smaller_threshold"])
        else:
            msg = "Unknown condition in l3 mask: %s" % condition
            self.error.add_error("invalid-l3mask-condition", msg)
            self.error.raise_on_error()
