# -*- coding: utf-8 -*-

"""
"""

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"


import re
import numpy as np
import xarray as xr
from loguru import logger
from pathlib import Path
from typing import Any, Dict, Union

from pysiral.l1data import Level1bData
from pysiral.l1preproc.procitems import L1PProcItem


class CryoTEMPOExternalTides(L1PProcItem):
    """
    A class to update tides specifically for the production of Cryo-TEMPO data.
    Here, the tide information is extracted from the FES2022b tide model files with identical
    temporal coverage as CryoSat-2 ICE Level-1 files. Thus, only the correct file has to
    be identified and then mapped to the Level-1 data object.

    NOTES:
        - This Level-1 processor item is purpose built for the Cryo-TEMPO
          production workflow.
        - This Level-1 processor item can only be run at stage `post_source_file`,
          because otherwise the mapping strategy does not work.
    """

    def __init__(self, **cfg: Dict[str, Any]):
        super(CryoTEMPOExternalTides, self).__init__(**cfg)

    def apply(self, l1: Level1bData) -> None:
        """
        Replace the tide values in the l1b data object with the FES2022b tide model.
        This is done by reading tide data extracted from a file with the same coverage
        as the

        :param l1: The Level-1 data container

        :raises ValueError: Class not applied to CryoSat-2 L1B netcdf
        :raises FileNotFoundError: Expected FES2022b tide file not found and setting
            `raise_on_missing_file` is set to True in the configuration

        :return: None, Level-1 data object is updated in place
        """

        # Ensure that the Level-1 data is only based on one CryoSat-2 file
        if len(l1.source_filepaths) != 1 and l1.info.mission != "cryosat2":
            msg = "CryoTEMPOExternalTides can only be applied to Level-1 data objects based on one CryoSat-2 file"
            logger.error(msg)
            if self.cfg.get("raise_on_missing_file", False):
                raise ValueError(msg)
            return

        # Construct the expected filename
        tides_file_path = self.get_tide_filepath(l1.source_filepaths[0])

        # Test if file exists
        logger.info(f"- Looking for tide file {tides_file_path} to update tide values in Level-1 data object")
        if not tides_file_path.is_file():
            msg = f"- Expected tide file {tides_file_path} not found. Cannot update tide values."
            logger.error(msg)
            if self.cfg.get("raise_on_missing_file", False):
                raise FileNotFoundError(msg)
            return

        # Map the tide variables
        self.map_external_tides(tides_file_path, l1)

    def get_tide_filepath(self, filepath: Union[Path, str]) -> Path:
        """
        Construct the file path for the FES2022 data from the stored L1b file path.
        This is done by replacing parts of the file name and the necessary information
        needs to be specified in the Level-1 pre-processor configuration file.

        :param filepath: Full file path of the CryoSat-2 L1B netCDF file

        :return: Fulle file path for the FES2022 tide data file to be mapped the L1 data object
        """
        expected_fes_file_path = str(filepath)
        for tag_orig, tag_replacement in self.cfg.get("file_name_mapping", []):
            p = re.compile(tag_orig)
            expected_fes_file_path = p.sub(tag_replacement, expected_fes_file_path)
        return Path(str(expected_fes_file_path))

    def map_external_tides(self, fes_file_path: Path, l1: Level1bData) -> None:
        """
        Reads the netCDF file with the FES2022b tide data and maps the tide values to the Level-1 data object.

        :param fes_file_path: The FES2022 netCDF file path
        :param l1: The Level-1 data object

        :raises ValueError: Mismatch of dimensions between Level-1 data object and FES file

        :return: None, Level-1 data object is changed in place
        """
        ds = xr.open_dataset(fes_file_path, decode_times=False, mask_and_scale=True)
        tide_variable_name_mapping = self.cfg.get("tide_variable_name_mapping", {})

        for l1_var_name, external_var_name in tide_variable_name_mapping.items():

            # Get both tides
            external_tide_value = ds.variables[external_var_name].values
            l1_tide_value = l1.get_parameter_by_name("correction", l1_var_name)

            # Dimensions must match (no full check of actual times)
            if external_tide_value.size != l1_tide_value.size:
                msg = (
                    f"Dimension mismatch between external tide variable {external_tide_value} "
                    f"({external_tide_value.size}) and Level-1 tide variable {l1_var_name} "
                    f"({l1_tide_value.size}), Skipping tide update"
                )
                logger.error(msg)
                if self.cfg.get("raise_on_missing_file", False):
                    raise ValueError(msg)
                return

            # Check if NaNs exist in external data
            nans_indices = np.where(np.isnan(external_tide_value))[0]
            if len(nans_indices) > 0:
                external_tide_value = self.handle_nan_values(external_tide_value, nans_indices, l1_tide_value)

            # All checks complete
            l1.set_parameter_by_name("correction", l1_var_name, external_tide_value)

    def handle_nan_values(
            self,
            tide_file_value: np.ndarray,
            nans_indices: np.ndarray,
            l1_tide_value: np.ndarray
    ) -> np.ndarray:
        """
        Handle NaN values in the external tide data. This is done by replacing the NaN values with the corresponding
        values from the Level-1 data object.

        :param tide_file_value: The tide values from the external file
        :param nans_indices: The indices of the NaN values in the external tide data
        :param l1_tide_value: The tide values from the Level-1 data object

        :return: The updated tide values with NaNs handled
        """
        nan_policy = self.cfg.get("tide_file_nan_policy", "ignore")
        match nan_policy:
            case "ignore": return tide_file_value
            case "fill": tide_file_value[nans_indices] = l1_tide_value[nans_indices]
            case "interpolate": raise NotImplementedError("Interpolation of NaN values in tide file not implemented")
            case _: raise ValueError(f"Unknown NaN policy {nan_policy} for tide file")
        return tide_file_value
