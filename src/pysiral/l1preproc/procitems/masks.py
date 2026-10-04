# -*- coding: utf-8 -*-
"""

"""

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"

from pathlib import Path
import numpy as np
import xarray as xr
import numpy.typing as npt
from loguru import logger
from typing import Tuple, Union, Optional, Literal

from pydantic import BaseModel, computed_field

from pysiral import psrlcfg

from pysiral.core.flags import SURFACE_TYPE_DICT
from pysiral.grid import GridTrajectoryExtract

from pysiral.l1data import Level1bData
from pysiral.l1preproc.procitems import L1PProcItem


class L1PHighResolutionLandMask(L1PProcItem):
    """
    Level-1 processor item providing access to a high resolution
    land mask and distance to land fields
    """

    def __init__(self, **cfg):
        """
        Initialize the class. This step includes parsing the static mask
        and keeping it in memory
        :param cfg:
        """
        super(L1PHighResolutionLandMask, self).__init__(**cfg)

        # Get the local file path
        type_ = self.cfg.get("local_machine_def_auxclass")
        tag = self.cfg.get("local_machine_def_tag")
        filename = self.cfg.get("filename")
        lookup_directory = psrlcfg.local_machine.auxdata_repository[type_][tag]

        # Set the file for each hemisphere type
        mask_filepath = Path(lookup_directory) / str(filename)
        with xr.open_dataset(mask_filepath) as nc:
            self.grid_def = self.cfg.get("grid_def")
            self.land_ocean_flag_grid = nc.land_ocean_flag.values
            self.distance_to_coast_grid = nc.distance_to_coast.values

    def apply(self, l1: Level1bData) -> None:
        """
        Extract land/ocean flag and distance to coast along the l1p trajectory if a mask exists
        for the corresponding hemisphere of the l1p data object.

        The parameters are stored in the classifier data group among the original surface type
        value, which is then update for the mask coverage

        :param l1:
        :return: None
        """

        # Get the land/ocean flag from the grid
        land_ocean_flag, distance_to_coast = self.get_trajectory(l1.time_orbit.longitude, l1.time_orbit.latitude)

        # TODO: Update when global land/ocean flag becomes available
        # NOTE: Currently only a grid for the northern hemisphere is available, thus
        #       there will be regions not covered by the gridded land/ocean flag.
        #       That's why the built/in land/ocean flag will only be updated for
        #       region with coverage and be left untouched everywhere else.
        #       The arrays are added to the l1 classifier container regardless
        #       for consistency.

        # --- Update the L1 data container ---
        # 1. Save both extracted variables to classifier data groups
        l1.classifier.add(land_ocean_flag, "hr_land_ocean_flag")
        l1.classifier.add(distance_to_coast, "distance_to_coast")

        # 2. Save original ESA surface type variable to classifier data group
        #    and update the surface type flag in the surface type data group
        l1.classifier.add(l1.surface_type.flag, "orig_land_ocean_flag")

        # 3. Update the surface type instance
        valid_mask_indices = land_ocean_flag != self.cfg.get("dummy_val")["land_ocean_flag"]
        flag_update = np.full(l1.n_records, SURFACE_TYPE_DICT["invalid"])
        flag_update[land_ocean_flag == 1] = SURFACE_TYPE_DICT["land"]
        flag_update[land_ocean_flag == 0] = SURFACE_TYPE_DICT["ocean"]
        updated_surface_type_flag = l1.surface_type.flag.copy()
        updated_surface_type_flag[valid_mask_indices] = flag_update[valid_mask_indices]
        l1.surface_type.set_flag(updated_surface_type_flag)

    def get_trajectory(self, longitude: npt.NDArray, latitude: npt.NDArray) -> Tuple[npt.NDArray, npt.NDArray]:
        """
        Get extract of the land/ocean flag and the distance to coast value for an
        array of (longitude, latitude) positions.

        If longitude, latitude is outside the grid, a pre-defined dummy value will
        be returned.

        :param longitude: Longitude values in degrees
        :param latitude: latitude values in degrees

        :raises None:

        :return: land ocean flag & distance to coast values for longitude, latitude positions
        """
        grid2track = GridTrajectoryExtract(longitude, latitude, self.grid_def)
        land_ocean_flag = grid2track.get_from_grid_variable(
            self.land_ocean_flag_grid,
            outside_value=self.dummy_val["land_ocean_flag"]
        )

        distance_to_coast = grid2track.get_from_grid_variable(
            self.distance_to_coast_grid,
            outside_value=self.dummy_val["distance_to_coast"]
        )
        return land_ocean_flag, distance_to_coast



class GridDefinition(BaseModel):
    center_latitude: float

    @computed_field
    def projection(self) -> dict[str, Union[str, float]]:
        return {"proj": "laea", "lon_0": 0.0, "lat_0": self.center_latitude}

    @property
    def dimension(self) -> dict[str, Union[int, float]]:
        return {"n_cols": 43200, "n_lines": 43200, "dx": 250, "dy": 250}


class SurfaceClassificationMaskConfig(BaseModel):
    local_machine_def_auxclass: str = "mask"
    local_machine_def_tag: str = "cryotempo_surface_classification"
    version: str = "v1p0"
    dummy_val: dict[str, int] = {"land_ocean_flag": -1}
    filename_template: str = "cryotempo-surface_classification-{hemisphere}-250m-{version}.nc"
    grid_def: dict[Literal["nh", "sh"], GridDefinition] = {
        "nh": GridDefinition(center_latitude=90),
        "sh": GridDefinition(center_latitude=-90)
    }

    @computed_field
    def lookup_directory(self) -> Optional[Path]:
        try:
            directory = psrlcfg.local_machine.auxdata_repository[self.local_machine_def_auxclass][self.local_machine_def_tag]
        except KeyError:
            logger.error(
                f"Could not find lookup directory for auxclass '{self.local_machine_def_auxclass}' and tag "
                f"'{self.local_machine_def_tag}' in local machine configuration.")
            return None
        directory = Path(directory)
        if not directory.exists():
            logger.error(f"Lookup directory '{directory}' does not exist.")
            return None
        return directory

class SurfaceClassificationData(object):

    def __init__(self, ds: xr.Dataset, grid_def: GridDefinition) -> None:
        self.grid_def = grid_def
        self.ds = ds

    @classmethod
    def from_file(cls, filepath: Path, grid_def: GridDefinition) -> "SurfaceClassificationData":
        ds = xr.load_dataset(filepath)
        return cls(ds=ds, grid_def=grid_def)


class L1PCryoTEMPOSurfaceClassification(L1PProcItem):
    """
    Level-1 processor item providing access to a high resolution
    surface classification mask for both Arctic and Antarctic hemispheres.

    The mask provides a land/ocean flag, a monthly surface classification flag,
    that accounts for climatological sea ice extent and sea ice region codes.
    """

    def __init__(self, **cfg: dict) -> None:
        """
        Initialize the class. This step includes parsing the static mask
        and keeping it in memory

        :param cfg: The configuration dictionary (usually from L1 preprocessor setting file)
        """
        super(L1PCryoTEMPOSurfaceClassification, self).__init__(**{})

        # Overwrite the configuration with the provided cfg
        self.cfg = SurfaceClassificationMaskConfig(**cfg)

        # Set the file for each hemisphere type
        # TODO: At the moment both masks are loaded in the memory, as pre-processor settings are unknown at this point.
        #       This could be changed by giving preprocessor context to the L1PProcItem class initialization,
        #       so that only needed data is loaded.
        self.data: dict[str, SurfaceClassificationData] = {}
        for hemisphere in ["nh", "sh"]:
            mask_filepath = self.cfg.lookup_directory / str(self.cfg.filename_template).format(
                hemisphere=hemisphere, version=self.cfg.version
            )
            logger.debug(f"Loading mask from {mask_filepath}")
            self.data[hemisphere] = SurfaceClassificationData.from_file(mask_filepath, self.cfg.grid_def[hemisphere])


    def apply(self, l1: Level1bData) -> None:
        """
        Extract land/ocean flag and distance to coast along the l1p trajectory if a mask exists
        for the corresponding hemisphere of the l1p data object.

        The parameters are stored in the classifier data group among the original surface type
        value, which is then update for the mask coverage

        :param l1:
        :return: None
        """

        # Get the land/ocean flag from the grid
        land_ocean_flag, distance_to_coast = self.get_trajectory(
            l1.time_orbit.longitude,
            l1.time_orbit.latitude,
            l1.time_orbit.timestamp[0].month
        )

        # TODO: Update when global land/ocean flag becomes available
        # NOTE: Currently only a grid for the northern hemisphere is available, thus
        #       there will be regions not covered by the gridded land/ocean flag.
        #       That's why the built/in land/ocean flag will only be updated for
        #       region with coverage and be left untouched everywhere else.
        #       The arrays are added to the l1 classifier container regardless
        #       for consistency.

        # --- Update the L1 data container ---
        # 1. Save both extracted variables to classifier data groups
        l1.classifier.add(land_ocean_flag, "hr_land_ocean_flag")
        l1.classifier.add(distance_to_coast, "distance_to_coast")

        # 2. Save original ESA surface type variable to classifier data group
        #    and update the surface type flag in the surface type data group
        l1.classifier.add(l1.surface_type.flag, "orig_land_ocean_flag")

        # 3. Update the surface type instance
        valid_mask_indices = land_ocean_flag != self.cfg.get("dummy_val")["land_ocean_flag"]
        flag_update = np.full(l1.n_records, SURFACE_TYPE_DICT["invalid"])
        flag_update[land_ocean_flag == 1] = SURFACE_TYPE_DICT["land"]
        flag_update[land_ocean_flag == 0] = SURFACE_TYPE_DICT["ocean"]
        updated_surface_type_flag = l1.surface_type.flag.copy()
        updated_surface_type_flag[valid_mask_indices] = flag_update[valid_mask_indices]
        l1.surface_type.set_flag(updated_surface_type_flag)

    def get_trajectory(
            self,
            longitude: npt.NDArray,
            latitude: npt.NDArray,
            month_number: int,
    ) -> Tuple[npt.NDArray, npt.NDArray]:
        """
        Get extract of the land/ocean flag and the distance to coast value for an
        array of (longitude, latitude) positions.

        If longitude, latitude is outside the grid, a pre-defined dummy value will
        be returned.

        :param longitude: Longitude values in degrees
        :param latitude: latitude values in degrees

        :raises None:

        :return: land ocean flag & distance to coast values for longitude, latitude positions
        """

        hemisphere = "nh" if np.mean(latitude) > 0 else "sh"
        grid2track = GridTrajectoryExtract(longitude, latitude, self.cfg.grid_def[hemisphere])

        var = self.data[hemisphere].ds.surface_classification_flag.sel(month=month_number).values
        land_ocean_flag = grid2track.get_from_grid_variable(
            var,
            outside_value=self.dummy_val["surface_classification_flag"]
        )

        distance_to_coast = grid2track.get_from_grid_variable(
            self.distance_to_coast_grid,
            outside_value=self.dummy_val["distance_to_coast"]
        )
        return land_ocean_flag, distance_to_coast