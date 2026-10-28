# -*- coding: utf-8 -*-
"""

"""

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"

from pathlib import Path
import numpy as np
import xarray as xr
import numpy.typing as npt
from loguru import logger
from typing import List, Tuple, Union, Optional, Literal

from pydantic import BaseModel, computed_field
from scipy import ndimage

from pysiral import psrlcfg

from pysiral.core.legacy_classes import AttrDict
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

    @computed_field
    def dimension(self) -> dict[str, Union[int, float]]:
        return {"n_cols": 43200, "n_lines": 43200, "dx": 250, "dy": 250}

    def to_attrdict(self) -> AttrDict:
        grid_def_dict = self.model_dump()
        grid_def_dict.pop("center_latitude")
        return AttrDict(**grid_def_dict)


class SurfaceClassificationMaskConfig(BaseModel):
    local_machine_def_auxclass: str = "mask"
    local_machine_def_tag: str = "cryotempo_surface_classification"
    version: str = "v1p0"
    dummy_val: dict[str, int] = {"land_ocean_flag": -1, "region_code": -1, "surface_classification": -1}
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
    """
    Lazy access to the surface classification mask of one hemisphere.

    The full mask file decompresses to ~26 GB (12 monthly 43200 x 43200 surface
    classification grids, plus the sea ice region and land/ocean grids). Only the
    variables needed for the trajectory extraction are read from file and only
    when first requested: the surface classification of the current month and the
    sea ice region (~1.9 GB each). The monthly slice is replaced when the month changes.
    """

    def __init__(self, filepath: Path, grid_def: GridDefinition) -> None:
        self.filepath = filepath
        self.grid_def = grid_def
        self._month_number: Optional[int] = None
        self._surface_classification: Optional[npt.NDArray] = None
        self._sea_ice_region: Optional[npt.NDArray] = None
        self._dtypes: dict[str, np.dtype] = {}

    @classmethod
    def from_file(cls, filepath: Path, grid_def: GridDefinition) -> "SurfaceClassificationData":
        return cls(filepath=filepath, grid_def=grid_def)

    def get_surface_classification(self, month_number: int) -> npt.NDArray:
        """
        Return the surface classification grid for the given month (1-12)

        :param month_number: Month number (1-12)

        :return: surface classification grid (yc, xc)
        """
        if month_number != self._month_number:
            # Release the previous month before reading the next one
            self._surface_classification = None
            logger.debug(f"Reading surface classification (month={month_number}) from {self.filepath}")
            with xr.open_dataset(self.filepath) as ds:
                self._surface_classification = ds.surface_classification.isel(month=month_number - 1).values
            self._month_number = month_number
        return self._surface_classification

    def get_sea_ice_region(self) -> npt.NDArray:
        """
        Return the sea ice region grid

        :return: sea ice region grid (yc, xc)
        """
        if self._sea_ice_region is None:
            logger.debug(f"Reading sea ice region from {self.filepath}")
            with xr.open_dataset(self.filepath) as ds:
                self._sea_ice_region = ds.sea_ice_region.values
        return self._sea_ice_region

    def get_dtype(self, variable_name: str) -> np.dtype:
        """
        Return the dtype of a mask variable (reads only the file metadata)

        :param variable_name: Name of the variable in the mask file

        :return: dtype of the variable
        """
        if variable_name not in self._dtypes:
            with xr.open_dataset(self.filepath) as ds:
                self._dtypes[variable_name] = ds[variable_name].dtype
        return self._dtypes[variable_name]


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
        # NOTE: No mask data is read here. The masks are read on first use and only for the
        #       target hemisphere(s) of the pre-processor (see `set_target_hemisphere`).
        self.target_hemisphere = ["nh", "sh"]
        self.data: dict[str, SurfaceClassificationData] = {}
        for hemisphere in ["nh", "sh"]:
            mask_filepath = self.cfg.lookup_directory / str(self.cfg.filename_template).format(
                hemisphere=hemisphere, version=self.cfg.version
            )
            self.data[hemisphere] = SurfaceClassificationData.from_file(mask_filepath, self.cfg.grid_def[hemisphere])

    def set_target_hemisphere(self, target_hemisphere: List[str]) -> None:
        """
        Restrict the mask lookup to the target hemisphere(s) of the pre-processor.
        Records in other hemispheres are treated as outside the mask, and the mask
        file of a non-target hemisphere is never read.

        :param target_hemisphere: List of hemisphere names (`north`, `south`)
        """
        names = {"north": "nh", "south": "sh"}
        self.target_hemisphere = [names[hemisphere] for hemisphere in target_hemisphere]
        logger.debug(f"Surface classification mask restricted to hemisphere(s): {self.target_hemisphere}")

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
        surface_classification, region_code = self.get_trajectory(
            l1.time_orbit.longitude,
            l1.time_orbit.latitude,
            l1.time_orbit.timestamp[0].month
        )

        # --- Update the L1 data container ---
        # 1. Save both extracted variables to classifier data groups
        l1.classifier.add(surface_classification, "surface_classification")
        l1.classifier.add(region_code, "region_code")

        # 2. Save original ESA surface type variable to classifier data group
        #    and update the surface type flag in the surface type data group
        l1.classifier.add(l1.surface_type.flag, "orig_land_ocean_flag")

        # 3. Fill the surface type instance with the new values
        valid_mask_indices = surface_classification != self.cfg.dummy_val.get("surface_classification")
        flag_update = np.full(l1.n_records, SURFACE_TYPE_DICT["invalid"])

        # TODO: harc-coded flags
        # The `land` surface type flag is either 2 (continents) or 3 (islands) in the surface classification
        is_land = np.isin(surface_classification, [2, 3])  # continents & islands
        flag_update[is_land] = SURFACE_TYPE_DICT["land"]

        # The `ocean` surface type flag is either 0 (ocean) or 1 (sea ice climatology) in the surface classification
        is_ocean = np.isin(surface_classification, [0, 1])
        flag_update[is_ocean] = SURFACE_TYPE_DICT["ocean"]

        # The `land_ice` surface type flag is 4 in the surface classification
        flag_update[surface_classification == 4] = SURFACE_TYPE_DICT["land_ice"]

        # Now update the surface type flag in the l1 surface type instance
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
        :param month_number: Month number (1-12) for which the surface classification is requested

        :raises None:

        :return: surface classification and sea ice region id for longitude, latitude positions
        """

        # Records outside the target hemisphere(s) get the value for positions outside the grid
        data = self.data[self.target_hemisphere[0]]
        surface_classification = self._get_outside_array(
            data.get_dtype("surface_classification"), longitude.size, "surface_classification"
        )
        sea_ice_region = self._get_outside_array(
            data.get_dtype("sea_ice_region"), longitude.size, "region_code"
        )

        for hemisphere in self.target_hemisphere:

            # Each record is looked up in the mask of its own hemisphere
            is_hemisphere = latitude >= 0.0 if hemisphere == "nh" else latitude < 0.0
            if not is_hemisphere.any():
                continue
            data = self.data[hemisphere]

            grid_def = self.cfg.grid_def[hemisphere].to_attrdict()
            grid2track = GridTrajectoryExtract(longitude[is_hemisphere], latitude[is_hemisphere], grid_def)

            surface_classification[is_hemisphere] = grid2track.get_from_grid_variable(
                data.get_surface_classification(month_number),
                outside_value=self.cfg.dummy_val["surface_classification"]
            )

            sea_ice_region[is_hemisphere] = grid2track.get_from_grid_variable(
                data.get_sea_ice_region(),
                outside_value=self.cfg.dummy_val["region_code"]
            )

        return surface_classification, sea_ice_region

    def _get_outside_array(self, dtype: np.dtype, n_records: int, dummy_val_name: str) -> npt.NDArray:
        """
        Return an array for all records, filled with the value that the grid extraction
        assigns to positions outside the grid (the dummy value cast to the grid dtype
        by `scipy.ndimage.map_coordinates`).

        :param dtype: dtype of the grid variable
        :param n_records: Number of records
        :param dummy_val_name: Name of the dummy value in the configuration

        :return: Array of size n_records
        """
        outside_value = ndimage.map_coordinates(
            np.zeros((1, 1), dtype=dtype),
            [[-1.0], [-1.0]],
            order=0,
            mode="constant",
            cval=self.cfg.dummy_val[dummy_val_name]
        )[0]
        return np.full(n_records, outside_value, dtype=dtype)
