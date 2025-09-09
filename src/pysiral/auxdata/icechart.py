# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 13:57:56 2016

@author: Stefan

Module created for FMI version of pysiral
"""

__all__ = ["IC", "ICA"]


import contextlib
import copy
import datetime
from dataclasses import dataclass, field, InitVar
from datetime import date, timedelta
from pathlib import Path
from typing import Dict, Optional, Tuple, Literal, List
from itertools import product
from loguru import logger

import pandas as pd
import geopandas as gpd
import numpy as np
import pyproj
import xarray as xr
from parse import parse
from PIL import Image
from pyproj import Proj
from returns.result import Failure, Success, Result
from shapely import LineString, MultiPoint
from shapely.strtree import STRtree

from pysiral.auxdata import AuxdataBaseClass
from pysiral.l2data import Level2Data
from pysiral.core.iotools import ReadNC
from pysiral.auxdata import GridTrackInterpol

# Sea Ice Concentration (SIC) code to class conversion lookup table.
SIC_LOOKUP = {
    'polygon_idx': 0,  # Index of polygon number.
    'total_sic_idx': 1,  # Total Sea Ice Concentration Index, CT.
    'sic_partial_idx': [2, 5, 8],  # Partial SIC polygon code index. CA, CB, CC.
    0: 0,
    1: 0,
    2: 0,
    55: 0,
    10: 1,  # 10 %
    20: 2,  # 20 %
    30: 3,  # 30 %
    40: 4,  # 40 %
    50: 5,  # 50 %
    60: 6,  # 60 %
    70: 7,  # 70 %
    80: 8,  # 80 %
    90: 9,  # 90 %
    91: 10,  # 100 %
    92: 10,  # Fast ice

    89: 8.5,  # 8/10 – 9/10
    81: 9,    # 8/10 – 10/10
    79: 8,    # 7/10 – 9/10
    78: 7.5,  # 7/10 – 8 /10
    68: 7,    # 6/10 – 8/10
    67: 6.5,  # 6/10 – 7/10
    57: 6,    # 5/10 – 7/10
    56: 5.5,  # 5/10 – 6/10
    46: 5,    # 4/10 – 6/10
    45: 4.5,  # 4/10 – 5/10
    35: 4,    # 3/10 – 5/10
    34: 3.5,  # 3/10 – 4/10
    24: 3,    # 2/10 – 4/10
    23: 2.5,  # 2/10 – 3/10
    13: 2,    # 1/10 – 3/10
    12: 1.5,  # 1/10 – 2/10

    'mask': 255,
    'n_classes': 12
}

# Names of the SIC classes.
SIC_GROUPS = {
    0: 0,
    1: 10,
    2: 20,
    3: 30,
    4: 40,
    5: 50,
    6: 60,
    7: 70,
    8: 80,
    9: 90,
    10: 100
}


# # Stage of Development code to class conversion lookup table.
# SOD_LOOKUP = {
#     'sod_partial_idx': [3, 6, 9],  # Partial SIC polygon code index. SA, SB, SC.
#     'threshold': 0.65,  # < 1. Minimum partial percentage SIC of total SIC to select SOD. Otherwise ambiguous polygon.
#                        # larger than threshold.
#     'invalid': -9,  # Value for polygons where the SOD is ambiguous or not filled.
#     'water': 0,
#     0: 0,
#     80: 0,  # No stage of development
#     81: 1,  # New ice
#     82: 2,  # Nilas, ring ice
#     83: 3,  # Young ice
#     84: 4,  # Grey ice
#     85: 5,  # White ice
#     86: 6,  # First-year ice, overall categary
#     87: 6,  # Thin first-year ice
#     88: 6,  # Thin first-year ice, stage 1
#     89: 6,  # Thin first-year ice, stage 2
#     91: 6,  # Medium first-year ice
#     93: 6,  # Thick first-year ice
#     95: 7,  # Old ice
#     96: 7,  # Second year ice
#     97: 7,  # Multi-year ice
#     98: 255,  # Glacier ice
#     99: 255,
#     'mask': 255,
#     'n_classes': 9
# }

# Stage of Development code to class conversion lookup table.
SOD_LOOKUP = {
    'sod_partial_idx': [3, 6, 9],  # Partial SIC polygon code index. SA, SB, SC.
    'threshold': 0.65,  # < 1. Minimum partial percentage SIC of total SIC to select SOD. Otherwise, ambiguous polygon.
                        # larger than threshold.
    'invalid': -9,  # Value for polygons where the SOD is ambiguous or not filled.
    'water': 0,
    0: 0,
    80: 0,  # No stage of development
    81: 1,  # New ice
    82: 2,  # Nilas, ring ice
    83: 3,  # Young ice
    84: 4,  # Grey ice
    85: 5,  # White ice
    86: 6,  # First-year ice, overall category
    87: 7,  # Thin first-year ice
    88: 8,  # Thin first-year ice, stage 1
    89: 9,  # Thin first-year ice, stage 2
    91: 10,  # Medium first-year ice
    93: 11,  # Thick first-year ice
    95: 12,  # Old ice
    96: 13,  # Second year ice
    97: 14,  # Multi-year ice
    98: 255,  # Glacier ice
    99: 255,
    'mask': 255,
    'n_classes': 16
}

# Names of the SOD classes.
SOD_GROUPS = {

    # 0: 'Open water',
    # 1: 'New Ice',
    # 2: 'Nilas, ring ice',
    # 3: 'Young ice',
    # 4: 'Grey ice',
    # 5: 'White ice',
    # 6: 'First-year ice',
    # 7: 'Multi-year ice'

    0: 'Open water',
    1: 'New Ice',
    2: 'Nilas, ring ice',
    3: 'Young ice',
    4: 'Grey ice',  # Grey ice
    5: 'White ice',  # White ice
    6: 'First-year ice, overall categary',
    7: 'Thin first-year ice',
    8: 'Thin first-year ice, stage 1',
    9: 'Thin first-year ice, stage 2',
    10: 'Medium first-year ice',
    11: 'Thick first-year ice',
    12: 'Old ice',
    13: 'Second year ice',
    14: 'Multi-year ice'
}


# Ice floe/form code to class conversion lookup table.
FLOE_LOOKUP = {
    'floe_partial_idx': [4, 7, 10],  # Partial SIC polygon code index. FA, FB, FC.
    'threshold': 0.65,  # < 1. Minimum partial concentration to select floe. Otherwise, polygon may be ambiguous.
    'invalid': -9,  # Value for polygons where the floe is ambiguous or not filled.
    'water': 0,
    0: 0,
    22: 7,  # Pancake ice
    1: 8,  # Shuga / small ice cake
    2: 1,  # Ice cake
    3: 2,  # Small floe
    4: 3,  # Medium floe
    5: 4,  # Big floe
    6: 5,  # Vast fæpe
    7: 9,  # Gian floe
    8: 255,  # Fast ice
    9: 6,  # Growlers, floebergs or floebits
    10: 6,  # Icebergs
    21: 255,  # Level ice
    'fastice_class': 255,
    'mask': 255,
    'n_classes': 11
}

# Names of the FLOE classes.
FLOE_GROUPS = {
    0: 'Open water',
    1: 'Cake Ice',
    2: 'Small floe',
    3: 'Medium floe',
    4: 'Big floe',
    5: 'Vast floe',
    6: 'Bergs',
    7: 'Pancake ice',
    8: 'Shuga / small ice cake',
    9: 'Gian floe'
}

CHARTS = ['SIC', 'SOD', 'FLOE']

LOOKUP_NAMES = {
    'SIC': SIC_LOOKUP,
    'SOD': SOD_LOOKUP,
    'FLOE': FLOE_LOOKUP
}


ICECHART_NOT_FILLED_VALUE = -9  # Used in converting polygon codes into the preprocessed ice charts.
ICECHART_UNKNOWN = 99  # Used in converting polygon codes into the preprocessed ice charts.


@dataclass
class NSIDCIceChartFileEntry:
    file_path: Path
    region: str
    issue_date: datetime.date
    validity_start_date: datetime.date
    validity_end_date: datetime.date

    @property
    def id(self) -> str:
        return self.file_path.name


class NSIDCIceChartFileCatalog(object):

    def __init__(self, lookup_directory: Path, file_validity_period_days_default: int = 7) -> None:
        self.lookup_directory = lookup_directory
        self.file_validity_period_days_default = file_validity_period_days_default
        self.filepaths = sorted(list(self.lookup_directory.rglob("*.shp")))
        self._ctlg = self._catalog_files()

    def _catalog_files(self) -> Dict[int, NSIDCIceChartFileEntry]:
        """
        Create a catolog of ice charts files
        TODO: So far it is assumed that ice charts are valid 7 days with issue date (filename) as last date
        :return:
        """

        icechart_ctlg = {}
        for filepath in self.filepaths:

            filename_attributes = parse(self.filename_parser, filepath.name)
            year_offset = 2000 if filename_attributes.named["year"] < 30 else 1900
            year = year_offset + filename_attributes.named["year"]

            issue_date = date(year, filename_attributes.named["month"], filename_attributes.named["day"])
            validity_start_date = issue_date - timedelta(days=self.file_validity_period_days_default - 1)
            validity_end_date = issue_date

            # Using integer (days since 1/1/1) as key for fast lookup
            icechart_ctlg[issue_date.toordinal()] = NSIDCIceChartFileEntry(
                filepath,
                filename_attributes.named["region_code"],
                issue_date,
                validity_start_date,
                validity_end_date
            )
        return icechart_ctlg

    def query_icechart_files(self, target_date: datetime.date) -> Optional[NSIDCIceChartFileEntry]:
        """
        Select icechart file for specified date

        :param target_date: Target date

        :raises ValueError: Multiple ice chart files found

        :return: ice chart catalog entry or None
        """

        # Get matches by checking if the target date number is within ice chart validity range
        target_date_ordinal = target_date.toordinal()
        condition = np.logical_and(
            target_date_ordinal <= self.issue_dates_ordinal,
            target_date_ordinal >= self.issue_dates_ordinal - self.file_validity_period_days_default - 1
        )
        try:
            match_index = np.argwhere(condition)[0]
        except IndexError:
            return None

        # No ice chart found for target date
        if match_index.size == 0:
            return None

        # If not 0, there should only be one match
        if match_index.size > 1:
            raise ValueError(f"Found more than one ice chart for {target_date}.")

        return self._ctlg[int(self.issue_dates_ordinal[match_index[0]])]

    @property
    def filename_parser(self) -> str:
        return "{region_code:D}{year:2d}{month:2d}{day:2d}.shp"

    @property
    def issue_dates_ordinal(self) -> np.ndarray:
        return np.array(sorted(list(self._ctlg.keys())))


class NSIDCSeaIceChartsSIGRID3(AuxdataBaseClass):

    def __init__(self, *args, **kwargs):
        super(NSIDCSeaIceChartsSIGRID3, self).__init__(*args, **kwargs)
        self.ctlg = self.get_file_catalog()
        self.loaded_file = None
        self.ice_chart = None
        self.ice_chart_bounds = None

    def get_l2_track_vars(self, l2: Level2Data) -> None:
        """
        Add ice chart parameters to Level-2 data object

        :param l2: Level-2 data object

        :return: None, l2 is changed in place
        """

        # Get the target file
        self.set_requested_date_from_l2(l2)
        icechart_file = self.ctlg.query_icechart_files(self.requested_date)

        # No action if no corresponding ice chart file is found
        if icechart_file is None:
            return

        # Load icechart only if required
        if self.loaded_file != icechart_file.id:
            self.load_ice_chart_file(icechart_file)

        # Get icechart data frame for trajectory
        ice_chart_l2_track = self.extract_track(l2.longitude, l2.latitude)

        # Pre-process parameters and set to l2 object
        self._set_l2_parameters(l2, ice_chart_l2_track)

    def get_file_catalog(self) -> NSIDCIceChartFileCatalog:
        lookup_directory = Path(self.cfg.local_repository) / self.cfg.options.hemisphere
        return NSIDCIceChartFileCatalog(lookup_directory)

    def load_ice_chart_file(self, icechart_file: NSIDCIceChartFileEntry) -> None:
        """
        Load and convert the ice chat data and extent

        :param icechart_file: ice chart catalog entry

        :return:
        """
        icechart_raw = gpd.read_file(icechart_file.file_path)
        self.ice_chart_bounds = icechart_raw.total_bounds
        self.ice_chart = convert_polygon_icechart(icechart_raw)

    def extract_track(self, longitude: np.ndarray, latitude: np.ndarray) -> xr.Dataset:
        """
        Extract all relevant ice chart variables along the Level-2 track

        :param longitude: Level-2 geodetic longitude
        :param latitude: Level-2 geodetic latitude

        :return: Dataset with ice chart variables for each Level-2 track record
        """

        # Dataset with ice chart variables
        ice_chart_dataset = self.get_ice_chart_dataset(longitude.size)

        # Get the Level-2 track as geometry in ice chart projection coordinates
        track_linestring, track_multipoint = self.get_l2_track_projected_geometry(longitude, latitude)

        # Get the index list of polygons that overlap with the track
        poly_overlap_idxs = self.get_polygons_overlapping_with_track(track_linestring)
        if poly_overlap_idxs.size == 0:
            return ice_chart_dataset

        # Get the polygon index for each track record
        track_poly_idx = self.get_polygon_index_per_record(track_multipoint, poly_overlap_idxs)

        # Map ice chart data to track
        for index in np.arange(longitude.size):
            polygon_index = track_poly_idx[index]
            ice_chart_dataset.polygon_index[index] = polygon_index

            ice_chart_dataset.SIC_T[index] = self.ice_chart.iloc[polygon_index]['SIC_T']

            ice_chart_dataset.SIC_A[index] = self.ice_chart.iloc[polygon_index]['SIC_A']
            ice_chart_dataset.SOD_A[index] = self.ice_chart.iloc[polygon_index]['SOD_A']
            ice_chart_dataset.Floe_A[index] = self.ice_chart.iloc[polygon_index]['Floe_A']

            ice_chart_dataset.SIC_B[index] = self.ice_chart.iloc[polygon_index]['SIC_B']
            ice_chart_dataset.SOD_B[index] = self.ice_chart.iloc[polygon_index]['SOD_B']
            ice_chart_dataset.Floe_B[index] = self.ice_chart.iloc[polygon_index]['Floe_B']

            ice_chart_dataset.SIC_C[index] = self.ice_chart.iloc[polygon_index]['SIC_C']
            ice_chart_dataset.SOD_C[index] = self.ice_chart.iloc[polygon_index]['SOD_C']
            ice_chart_dataset.Floe_C[index] = self.ice_chart.iloc[polygon_index]['Floe_C']

        return ice_chart_dataset

    @staticmethod
    def get_ice_chart_dataset(dim_size: int, dim_name: str = "time") -> xr.Dataset:
        """
        Datamodel of Level-2 track data extracted from the ice chart

        :param dim_size: Number of Level-2 data records
        :param dim_name: The name of the dimension

        :return: xarray.Dataset
        """
        return xr.Dataset({
            "SIC_T": ([dim_name], np.full((dim_size,), np.nan)),
            'SIC_A': ([dim_name], np.full((dim_size,), np.nan)),
            'SOD_A': ([dim_name], np.full((dim_size,), -9999)),
            'Floe_A': ([dim_name], np.full((dim_size,), -9999)),
            'SIC_B': ([dim_name], np.full((dim_size,), np.nan)),
            'SOD_B': ([dim_name], np.full((dim_size,), -9999)),
            'Floe_B': ([dim_name], np.full((dim_size,), -9999)),
            'SIC_C': ([dim_name], np.full((dim_size,), np.nan)),
            'SOD_C': ([dim_name], np.full((dim_size,), -9999)),
            'Floe_C': ([dim_name], np.full((dim_size,), -9999)),
            'polygon_index': ([dim_name], np.full((dim_size,), -9999))
        },
            coords={dim_name: ([dim_name], np.arange(dim_size))},
        )

    def get_l2_track_projected_geometry(
            self,
            longitude: np.ndarray,
            latitude: np.ndarray
    ) -> Tuple[LineString, MultiPoint]:
        """
        Converts the geodetic longitude, latitude values of the Level-2 track in
        shapely geometries in the Ice Chart coordinate reference system.

        The linestring will be used for fast icechart polygon subsetting
        and the MultiPoint representation for polygon index lookup.

        :param longitude: Level-2 geodetic longitude
        :param latitude: Level-2 geodetic latitude

        :return: Linestring and Multipoint objects.
        """

        proj = Proj(self.ice_chart.crs)
        xc, yc = proj(longitude, latitude)
        track_linestring = LineString(list(zip(xc, yc)))
        track_multipoint = MultiPoint(track_linestring.coords)
        return track_linestring, track_multipoint

    def get_polygons_overlapping_with_track(self, track_linestring: LineString) -> np.ndarray:
        """
        Get the index of polygons overlapping with the track.

        :param track_linestring: Level-2 track geometry in ice chart coordinate reference system

        :return: List of ice chart polygons indices overlapping with track
        """
        ice_chart_polygon_strtree = STRtree(self.ice_chart.geometry)
        return ice_chart_polygon_strtree.query(track_linestring)

    def get_polygon_index_per_record(
            self,
            track_multipoint,
            poly_overlap_idxs,
            missing_value: int = -1
    ) -> np.ndarray:
        """
        Identify the corresponding polygon index for each data record. This is done with by
        querying the nearest geometry with a strict maximum distance threshold.

        :param track_multipoint: Level-2 track geometry in ice chart coordinate reference system
        :param poly_overlap_idxs: List of ice chart polygons indices overlapping with track
        :param missing_value: Default value when no overlapping polygon is found

        :return: Index of overlapping polygon (or missing value) for each Level-2 record
        """
        ice_chart_polygon_subset_strtree = STRtree(self.ice_chart.geometry[poly_overlap_idxs])
        track_poly_idx = np.full(len(track_multipoint.geoms), missing_value)
        for idx, geom in enumerate(track_multipoint.geoms):
            sub_idx = ice_chart_polygon_subset_strtree.query_nearest(geom, max_distance=1.0)
            # in case of missing polygon, statement below will raise an ValueError
            # -> ignoring exception will lead to missing value in data record (as intended)
            with contextlib.suppress(ValueError):
                track_poly_idx[idx] = poly_overlap_idxs[sub_idx]
        return track_poly_idx

    def _set_l2_parameters(self, l2: Level2Data, ice_chart_l2_track: xr.Dataset) -> None:
        """
        Set the parameter sea ice concentrations, stage of developement and floe
        for the three categories A, B, C as multidim parameters and the total
        concentration as single parameter
        """

        # Only total sea ice concentration is single dimension parameter
        self.register_auxvar(
            "ictsic", "ice_chart_sea_ice_concentration_total",
            ice_chart_l2_track.SIC_T.values
        )

        dims = {"new_dims": (("ice_chart_class", 3),),
                "dimensions": ("time", "ice_chart_class"),
                "add_dims": (("ice_chart_class", np.arange(3)),)}

        classes_sic = np.column_stack([
            ice_chart_l2_track.SIC_A.values,
            ice_chart_l2_track.SIC_B.values,
            ice_chart_l2_track.SIC_C.values,
        ])
        l2.set_multidim_auxiliary_parameter(
            "icsicc", "ice_chart_sea_ice_concentration_classes",
            classes_sic, dims, update=True
        )
        classes_sod = np.column_stack([
            ice_chart_l2_track.SOD_A.values,
            ice_chart_l2_track.SOD_B.values,
            ice_chart_l2_track.SOD_C.values,
        ])
        l2.set_multidim_auxiliary_parameter(
            "icsodc", "ice_chart_stage_of_development_classes",
            classes_sod, dims, update=True
        )
        classes_floe = np.column_stack([
            ice_chart_l2_track.Floe_A.values,
            ice_chart_l2_track.Floe_B.values,
            ice_chart_l2_track.Floe_C.values,
        ])
        l2.set_multidim_auxiliary_parameter(
            "icfloec", "ice_chart_floe_parameter_classes",
            classes_floe, dims, update=True
        )


@dataclass
class USNICGridFileEntry:
    """
    A single US National Ice Center sea ice chart grid file and its properties
    derived from the filename.
    """
    file_path: Path | None
    grid: str = field(init=False)
    region: str = field(init=False)
    version: str = field(init=False)
    validity_start_date: datetime.date = field(init=False)
    validity_end_date: datetime.date = field(init=False)
    issue_date: datetime.date = field(init=False)
    filename_parser: InitVar[str] = "usnic-icechart-{grid}-{start_date}_{end_date}-{version}.nc"

    def __post_init__(self, filename_parser: str) -> None:
        """
        Compute file properties from filename. The file path must match the filename
        parser, otherwise the file is considered invalid and the file path is set to None.

        :param filename_parser: InitVar object containing the filename parser

        :return: None
        """

        # Parse filename and handle invalid files
        filename_attributes = parse(filename_parser, self.file_path.name)
        if not hasattr(filename_attributes, "named"):
            self.file_path = None
            return

        self.grid = filename_attributes.named["grid"]
        self.region = filename_attributes.named["grid"].split("_")[0]
        self.version = filename_attributes.named["version"]
        self.validity_start_date = date.fromisoformat(filename_attributes.named["start_date"])
        self.validity_end_date = date.fromisoformat(filename_attributes.named["end_date"])
        self.issue_date = self.validity_end_date - timedelta(days=1)

    @property
    def is_valid(self) -> bool:
        return self.file_path is not None


class USNICGridFileCatalog(object):
    """
    A catalog of USNIC grid files in a given directory and
    the ability to query the closest file for a given date.
    """

    def __init__(self, lookup_directory: Path) -> None:
        """
        Initialize the catalog by scanning the given directory for files.

        :param lookup_directory: The lookup directory. May be root directory
            or hemisphere subdirectory (e.g. 'north' or 'south').
        """
        self.lookup_directory = lookup_directory
        self.ctlg = self._catalog_files()

    def get_closest(
            self,
            target_date_list: List[int],
            max_offset_days: int = 7,
            tie_breaker: Literal["before", "after"] = "before"
    ) -> Path | None:
        """
        Get the closest file for the given date within the maximum offset.

        :param target_date: The target date.
        :param max_offset_days: The maximum offset in days an ice chart period can
            be away from the target date. The default is 7 days, to allow for the
            periods where ice charts were only distributed bi-weekly.
        :param tie_breaker: If two files are equally close, choose either the one
            before or after the target date.

        :return: Filepath of the target file, or None if no file is found within
        """
        target_date = date(*target_date_list)
        result = self._get_closest(target_date, max_offset_days, tie_breaker)
        match result:
            case Success(value):
                logger.info(f"Found {value} on {target_date}.")
                return value[0]
            case Failure(_):
                return None
        return None

    def _get_closest(
            self,
            target_date: date,
            max_offset_days: int,
            tie_breaker: Literal["before", "after"]
    ) -> Result[Tuple[Path, int], Tuple[str, int]]:
        """
        Get the closest file for the given date within the maximum offset.

        :param target_date: The target date.
        :param max_offset_days: The maximum offset in days an ice chart period can
            be away from the target date. The default is 7 days, to allow for the
            periods where ice charts were only distributed bi-weekly.
        :param tie_breaker: If two files are equally close, choose either the one
            before or after the target date.

        :return: A tuple of the file path and the offset in days, or (None, None)
            if no file is found within the maximum offset.
        """
        if self.ctlg.empty:
            return Failure(("Catalog is empty", -1))

        # Check if there is any overlap with the 7-day validity periods of ice charts
        is_in_period = (
            (self.ctlg["validity_start_date"] <= target_date) &
            (self.ctlg["validity_end_date"] >= target_date)
        )
        if is_in_period.any():
            matched_entry = self.ctlg[is_in_period].iloc[0]
            return Success((matched_entry.file_path, 0))

        # Find the closest ice chart (either start or end)
        closest_before = self._get_days_offset(self.ctlg["validity_start_date"], target_date)
        closest_after =  self._get_days_offset(self.ctlg["validity_end_date"], target_date)

        closest_before_idx = np.argmin(closest_before)
        closest_before_value = closest_before[closest_before_idx]
        closest_after_idx = np.argmin(closest_after)
        closest_after_value = closest_after[closest_after_idx]

        if closest_before_value == closest_after_value:
            if tie_breaker == "before":
                matched_entry = self.ctlg.iloc[closest_before_idx]
                offset_days = closest_before_value
            else:
                matched_entry = self.ctlg.iloc[closest_after_idx]
                offset_days = closest_after_value
        elif closest_before_value < closest_after_value:
            matched_entry = self.ctlg.iloc[closest_before_idx]
            offset_days = closest_before_value
        elif closest_after_value < closest_before_value:
            matched_entry = self.ctlg.iloc[closest_after_idx]
            offset_days = closest_after_value
        else:
            raise ValueError("should not happen")

        return (
            Success((matched_entry.file_path, int(offset_days))) if offset_days <= max_offset_days
            else Failure((f"Closest match ({int(offset_days)}) exceeding {max_offset_days=}", int(offset_days)))
        )

    @staticmethod
    def _get_days_offset(date_series: pd.Series, target_date: date) -> np.ndarray:
        """
        Get the offset in days between the target date and the issue date of each file.

        :param target_date: The target date.

        :return: An array of offsets in days.
        """
        return np.abs(np.fromiter([d.days for d in (date_series - target_date).values], int))

    def _catalog_files(self) -> pd.DataFrame:
        """
        Create catalog of file properties as a pandas DataFrame

        :return: Dataframe with file catalog entries
        """
        filepaths = sorted(list(self.lookup_directory.rglob("*.nc")))
        entries = [USNICGridFileEntry(fp) for fp in filepaths]
        entries = [e for e in entries if e.is_valid]
        return pd.DataFrame(entries)


class USNICGrid(AuxdataBaseClass):
    """
    Class for US National Ice Center sea ice charts in gridded netCDF format.
    Source data is available from the NSIDC (https://nsidc.org/data/g02135)
    Software to convert into gridded format: https://gitlab.awi.de/siteo/cryotempo/icgrid
    """

    def __init__(self, *args, **kwargs) -> None:
        """
        Init the class.
        NOTE: The options template can be different for continuous data sets (is_cdr_icdr: False)
              and those with a dedicated split into a climate data record (cdr) and an interim
              climate data record (icdr). A pre-processing of the options dictionary is therefore necessary
              to follow the mechanics of the auxiliary data class.
        :param args:
        :param kwargs:
        """
        super(USNICGrid, self).__init__(*args, **kwargs)
        self.ctlg = self.get_file_catalog()
        self._data = None
        self.start_time = None
        self.hemisphere_code = None

    def get_file_catalog(self) -> USNICGridFileCatalog:
        """
        Generate the file catalog object for the specific hemisphere

        :return: USNICGridFileCatalog object
        """
        lookup_directory = Path(self.cfg.local_repository) / self.cfg.options.hemisphere
        return USNICGridFileCatalog(lookup_directory)

    def get_l2_track_vars(self, l2):
        """
        Mandadory method of AuxdataBaseClass subclass.
        Registers two variables to the Level-2 data container:
            - MYI fraction (id: sitype, name: sea_ice_type)
            - MYI fraction uncertainty (id: sitype.uncertainty, name: sea_ice_type_uncertainty)
        :param l2: Level-2 Data object
        :return: None
        """

        # Set the requested data
        self.set_requested_date_from_l2(l2)

        # Update the external data
        self.update_external_data()

        # Check if error with file I/O
        dataset = self.get_empty_dataset(l2.time)
        if self.error.status or self._data is None:
            ice_chart_l2_track = dataset.copy()
        else:
            # Get icechart data frame for trajectory
            ice_chart_l2_track = self.extract_track(l2.longitude, l2.latitude, dataset)

        # Pre-process parameters and set to l2 object
        self.set_l2_parameters(l2, ice_chart_l2_track)

    @staticmethod
    def get_empty_dataset(time: np.ndarray) -> xr.Dataset:
        """
        Create a data structure with the correct dimensions and fill values
        for the ice chart parameters. This dataset is used if to map
        the ice charts parameters onto, as well as dummy dataset if no
        ice chart file is found.

        :param time:

        :return: Empty trajectory ice chart dataset
        """
        time_dims = time.shape
        class_dims = (3,) + time_dims

        # Rer
        attrs = {
            "time_dim_parameter": [
                "sea_ice_concentration_total",
                "stage_of_development_highest_concentration",
                "fraction_thin_ice",
                "fraction_first_year_ice",
                "fraction_multi_year_ice",
                "number_of_seaice_classes",
            ],
            "class_time_dim_parameter": [
                "sea_ice_concentration_partial",
                "stage_of_development_partial",
                "stage_of_development_partial_is_overall_class",
                "form_of_ice_partial",
            ],
        }
        coords = {
            "time": xr.Variable(("time",), time),
            "n_classes" : xr.Variable(("n_classes",), np.arange(3))
        }
        data_vars = {
            "sea_ice_concentration_total": (("time",), np.full(time_dims, np.nan)),
            "stage_of_development_highest_concentration": (("time",), np.full(time_dims, -1)),
            "fraction_thin_ice": (("time",), np.full(time_dims, np.nan)),
            "fraction_first_year_ice": (("time",), np.full(time_dims, np.nan)),
            "fraction_multi_year_ice": (("time",), np.full(time_dims, np.nan)),
            "number_of_seaice_classes": (("time", ), np.full(time_dims, -1)),
            'sea_ice_concentration_partial': (("n_classes", "time",), np.full(class_dims, np.nan)),
            'stage_of_development_partial': (("n_classes", "time",), np.full(class_dims, 0)),
            "stage_of_development_partial_is_overall_class": (("n_classes", "time",), np.full(class_dims, -1)),
            'form_of_ice_partial': (("n_classes", "time",), np.full(class_dims, 0)),
        }
        return xr.Dataset(attrs=attrs, coords=coords, data_vars=data_vars)

    def load_requested_auxdata(self) -> None:
        """
        Loads file from local repository only if needed
        :return:
        """

        # Retrieve the file path for the requested date from a property of the auxdata parent class
        path = Path(self.requested_filepath)

        # Validation
        if not path.is_file():
            msg = f"{self.__class__.__name__}: File not found: {path} "
            self.add_handler_message(msg)
            self.error.add_error("auxdata_missing_sitype", msg)
            return

        # Read and prepare input data
        self._data = ReadNC(path)

    def extract_track(self, longitude: float, latitude: float, data: xr.Dataset) -> xr.Dataset:
        """
        Map the ice chart parameters onto the Level-2 track

        :param longitude:
        :param latitude:

        :param data: Data structure with same variables and dimensions as the ice chart data

        :return: Filled data structure
        """

        # Set up the grid interpolator
        griddef = self.cfg.options[self.cfg.options.hemisphere]
        grid_lons, grid_lats = self._data.lon, self._data.lat
        grid2track = GridTrackInterpol(longitude, latitude, grid_lons, grid_lats, griddef)

        # First get the time-dim parameter
        for var_name in data.attrs["time_dim_parameter"]:
            data[var_name].values = grid2track.get_from_grid_variable(
                getattr(self._data, var_name)[0, :, :], flipud=True
            )

        for idx, var_name in product(np.arange(3), data.attrs["class_time_dim_parameter"]):
            data[var_name].values = grid2track.get_from_grid_variable(
                getattr(self._data, var_name)[0, idx, :, :], flipud=True
            )

        return data

    def set_l2_parameters(self, l2: Level2Data, ice_chart_l2_track: xr.Dataset) -> None:
        """
        Set the parameter sea ice concentrations, stage of developement and floe
        for the three categories A, B, C as multidim parameters and the total
        concentration as single parameter
        """

        var_id_dict = {
            "sea_ice_concentration_total": "ictsic",
            "stage_of_development_highest_concentration": "icsodhc",
            "stage_of_development_partial_is_overall_class": "icsodioc",
            "fraction_thin_ice": "icfthi",
            "fraction_first_year_ice": "icffyi",
            "fraction_multi_year_ice": "icfmyi",
            "number_of_seaice_classes": "icncls",
            "sea_ice_concentration_partial": "icsicp",
            "stage_of_development_partial": "icsodp",
            "form_of_ice_partial": "icfloep",
        }

        for var_name in ice_chart_l2_track.attrs["time_dim_parameter"]:
            self.register_auxvar(var_id_dict[var_name], var_name, ice_chart_l2_track[var_name].values)

        # # Only total sea ice concentration is single dimension parameter
        # self.register_auxvar(
        #     "ictsic", "ice_chart_sea_ice_concentration_total",
        #     ice_chart_l2_track.SIC_T.values
        # )
        #
        dims = {"new_dims": (("ice_chart_class", 3),),
                "dimensions": ("time", "ice_chart_class"),
                "add_dims": (("ice_chart_class", np.arange(3)),)}
        for name in ice_chart_l2_track.attrs["class_time_dim_parameter"]:
            l2.set_multidim_auxiliary_parameter(
                var_id_dict[name], name,
                ice_chart_l2_track[name].values, dims, update=True
            )


    @property
    def requested_filepath(self) -> Path | None:
        """
        Note: this overwrites the property in the super class due to some
        peculiarities with the filenaming (auto product changes etc)
        :return: The filepath to the target file
        """

        # The path needs to be completed if two products shall be used
        opt = self.cfg.options
        return self.ctlg.get_closest(
            self._requested_date,
            max_offset_days=opt.max_offset_days,
            tie_breaker=opt.tie_breaker
        )


class IC(AuxdataBaseClass):

    def __init__(self, *args, **kwargs):
        super(IC, self).__init__(*args, **kwargs)
        self._data_ct = None
        self._data_ca = None
        self._data_cb = None
        self._data_cc = None
        self._data_sa = None
        self._data_sb = None
        self._data_sc = None
        self._current_date = [0, 0, 0]
        self._requested_date = [-1, -1, -1]
        self.error.caller_id = self.__class__.__name__

    def get_l2_track_vars(self, l2):
        self._get_requested_date(l2)
        self._get_data(l2)
        if not self.error.status:
            ic_sics = self._get_sic_ic_track(l2)
            self.register_auxvar("ic_ct", "icechart_ct", ic_sics[0], None)
            self.register_auxvar("ic_ca", "icechart_ca", ic_sics[1], None)
            self.register_auxvar("ic_cb", "icechart_cb", ic_sics[2], None)
            self.register_auxvar("ic_cc", "icechart_cc", ic_sics[3], None)
            self.register_auxvar("ic_sa", "icechart_sa", ic_sics[4], None)
            self.register_auxvar("ic_sb", "icechart_sb", ic_sics[5], None)
            self.register_auxvar("ic_sc", "icechart_sc", ic_sics[6], None)
        else:
            self.register_auxvar("ic_ct", "icechart_ct", None, None)
            self.register_auxvar("ic_ca", "icechart_ca", None, None)
            self.register_auxvar("ic_cb", "icechart_cb", None, None)
            self.register_auxvar("ic_cc", "icechart_cc", None, None)
            self.register_auxvar("ic_sa", "icechart_sa", None, None)
            self.register_auxvar("ic_sb", "icechart_sb", None, None)
            self.register_auxvar("ic_sc", "icechart_sc", None, None)

    def _get_requested_date(self, l2):
        """ Use first timestamp as reference, date changes are ignored """
        year = l2.track.timestamp[0].year
        month = l2.track.timestamp[0].month
        day = l2.track.timestamp[0].day
        self._requested_date = [year, month, day]

    def _get_data(self, l2):
        """ Loads file from local repository only if needed """
        if self._requested_date == self._current_date:
            # Data already loaded, nothing to do
            self.add_handler_message("IC: tif already present")
            return

        paths, timedelta = self._get_local_repository_filename(l2)
        self.add_handler_message('IN ICECHART _GET_DATA, PATHS 0: %s' % paths[0])
        # print os.path.isfile(paths[0])

        # Validation
        if not Path(paths[0]).is_file():
            msg = "IC: File not found: %s " % paths[0]
            self.error.add_error("auxdata_missing_icechart", msg)
            return

        self._data_ct = get_tif_image_data(paths[0])
        self._data_ca = get_tif_image_data(paths[1])
        self._data_cb = get_tif_image_data(paths[2])
        self._data_cc = get_tif_image_data(paths[3])
        self._data_sa = get_tif_image_data(paths[4])
        self._data_sb = get_tif_image_data(paths[5])
        self._data_sc = get_tif_image_data(paths[6])

        # to flag the CT values =-9 (corresponding to no data in the tif)
        flagged = np.where(self._data_ct == -9)
        self._data_ct[flagged] = np.nan
        self._data_ca[flagged] = np.nan
        self._data_cb[flagged] = np.nan
        self._data_cc[flagged] = np.nan
        self._data_sa[flagged] = np.nan
        self._data_sb[flagged] = np.nan
        self._data_sc[flagged] = np.nan

        # This step is important for calculation of image coordinates
        # self._data.ice_conc = np.flipud(self._data.ice_conc)
        # self._data_ct = np.flipud(self._data_ct)
        # self._data_ca = np.flipud(self._data_ca)
        # self._data_cb = np.flipud(self._data_cb)
        # self._data_cc = np.flipud(self._data_cc)
        # self._data_sa = np.flipud(self._data_sa)
        # self._data_sb = np.flipud(self._data_sb)
        # self._data_sc = np.flipud(self._data_sc)

        self.add_handler_message("IC: Loaded IC file: %s (and corresponding CABC, SABC" % paths[0])
        self._current_date = self._requested_date

    def _get_local_repository_filename(self, l2):

        time = datetime.datetime(int(self.year), int(self.month), int(self.day))

        path = Path(self.cfg.local_repository)

        other_icevars = ['CA', 'CB', 'CC', 'SA', 'SB', 'SC']
        for delta in [0, -1, 1, -2, 2, -3, 3, 4, -4, 5, -5, 6, -6, 7, -7]:
            fnames = []
            datestring = (time + datetime.timedelta(delta)).strftime('%Y%m%d')
            year = datestring[:4]
            month = datestring[4:6]
            str_filename = path / str(year) / str(month) / "merged_" + datestring + "_CT.tif"
            if str_filename.is_file():
                timedelta = delta
                fnames.append(str_filename)
                for o in other_icevars:
                    fnames.append(path / str(year) / str(month) / 'merged_' + datestring + '_' + o + '.tif')
                return fnames, timedelta
                break
            else:
                fnames = ['', '', '', '', '', '']
                timedelta = np.nan
        return fnames, timedelta

    def XXYYGrids(self):
        # Returns 2 km EASE2 grid XX and YY value

        vec_Y = np.arange(5400000, -5400000, -2000) - 1000
        vec_X = np.arange(-5400000, 5400000, 2000) + 1000
        [XX, YY] = np.meshgrid(vec_X, vec_Y)
        return XX,YY

    def _get_sic_ic_track(self, l2):
        # Convert grid/track coordinates to grid projection coordinates
        kwargs = self.cfg.options[l2.hemisphere].projection
        # p = Proj(**kwargs)
        data_ct = list()
        data_ca = list()
        data_cb = list()
        data_cc = list()
        data_sa = list()
        data_sb = list()
        data_sc = list()

        EASE_proj = pyproj.Proj('+proj=laea +lat_0=90 +lon_0=0 +ellps=WGS84 +datum=WGS84 +units=m')
        latlon_proj = pyproj.Proj('+proj=latlong +ellps=WGS84 +datum=WGS84 +no_defs')

        vec_x, vec_y = pyproj.transform(latlon_proj, EASE_proj, l2.track.longitude, l2.track.latitude)

        # Convert track projection coordinates to image coordinates
        XX, YY = self.XXYYGrids()
        xOrigin = XX[0][0]
        yOrigin = YY[0][0]
        pixelWidth = XX[0][1]-XX[0][0]
        pixelHeight = YY[1][0]-YY[0][0]

        for coords in zip(vec_x,vec_y):
            x = coords[0]
            y = coords[1]
            xOffset = int((x - xOrigin) / pixelWidth)
            yOffset = int((y - yOrigin) / pixelHeight)

            stuf = self._data_ct[yOffset][xOffset]
            dct = str(self._data_ct[yOffset][xOffset])  # AMANDINE: why this conversion to string?
            dca = self._data_ca[yOffset][xOffset]
            dcb = self._data_cb[yOffset][xOffset]
            dcc = self._data_cc[yOffset][xOffset]
            dsa = self._data_sa[yOffset][xOffset]
            dsb = self._data_sb[yOffset][xOffset]
            dsc = self._data_sc[yOffset][xOffset]

            # No value in the tif corresponds to -9
            if (dca == -9) and (dcb == -9) and (dcc == -9):
                #print "ICECHART, all is -9"
                dca = dct
            data_ct.append(dct)
            data_ca.append(dca)
            data_cb.append(dcb)
            data_cc.append(dcc)
            data_sa.append(dsa)
            data_sb.append(dsb)
            data_sc.append(dsc)

        return [data_ct, data_ca, data_cb, data_cc, data_sa, data_sb, data_sc]


class ICA(AuxdataBaseClass):

    def __init__(self, *args, **kwargs):

        super(ICA, self).__init__(*args, **kwargs) #MUOKS20170517
        self._data_ct = None
        self._data_ca = None
        self._data_cb = None
        self._data_cc = None
        self._data_sa = None
        self._data_sb = None
        self._data_sc = None
        self.error.caller_id = self.__class__.__name__

    def get_l2_track_vars(self, l2):
        self._get_requested_date(l2)
        self._get_data(l2)
        if not self.error.status:
            ic_sics = self._get_sic_ic_track(l2)
            self.register_auxvar("ica_ct", "icechart_aari_ct", ic_sics[0], None)
            self.register_auxvar("ica_ca", "icechart_aari_ca", ic_sics[1], None)
            self.register_auxvar("ica_cb", "icechart_aari_cb", ic_sics[2], None)
            self.register_auxvar("ica_cc", "icechart_aari_cc", ic_sics[3], None)
            self.register_auxvar("ica_sa", "icechart_aari_sa", ic_sics[4], None)
            self.register_auxvar("ica_sb", "icechart_aari_sb", ic_sics[5], None)
            self.register_auxvar("ica_sc", "icechart_aari_sc", ic_sics[6], None)
        else:
            self.register_auxvar("ica_ct", "icechart_aari_ct", None, None)
            self.register_auxvar("ica_ca", "icechart_aari_ca", None, None)
            self.register_auxvar("ica_cb", "icechart_aari_cb", None, None)
            self.register_auxvar("ica_cc", "icechart_aari_cc", None, None)
            self.register_auxvar("ica_sa", "icechart_aari_sa", None, None)
            self.register_auxvar("ica_sb", "icechart_aari_sb", None, None)
            self.register_auxvar("ica_sc", "icechart_aari_sc", None, None)

    def _get_requested_date(self, l2):
        """ Use first timestamp as reference, date changes are ignored """
        year = l2.track.timestamp[0].year
        month = l2.track.timestamp[0].month
        day = l2.track.timestamp[0].day
        self._requested_date = [year, month, day]

    def _get_data(self, l2):
        """ Loads file from local repository only if needed """
        if self._requested_date == self._current_date:
            # Data already loaded, nothing to do
            self.add_handler_message("ICA: tif already present")
            return
        paths, timedelta = self._get_local_repository_filename(l2)
        # print 'AARI PATHS 0', paths[0]

        # Validation
        if not Path(paths[0]).is_file():
            self._msg = "ICA: File not found: %s " % paths[0]
            self.error.add_error("auxdata_missing_sic", self._msg)
            return

        self._data_ct = get_tif_image_data(paths[0])
        self._data_ca = get_tif_image_data(paths[1])
        self._data_cb = get_tif_image_data(paths[2])
        self._data_cc = get_tif_image_data(paths[3])
        self._data_sa = get_tif_image_data(paths[4])
        self._data_sb = get_tif_image_data(paths[5])
        self._data_sc = get_tif_image_data(paths[6])
        flagged = np.where(self._data_ct == -5)
        # NOTE August 2017 exception added for missing AARI years
        try:
            self._data_ct[flagged] = np.nan
            self._data_ca[flagged] = np.nan
            self._data_cb[flagged] = np.nan
            self._data_cc[flagged] = np.nan
            self._data_sa[flagged] = np.nan
            self._data_sb[flagged] = np.nan
            self._data_sc[flagged] = np.nan
        except TypeError:
            # print 'No AARI icechart for this day'
            pass
        # This step is important for calculation of image coordinates
        # self._data.ice_conc = np.flipud(self._data.ice_conc)
        # self._data_ct = np.flipud(self._data_ct)
        # self._data_ca = np.flipud(self._data_ca)
        # self._data_cb = np.flipud(self._data_cb)
        # self._data_cc = np.flipud(self._data_cc)
        # self._data_sa = np.flipud(self._data_sa)
        # self._data_sb = np.flipud(self._data_sb)
        # self._data_sc = np.flipud(self._data_sc)
        self.add_handler_message("ICA: Loaded IC file: %s (and corresponding CABC, SABC" % paths[0])
        self._current_date = self._requested_date

    def _get_local_repository_filename(self, l2):

        time = datetime.datetime(int(self.year), int(self.month), int(self.day))

        path = Path(self.cfg.local_repository)

        other_icevars = ['CA', 'CB', 'CC', 'SA', 'SB', 'SC']
        for delta in [0, -1, 1, -2, 2, -3, 3, 4, -4, 5, -5, 6, -6, 7, -7]:
            fnames = []
            datestring = (time + datetime.timedelta(delta)).strftime('%Y%m%d')
            year = datestring[:4]
            month = datestring[4:6]
            str_filename = path / str(year) / str(month) / "merged_aari_" + datestring+"_CT.tif"
            if str_filename.is_file():
                timedelta = delta
                fnames.append(str_filename)
                for o in other_icevars:
                    fnames.append(path / str(year) / str(month) / 'merged_aari_' + datestring + '_' + o +'.tif')
                return fnames, timedelta
                break
            else:
                fnames = ['', '', '', '', '', '']
                timedelta = np.nan
        return fnames, timedelta

    def XXYYGrids(self):
        # Returns 2 km EASE2 grid XX and YY values

        import numpy as np

        vec_Y = np.arange(5400000, -5400000, -2000) - 1000
        vec_X = np.arange(-5400000, 5400000, 2000) + 1000
        [XX, YY] = np.meshgrid(vec_X, vec_Y)
        return XX, YY

    def _get_sic_ic_track(self, l2):
        # Convert grid/track coordinates to grid projection coordinates
        # kwargs = self.cfg.option[l2.hemisphere].projection
        #p = Proj(**kwargs)
        data_ct = list()
        data_ca = list()
        data_cb = list()
        data_cc = list()
        data_sa = list()
        data_sb = list()
        data_sc = list()

        EASE_proj = pyproj.Proj('+proj=laea +lat_0=90 +lon_0=0 +ellps=WGS84 +datum=WGS84 +units=m')
        latlon_proj = pyproj.Proj('+proj=latlong +ellps=WGS84 +datum=WGS84 +no_defs')

        vec_x, vec_y = pyproj.transform(latlon_proj, EASE_proj, l2.track.longitude, l2.track.latitude)

        # Convert track projection coordinates to image coordinates
        XX, YY = self.XXYYGrids()
        xOrigin = XX[0][0]
        yOrigin = YY[0][0]
        pixelWidth = XX[0][1]-XX[0][0]
        pixelHeight = YY[1][0]-YY[0][0]

        for coords in zip(vec_x,vec_y):
            x = coords[0]
            y = coords[1]
            xOffset = int((x - xOrigin) / pixelWidth)
            yOffset = int((y - yOrigin) / pixelHeight)

            try:
                stuf =  self._data_ct[yOffset][xOffset]
                data_ct.append(str(self._data_ct[yOffset][xOffset]))
                data_ca.append((self._data_ca[yOffset][xOffset]))
                data_cb.append((self._data_cb[yOffset][xOffset]))
                data_cc.append((self._data_cc[yOffset][xOffset]))
                data_sa.append((self._data_sa[yOffset][xOffset]))
                data_sb.append((self._data_sb[yOffset][xOffset]))
                data_sc.append((self._data_sc[yOffset][xOffset]))
            except TypeError:
                #print 'No AARI icechart for this day'
                data_ct.append(None)
                data_ca.append(None)
                data_cb.append(None)
                data_cc.append(None)
                data_sa.append(None)
                data_sb.append(None)
                data_sc.append(None)

        return [data_ct,data_ca,data_cb,data_cc,data_sa,data_sb,data_sc]


def get_tif_image_data(path):
    """ Open the tif file and return its content """
    im = Image.open(str(path))
    pixel_list = im.getdata()
    width, height = im.size
    data = np.reshape(pixel_list, (width, height))
    return data


def convert_polygon_icechart(ice_chart_gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Code of this function is adopted from the Auto Ice Challenge
    #  # -- File info -- #
    #  __author__ = 'Andreas R. Stokholm'
    #  __contributors__ = 'Tore Wulf'
    #  __copyright__ = ['Technical University of Denmark', 'European Space Agency']
    #  __contact__ = ['stokholm@space.dtu.dk']
    #  __version__ = '0.0.1'
    #  __date__ = '2022-09-20'

    Original polygon_icechart in ASIP3 scenes consists of codes to a lookup table `polygon_codes`.

    This function looks up codes and converts them. 3 variables in the xr scene are created; SIC, SOD and FLOE.
    For SOD and FLOE the partial sea ice concentration is used to determine whether there is a dominant category
    in a polygon. The SOD and FLOE are created using the lookup tables in utils, which dictate the conversion
    from ice code to class, As multiple codes can be converted into a single class, these partial concentrations
     must also be added. In addition, empty codes, 'not filled values' and unknowns are replaced appropriately.

    Parameters
    ----------
    ice_chart_gdf :
        ice chart in form of a pandas geo dataframe.
    """

    to_numpy_kwargs = dict(dtype=object, copy=False, na_value=-9)
    try:
        codes = ice_chart_gdf[['CT', 'CT', 'CA', 'SA', 'FA', 'CB', 'SB', 'FB', 'CC', 'SC', 'FC']].to_numpy(**to_numpy_kwargs)
        # ToDo: Fix Nan values get converted to very high numbers
        codes = codes.astype(float).astype(int)
        poly_type = ice_chart_gdf['POLY_TYPE'].to_numpy(**to_numpy_kwargs)
    except:
        codes = ice_chart_gdf[['CT', 'CT', 'CA', 'SA', 'FA', 'CB', 'SB', 'FB', 'CC', 'SC', 'FC']].to_numpy(**to_numpy_kwargs)

        codes[codes == '-'] = '-9'
        codes[codes == '5C'] = '05'
        codes[codes == '9-'] = '-9'
        codes[codes == '9C'] = '09'
        codes[codes == '9P'] = '09'

        # ToDo: Fix Nan values get converted to very high numbers
        codes = codes.astype(float).astype(int)
        poly_type = ice_chart_gdf['POLY_TYPE'].to_numpy(dtype=object, copy=False, na_value=-9)

    # Convert codes to classes for Total and Partial SIC. (SIGRID code is replaced by classes defined above
    converted_codes = copy.deepcopy(codes)
    for key, value in SIC_LOOKUP.items():
        if type(key) == int:
            for partial_idx in SIC_LOOKUP['sic_partial_idx']:
                tmp = converted_codes[:, partial_idx]
                if key in tmp:
                    converted_codes[:, partial_idx][np.where((tmp == key))] = value

            tmp = converted_codes[:, SIC_LOOKUP['total_sic_idx']]
            if key in tmp:
                converted_codes[:, SIC_LOOKUP['total_sic_idx']][np.where((tmp == key))[0]] = value

    # Find where partial concentration is empty but total SIC exist.
    # ToDo: Only set partial concentration to total concentration if there is only one partial sod and floe
    ice_ct_ca_empty = np.logical_and(
        converted_codes[:, SIC_LOOKUP['total_sic_idx']] > SIC_LOOKUP[0],
        converted_codes[:, SIC_LOOKUP['sic_partial_idx'][0]] == ICECHART_NOT_FILLED_VALUE)
    # Assign total SIC to partial concentration when empty.
    converted_codes[:, SIC_LOOKUP['sic_partial_idx'][0]][ice_ct_ca_empty] = \
        converted_codes[:, SIC_LOOKUP['total_sic_idx']][ice_ct_ca_empty]

    # Convert codes to classes for partial SOD.
    for key, value in SOD_LOOKUP.items():
        if type(key) == int:
            for partial_idx in SOD_LOOKUP['sod_partial_idx']:
                tmp = converted_codes[:, partial_idx]
                if key in tmp:
                    converted_codes[:, partial_idx][np.where((tmp == key))] = value

    # Convert codes to classes for partial FLOE.
    for key, value in FLOE_LOOKUP.items():
        if type(key) == int:
            for partial_idx in FLOE_LOOKUP['floe_partial_idx']:
                tmp = converted_codes[:, partial_idx]
                if key in tmp:
                    converted_codes[:, partial_idx][np.where((tmp == key))] = value

    sic_T = np.full(codes.shape[0], -99)
    sic_A = np.full(codes.shape[0], -99)
    sod_A = np.full(codes.shape[0], -99)
    floe_A = np.full(codes.shape[0], -99)
    sic_B = np.full(codes.shape[0], -99)
    sod_B = np.full(codes.shape[0], -99)
    floe_B = np.full(codes.shape[0], -99)
    sic_C = np.full(codes.shape[0], -99)
    sod_C = np.full(codes.shape[0], -99)
    floe_C = np.full(codes.shape[0], -99)

    # Find and replace all codes with SIC, SOD and FLOE.
    with np.errstate(divide='ignore', invalid='ignore'):
        for i in range(codes.shape[0]):
            sic_T[i] = converted_codes[i, SIC_LOOKUP['total_sic_idx']]
            sic_A[i] = converted_codes[i, SIC_LOOKUP['sic_partial_idx'][0]]
            sod_A[i] = converted_codes[i, SOD_LOOKUP['sod_partial_idx'][0]]
            floe_A[i] = converted_codes[i, FLOE_LOOKUP['floe_partial_idx'][0]]
            sic_B[i] = converted_codes[i, SIC_LOOKUP['sic_partial_idx'][1]]
            sod_B[i] = converted_codes[i, SOD_LOOKUP['sod_partial_idx'][1]]
            floe_B[i] = converted_codes[i, FLOE_LOOKUP['floe_partial_idx'][1]]
            sic_C[i] = converted_codes[i, SIC_LOOKUP['sic_partial_idx'][2]]
            sod_C[i] = converted_codes[i, SOD_LOOKUP['sod_partial_idx'][2]]
            floe_C[i] = converted_codes[i, FLOE_LOOKUP['floe_partial_idx'][2]]

    # Add masked pixels for ambiguous polygons.
    sod_A[sod_A == SOD_LOOKUP['invalid']] = SOD_LOOKUP['mask']
    floe_A[floe_A == FLOE_LOOKUP['invalid']] = FLOE_LOOKUP['mask']
    sod_B[sod_B == SOD_LOOKUP['invalid']] = SOD_LOOKUP['mask']
    floe_B[floe_B == FLOE_LOOKUP['invalid']] = FLOE_LOOKUP['mask']
    sod_C[sod_C == SOD_LOOKUP['invalid']] = SOD_LOOKUP['mask']
    floe_C[floe_C == FLOE_LOOKUP['invalid']] = FLOE_LOOKUP['mask']

    # Ensure water is identical across charts.
    # ToDo: Same for partial concentrations
    sod_A[sic_T == SIC_LOOKUP[0]] = SOD_LOOKUP['water']
    floe_A[sic_T == SIC_LOOKUP[0]] = FLOE_LOOKUP['water']
    sod_B[sic_T == SIC_LOOKUP[0]] = SOD_LOOKUP['water']
    floe_B[sic_T == SIC_LOOKUP[0]] = FLOE_LOOKUP['water']
    sod_C[sic_T == SIC_LOOKUP[0]] = SOD_LOOKUP['water']
    floe_C[sic_T == SIC_LOOKUP[0]] = FLOE_LOOKUP['water']

    sic_T[np.where(sic_T == ICECHART_UNKNOWN)] = SIC_LOOKUP['mask']
    sic_A[np.where(sic_A == ICECHART_UNKNOWN)] = SIC_LOOKUP['mask']
    sod_A[np.where(sod_A == ICECHART_UNKNOWN)] = SOD_LOOKUP['mask']
    floe_A[np.where(floe_A == ICECHART_UNKNOWN)] = FLOE_LOOKUP['mask']
    sic_B[np.where(sic_B == ICECHART_UNKNOWN)] = SIC_LOOKUP['mask']
    sod_B[np.where(sod_B == ICECHART_UNKNOWN)] = SOD_LOOKUP['mask']
    floe_B[np.where(floe_B == ICECHART_UNKNOWN)] = FLOE_LOOKUP['mask']
    sic_C[np.where(sic_C == ICECHART_UNKNOWN)] = SIC_LOOKUP['mask']
    sod_C[np.where(sod_C == ICECHART_UNKNOWN)] = SOD_LOOKUP['mask']
    floe_C[np.where(floe_C == ICECHART_UNKNOWN)] = FLOE_LOOKUP['mask']

    ice_chart_gdf['SIC_T'] = sic_T
    ice_chart_gdf['SIC_A'] = sic_A
    ice_chart_gdf['SOD_A'] = sod_A
    ice_chart_gdf['Floe_A'] = floe_A

    ice_chart_gdf['SIC_B'] = sic_B
    ice_chart_gdf['SOD_B'] = sod_B
    ice_chart_gdf['Floe_B'] = floe_B

    ice_chart_gdf['SIC_C'] = sic_C
    ice_chart_gdf['SOD_C'] = sod_C
    ice_chart_gdf['Floe_C'] = floe_C

    return ice_chart_gdf
