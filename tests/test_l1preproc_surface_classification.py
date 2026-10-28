# -*- coding: utf-8 -*-
"""
Tests for the lazy loading of the CryoTEMPO surface classification mask
(`pysiral.l1preproc.procitems.masks.L1PCryoTEMPOSurfaceClassification`).

The tests use small synthetic mask files with the same layout as the
250 m mask files (month x yc x xc uint8 grids on a polar LAEA grid with
+/- 5400 km extent), so that no auxiliary data is required.
"""

import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
import xarray as xr
from loguru import logger
from pydantic import computed_field

logger.disable("pysiral")

from pysiral.core.legacy_classes import AttrDict
from pysiral.grid import GridTrajectoryExtract
from pysiral.l1preproc import L1PreProcCustomOrbitSegment
from pysiral.l1preproc.procitems import masks

N_GRID, DX_GRID = 216, 50000

# Unpatched xarray readers (the readers used by the processor item are replaced by spies in the tests)
XR_OPEN_DATASET, XR_LOAD_DATASET = xr.open_dataset, xr.load_dataset


class SmallGridDefinition(masks.GridDefinition):
    """ Same extent as the 250 m mask grid, but with 50 km grid cells """

    @computed_field
    def dimension(self) -> dict:
        return {"n_cols": N_GRID, "n_lines": N_GRID, "dx": DX_GRID, "dy": DX_GRID}


def write_mask_file(filepath: Path, offset: int) -> None:
    """ Write a synthetic surface classification mask with a known (month dependent) pattern """
    month, yc, xc = np.meshgrid(np.arange(12), np.arange(N_GRID), np.arange(N_GRID), indexing="ij")
    surface_classification = ((month + yc + xc + offset) % 5).astype(np.uint8)
    sea_ice_region = ((yc * 3 + xc + offset) % 19).astype(np.uint8)
    land_ocean_classification = ((yc + xc) % 2).astype(np.uint8)
    ds = xr.Dataset(
        {
            "surface_classification": (("month", "yc", "xc"), surface_classification),
            "sea_ice_region": (("yc", "xc"), sea_ice_region[0]),
            "land_ocean_classification": (("yc", "xc"), land_ocean_classification[0]),
        },
        coords={"month": np.arange(1, 13)},
    )
    ds.to_netcdf(filepath)


def track(lat_start: float, lat_stop: float, n: int = 2000) -> tuple:
    """ A trajectory with varying longitude between two latitudes """
    latitude = np.linspace(lat_start, lat_stop, n)
    longitude = np.linspace(-170.0, 170.0, n)
    return longitude, latitude


class TestL1PCryoTEMPOSurfaceClassification(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.lookup_directory = Path(self.tmpdir.name)
        write_mask_file(self.lookup_directory / "cryotempo-surface_classification-nh-250m-v1p0.nc", offset=0)
        write_mask_file(self.lookup_directory / "cryotempo-surface_classification-sh-250m-v1p0.nc", offset=2)

        # Record the files read by the processor item
        self.opened_files = []

        def get_spy(reader):
            def spy(filepath, *args, **kwargs):
                self.opened_files.append(Path(filepath).name)
                return reader(filepath, *args, **kwargs)
            return spy

        patchers = [
            mock.patch.object(
                masks.SurfaceClassificationMaskConfig,
                "lookup_directory",
                property(lambda _: self.lookup_directory)
            ),
            mock.patch.object(masks.xr, "open_dataset", get_spy(XR_OPEN_DATASET)),
            mock.patch.object(masks.xr, "load_dataset", get_spy(XR_LOAD_DATASET)),
        ]
        for patcher in patchers:
            patcher.start()
            self.addCleanup(patcher.stop)
        self.addCleanup(self.tmpdir.cleanup)

    def get_procitem(self) -> masks.L1PCryoTEMPOSurfaceClassification:
        return masks.L1PCryoTEMPOSurfaceClassification(
            grid_def={"nh": SmallGridDefinition(center_latitude=90), "sh": SmallGridDefinition(center_latitude=-90)}
        )

    def get_reference(self, procitem, hemisphere: str, longitude, latitude, month_number: int) -> tuple:
        """ Previous implementation: lookup on the fully loaded mask of a single hemisphere """
        filepath = self.lookup_directory / f"cryotempo-surface_classification-{hemisphere}-250m-v1p0.nc"
        ds = XR_LOAD_DATASET(filepath)
        grid2track = GridTrajectoryExtract(longitude, latitude, procitem.cfg.grid_def[hemisphere].to_attrdict())
        surface_classification = grid2track.get_from_grid_variable(
            ds.surface_classification.values[month_number - 1, :, :], outside_value=-1
        )
        sea_ice_region = grid2track.get_from_grid_variable(ds.sea_ice_region.values, outside_value=-1)
        return surface_classification, sea_ice_region

    def testNoMaskDataReadAtInit(self):
        self.get_procitem()
        self.assertEqual(self.opened_files, [])

    def testNorthernTrajectoryMatchesFullMaskLookup(self):
        procitem = self.get_procitem()
        procitem.set_target_hemisphere(["north"])
        longitude, latitude = track(45.0, 88.0)
        for month_number in (10, 11):
            surface_classification, sea_ice_region = procitem.get_trajectory(longitude, latitude, month_number)
            ref_surface_classification, ref_sea_ice_region = self.get_reference(
                procitem, "nh", longitude, latitude, month_number
            )
            self.assertEqual(surface_classification.dtype, ref_surface_classification.dtype)
            np.testing.assert_array_equal(surface_classification, ref_surface_classification)
            np.testing.assert_array_equal(sea_ice_region, ref_sea_ice_region)

    def testOnlyCurrentMonthIsKept(self):
        procitem = self.get_procitem()
        longitude, latitude = track(45.0, 88.0)
        procitem.get_trajectory(longitude, latitude, 10)
        procitem.get_trajectory(longitude, latitude, 11)
        data = procitem.data["nh"]
        self.assertEqual(data._month_number, 11)
        self.assertEqual(data._surface_classification.shape, (N_GRID, N_GRID))
        self.assertEqual(data._sea_ice_region.shape, (N_GRID, N_GRID))

    def testNonTargetHemisphereIsNotRead(self):
        procitem = self.get_procitem()
        procitem.set_target_hemisphere(["north"])
        longitude, latitude = track(-80.0, 80.0)
        surface_classification, sea_ice_region = procitem.get_trajectory(longitude, latitude, 10)
        self.assertFalse(any("-sh-" in filename for filename in self.opened_files))
        self.assertIsNone(procitem.data["sh"]._surface_classification)

        # Southern records are treated as outside the mask
        is_south = latitude < 0.0
        outside_surface_classification, outside_sea_ice_region = self.get_reference(
            procitem, "nh", np.array([0.0]), np.array([-80.0]), 10
        )
        np.testing.assert_array_equal(surface_classification[is_south], outside_surface_classification[0])
        np.testing.assert_array_equal(sea_ice_region[is_south], outside_sea_ice_region[0])

        # Northern records are unaffected by the southern records in the same trajectory
        ref_surface_classification, ref_sea_ice_region = self.get_reference(procitem, "nh", longitude, latitude, 10)
        np.testing.assert_array_equal(surface_classification[~is_south], ref_surface_classification[~is_south])
        np.testing.assert_array_equal(sea_ice_region[~is_south], ref_sea_ice_region[~is_south])

    def testGlobalTargetUsesMaskOfEachRecordsHemisphere(self):
        procitem = self.get_procitem()
        procitem.set_target_hemisphere(["north", "south"])
        longitude, latitude = track(-80.0, 80.0)
        surface_classification, sea_ice_region = procitem.get_trajectory(longitude, latitude, 10)
        is_south = latitude < 0.0
        for hemisphere, is_hemisphere in (("nh", ~is_south), ("sh", is_south)):
            ref_surface_classification, ref_sea_ice_region = self.get_reference(
                procitem, hemisphere, longitude, latitude, 10
            )
            np.testing.assert_array_equal(
                surface_classification[is_hemisphere], ref_surface_classification[is_hemisphere]
            )
            np.testing.assert_array_equal(sea_ice_region[is_hemisphere], ref_sea_ice_region[is_hemisphere])

    def testPreProcessorPassesTargetHemisphere(self):
        cfg = AttrDict(**{
            "polar_ocean": {"target_hemisphere": ["north"]},
            "processing_items": [{
                "label": "CryoTEMPO Surface Classification",
                "stage": "post_source_file",
                "module_name": "pysiral.l1preproc.procitems.masks",
                "class_name": "L1PCryoTEMPOSurfaceClassification",
                "options": {"version": "v1p0"},
            }]
        })
        preproc = L1PreProcCustomOrbitSegment(None, None, cfg)
        procitem, _ = preproc.processor_item_dict["post_source_file"][0]
        self.assertEqual(procitem.target_hemisphere, ["nh"])
        self.assertEqual(self.opened_files, [])


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestL1PCryoTEMPOSurfaceClassification)
    unittest.TextTestRunner(verbosity=2).run(suite)
