# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 20:52:32 2016

@author: Stefan
"""

from pysiral.l1bpreproc import L1bPreProc
from pysiral.cryosat2.iotools import CryoSat2MonthlyFileListAllModes
from pysiral.l1bdata import L1bConstructor

import time


class CryoSat2PreProc(L1bPreProc):

    def __init__(self):

        super(CryoSat2PreProc, self).__init__(self.__class__.__name__)

    def _get_input_files_local_machine_def(self, time_range, version):
        """
        Gets the input data files from default data repository that is
        specified in `local_machine_def.yaml`
        """

        # Get the path defined in local_machine_def.yaml
        pathdef = self._pysiral_config.local_machine
        cryosat_l1b_repository = pathdef.l1b_repository.cryosat2[version]

        # Get the list of files (SAR and SIN in chronological order)
        # for the case of this test only for one month of data
        cryosat2_files = CryoSat2MonthlyFileListAllModes()
        cryosat2_files.folder_sar = cryosat_l1b_repository.sar
        cryosat2_files.folder_sin = cryosat_l1b_repository.sin
        cryosat2_files.search(time_range)

        # Transfer list of l1b input files
        self._l1b_file_list = cryosat2_files.sorted_list

    def _get_l1bdata_ocean_segments(self, filename):
        """
        Returns the source CryoSat-2 l1b data as a list of
        pysiral.L1bdata objects.
        """

        # Read CryoSat-2 Header
        l1b = L1bConstructor(self._pysiral_config)
        l1b.mission = "cryosat2"
        l1b.filename = filename
        l1b.get_header_info()

        # 1) Check if file has ocean data in the polar regions at all
        has_polar_ocean = self.region_is_arctic_or_antarctic_ocean(l1b)

        # 2) Check if desires region is matched
        if self._jobdef.hemisphere == "global":
            matches_region = True
        elif l1b.info.region_name == self._jobdef.hemisphere:
            matches_region = True
        else:
            matches_region = False

        if not has_polar_ocean or not matches_region:
            return None

        # Only now read the full data set
        t0 = time.time()
        l1b.construct()
        t1 = time.time()
        self.log.info("- Parsed source file in %.3g seconds" % (t1 - t0))

        # Reduce the waveform count of SIN waveforms to the number of SAR
        # waveforms (for data mergung)
        sar_bin_count = {"baseline-b": 128, "baseline-c": 256}
        bin_count = sar_bin_count[l1b.info.mission_data_version.lower()]
        if l1b.radar_modes == "sin":
            l1b.reduce_waveform_bin_count(bin_count)

        # This step is CryoSat-2 specific, since orbit segments do not cover
        # both hemisphere
        l1b_roi = self.trim_single_hemisphere_segment_to_polar_region(l1b)

        # Trim non ocean margins
        l1b_roi = self.trim_non_ocean_data(l1b_roi)

        # Split orbits at larger non-ocean segments
        l1b_list = self.split_at_large_non_ocean_segments(l1b_roi)

        return l1b_list
