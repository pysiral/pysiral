# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 20:52:32 2016

@author: Stefan
"""

from pysiral.l1bpreproc import L1bPreProc
from pysiral.l1bdata import L1bConstructor
from pysiral.sentinel3.iotools import Sentinel3FileList

import time


class Sentinel3PreProc(L1bPreProc):

    def __init__(self):

        super(Sentinel3PreProc, self).__init__(self.__class__.__name__)

    def _get_input_files_local_machine_def(self, time_range, version):
        """
        Gets the input data files from default data repository that is
        specified in `local_machine_def.yaml`
        """

        # Get the path defined in local_machine_def.yaml
        l1b_repository = self._pysiral_config.local_machine.l1b_repository
        s3_l1b_repository = l1b_repository[self.mission_id][version].sral

        # Get the list of files (SAR and SIN in chronological order)
        # for the case of this test only for one month of data
        s3_files = Sentinel3FileList()
        s3_files.folder = s3_l1b_repository
        s3_files.search(time_range)

        # Transfer list of l1b input files
        self._l1b_file_list = s3_files.sorted_list

    def _get_l1bdata_ocean_segments(self, filename):
        """ Returns the source Sentinel-3 l1b data as a list of
        pysiral.L1bdata objects. """

        # Read CryoSat-2 Header
        l1b = L1bConstructor(self._pysiral_config)
        l1b.mission = self.mission_id
        l1b.filename = filename

        t0 = time.time()
        l1b.construct()
        t1 = time.time()
        self.log.info("- Parsed source file in %.3g seconds" % (t1 - t0))

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

        # Sentinel-3a WAT files cover a half-orbit from pole to plot
        # (with gaps). This approach is adapted from pre-processing
        # the Envisat half-orbits
        l1b_list = []
        seg1, seg2 = self.extract_polar_segments_from_halforbit(l1b)

        # Threshold for splitting
        seconds_threshold = self._mdef.max_connected_files_timedelta_seconds

        for seg in [seg1, seg2]:

            # Get polar ocean segments
            if seg is not None:

                # Trim non ocean margins
                l1b_roi = self.trim_non_ocean_data(seg)

                # Test for unlikely case of no ocean data
                if l1b_roi is None:
                    return []

                # Split orbits at larger non-ocean segments
                l1b_segments = self.split_at_large_non_ocean_segments(l1b_roi)

                # There are data gaps in the L2WAT files, therefore split
                # at timestamp discontinuities
                l1b_segments = self.split_at_time_discontinuities(
                        l1b_segments, seconds_threshold, trim_non_ocean=True)

                # Some minor gaps seem to remain in the data
                # -> bring to regular grid
                for l1b_segment in l1b_segments:
                    try:
                        l1b_segment.detect_and_fill_gaps()
                    except AttributeError:
                        stop

                # Add to l1b stack
                try:
                    l1b_list.extend(l1b_segments)
                except TypeError:
                    self.log.info("- no ocean data in file")

        return l1b_list

    @property
    def mission_id(self):
        # mission_id can be sentinel3a, ...
        return self._jobdef.mission_id
