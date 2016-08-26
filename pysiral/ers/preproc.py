# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 20:52:32 2016

@author: Stefan
"""

from pysiral.l1bpreproc import L1bPreProc
from pysiral.ers.iotools import ERSFileList
from pysiral.l1bdata import L1bConstructor


class ERSPreProc(L1bPreProc):

    def __init__(self):

        super(ERSPreProc, self).__init__(self.__class__.__name__)

    def _get_input_files_local_machine_def(self, time_range, version):
        """
        Gets the input data files from default data repository that is
        specified in `local_machine_def.yaml`
        """

        # mission_id can be ers1 or ers2
        mission_id = self._jobdef.mission_id

        # Get the l1b data repository
        pathdef = self._pysiral_config.local_machine.l1b_repository
        ersx_l1b_repository = pathdef[mission_id][version].sgdr

        # Get the list of SGDR files for specific time range
        ers_files = ERSFileList()
        ers_files.folder = ersx_l1b_repository
        ers_files.search(time_range)
        self._l1b_file_list = ers_files.sorted_list

    def _get_l1bdata_ocean_segments(self, filename):
        """
        Returns the source ERS l1b data as a list of
        pysiral.L1bdata objects.
        """

        # Read the envisat SGDR file
        l1b = L1bConstructor(self._pysiral_config)
        l1b.mission = self._jobdef.mission_id
        l1b.filename = filename
        l1b.construct()

        # Extract relevant segments over ocean
        l1b_list = []
        seg1, seg2 = self.extract_polar_segments_from_halforbit(l1b)

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

                # Add to l1b stack
                try:
                    l1b_list.extend(l1b_segments)
                except TypeError:
                    self.log.info("- no ocean data in file")

        return l1b_list