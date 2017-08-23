# -*- coding: utf-8 -*-
"""
Created on Sun Aug 20 13:09:50 2017

@author: Stefan
"""


from pysiral.l1bpreproc import L1bPreProc
from pysiral.l1bdata import L1bConstructor
from pysiral.icesat.iotools import ICESatGLAH13Repository

import numpy as np
import time


class ICESatPreProc(L1bPreProc):

    def __init__(self):

        super(ICESatPreProc, self).__init__(self.__class__.__name__)

    def _get_input_files_local_machine_def(self, time_range, version):
        """
        Gets the input data files from default data repository that is
        specified in `local_machine_def.yaml`
        """

        # Get the path defined in local_machine_def.yaml
        pathdef = self._pysiral_config.local_machine
        local_icesat_repository = pathdef.l1b_repository.icesat[version]

        # Get the list of input files
        glah13_repo = ICESatGLAH13Repository(local_icesat_repository.glah13)
        glah13_files = glah13_repo.get_glah13_hdfs(time_range)

        # Transfer list of l1b input files
        self._l1b_file_list = glah13_files

    def _get_l1bdata_ocean_segments(self, filename):
        """ Returns the source ICESat GLAH13 data as a list of
        pysiral.L1bdata objects (one for each orbit segment) """

        try:
            t0 = time.time()
            l1b = L1bConstructor(self._pysiral_config)
            l1b.mission = "icesat"
            l1b.filename = filename
            l1b.construct()
            t1 = time.time()

        except SystemExit:
            return []

        self.log.info("- Parsed source file in %.3g seconds" % (t1 - t0))

        # The current l1b object contains Arctic and Antarctic data
        # from multiple orbits. We therefore only need to split the
        # segments that can be assumed to be all over marine areas

        # index of pre-defined orbit segments, which should be written
        # in one output file
        orbit_segment_id = l1b.classifier.get_parameter("orbit_segment_id")

#        import matplotlib.pyplot as plt
#        plt.figure()
#        plt.plot(orbit_segment_id)
#        plt.show()
#        stop

        # track id from glah13 input data that can be used to trace the
        # orbit segment back to the original orbit number
        track_id = l1b.classifier.get_parameter("track_id")
        track_ids = np.unique(track_id)

        # Get number of segments
        segment_ids = np.unique(orbit_segment_id)
        n_segments = len(segment_ids)
        self.log.info("- %g orbit segments in file" % n_segments)

        # Get the orbits (there should be a director orbit - track relation)
        orbit_list = l1b.info.orbit
        orbits = [int(o) for o in orbit_list.split(",")]
        track_orbit_dict = {}
        for i, track_id_number in enumerate(track_ids):
            track_orbit_dict[track_id_number] = orbits[i]

        l1b_list = []
        for i, segment_id in enumerate(segment_ids):

            msg = "- Extracting segment (%g of %g)" % (i+1, n_segments)
            self.log.info(msg)

            # Get the list of indices and create a new l1b object
            index_list = np.where(orbit_segment_id == segment_id)[0]
            l1b_subset = l1b.extract_subset(index_list)

            # Get the orbit number for this particalur subset
            segment_track_ids = np.unique(track_id[index_list])

            # Use only the first track id to the get orbit number
            segment_orbit = track_orbit_dict[segment_track_ids[0]]
            l1b_subset.info.set_attribute("orbit", segment_orbit)

            # Trim subset
            l1b_subset_trimmed = self.trim_non_ocean_data(l1b_subset)

#            import matplotlib.pyplot as plt
#
#            plt.figure()
#            plt.scatter(l1b_subset.time_orbit.timestamp,
#                     l1b_subset.time_orbit.longitude)
#            plt.show()
#            stop

            if l1b_subset_trimmed is None:
                continue

            # Store subset in output list
            l1b_list.append(l1b_subset_trimmed)

        return l1b_list
