# -*- coding: utf-8 -*-
"""
Created on Sun Aug 20 13:09:50 2017

@author: Stefan
"""


from pysiral.l1bpreproc import L1bPreProc
from pysiral.l1bdata import L1bConstructor
from pysiral.icesat.iotools import ICESatGLAH13Repository

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
        pysiral.L1bdata objects (one for each orbit) """

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