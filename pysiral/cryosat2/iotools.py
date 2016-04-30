# -*- coding: utf-8 -*-

import os

import numpy as np
from collections import deque
from pysiral.logging import DefaultLoggingClass


class CryoSat2FileListAllModes(DefaultLoggingClass):
    """
    Class for the construction of a list of CryoSat-2 SAR/SIN files
    sorted by acquisition time
    """

    def __init__(self):
        super(CryoSat2FileListAllModes, self).__init__(
            "cryosat2-file-list-all-modes")
        self.folder_sar = None
        self.folder_sin = None
        self.year = None
        self.month = None
        self.pattern = ".DBL"
        self._list = deque([])
        self._sorted_list = None

    def search(self):
        self._search_specific_mode_files("sar")
        self._search_specific_mode_files("sin")
        self._sort_mixed_mode_file_list()

    @property
    def sorted_list(self):
        return [item[0] for item in self._sorted_list]

    def _search_specific_mode_files(self, mode):
        search_toplevel_folder = self._get_toplevel_search_folder(mode)
        # walk through files
        for dirpath, dirnames, filenames in os.walk(search_toplevel_folder):
            self.log.info("Searching folder: %s" % dirpath)
            # Get the list of all dbl files
            cs2files = [fn for fn in filenames if self.pattern in fn]
            self.log.info("Found %g %s level-1b files" % (
                len(cs2files), mode))
            # reform the list that each list entry is of type
            # [full_path, identifier (start_date)] for later sorting
            # of SAR and SIN files
            sublist = [self._get_list_item(fn, dirpath) for fn in cs2files]
            self._list.extend(sublist)

    def _get_toplevel_search_folder(self, mode):
        folder = getattr(self, "folder_"+mode)
        if self.year is not None:
            folder = os.path.join(folder, "%4g" % self.year)
        if self.month is not None:
            folder = os.path.join(folder, "%02g" % self.month)
        return folder

    def _get_list_item(self, filename, dirpath):
        return (os.path.join(dirpath, filename), filename.split("_")[6])

    def _sort_mixed_mode_file_list(self):
        dtypes = [('path', object), ('start_time', object)]
        self._sorted_list = np.array(self._list, dtype=dtypes)
        self._sorted_list.sort(order='start_time')
