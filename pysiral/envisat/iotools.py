# -*- coding: utf-8 -*-

import os

from collections import deque


class EnvisatFileList(object):
    """
    Class for the construction of a list of Envisat N1 files
    sorted by acquisition time
    XXX: Currently only support order by month/date and not cycle
    """

    def __init__(self):
        self.folder = None
        self.log = None
        self.year = None
        self.month = None
        self.pattern = ".N1"
        self._list = deque([])
        self._sorted_list = []

    def search(self):
        self._get_file_listing()

    @property
    def sorted_list(self):
        return self._sorted_list

    def _get_file_listing(self):
        search_toplevel_folder = self._get_toplevel_search_folder()
        # walk through files
        for dirpath, dirnames, filenames in os.walk(search_toplevel_folder):
            self.log.info("Searching folder: %s" % dirpath)
            # Get the list of all dbl files
            files = [os.path.join(dirpath, fn) for fn in filenames
                     if self.pattern in fn]
            self.log.info("Found %g level-1b SGDR files" % len(files))
            # reform the list that each list entry is of type
            # [full_path, identifier (start_date)] for later sorting
            # of SAR and SIN files
            sublist = sorted(files)
            self._sorted_list.extend(sublist)

    def _get_toplevel_search_folder(self):
        folder = self.folder
        if self.year is not None:
            folder = os.path.join(folder, "%4g" % self.year)
        if self.month is not None and self.year is not None:
            folder = os.path.join(folder, "%02g" % self.month)
        return folder
