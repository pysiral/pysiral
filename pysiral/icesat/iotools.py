# -*- coding: utf-8 -*-

from pysiral.logging import DefaultLoggingClass
from pysiral.errorhandler import ErrorStatus

from glob import glob
from pathlib import Path


class ICESatGLAH13Repository(DefaultLoggingClass):

    _GLAH13_SEARCH = r"GLAH13_*.H5"

    def __init__(self, local_repository_path):

        # Init class and error handler
        class_name = self.__class__.__name__
        super(ICESatGLAH13Repository, self).__init__(class_name)
        self.error = ErrorStatus(caller_id=class_name)

        # Sanity check on path to local repository
        if Path(local_repository_path).is_dir():
            self._local_repository_path = local_repository_path
        else:
            msg = "Invalid GLAH13 directory: %s" % str(local_repository_path)
            self.error.add_error("invalid-dir", msg)
            self.error.raise_on_error()

    def get_glah13_hdfs(self, time_range):
        search_folder = self._get_full_path(time_range)
        return sorted(Path(search_folder).glob(self._GLAH13_SEARCH))

    def _get_full_path(self, time_range):
        """ Assuming the time range monthly """
        folder = self.local_repository_path
        return Path(folder) / "%04g" % time_range.start.year / "%02g" % time_range.start.month

    @property
    def local_repository_path(self):
        return self._local_repository_path
