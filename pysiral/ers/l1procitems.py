# -*- coding: utf-8 -*-

"""

"""

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"

import numpy as np
import numpy.typing as npt
from typing import List
from loguru import logger

from pysiral.l1bdata import Level1bData
from pysiral.l1preproc.procitems import L1PProcItem


class SGDRRemoveOverlaps(L1PProcItem):

    def __init__(self, **cfg):
        super(SGDRRemoveOverlaps, self).__init__(**cfg)
        for option_name in self.required_options:
            if option_name not in self.cfg.keys():
                logger.error(f"Missing option: {option_name} -> Leading Edge Quality will not be computed")

    def apply(self, l1: Level1bData):
        """
        Mandatory class of a L1 preproceessor item. Has not any effect for this class

        :param l1:
        :return:
        """
        pass

    def apply_list(self, l1_list: List[Level1bData]) -> List[Level1bData]:
        """
        Ensure that there is no overlap for the SGDR segments.

        :param l1_list:
        :return:
        """

        breakpoint()

    @property
    def required_options(self):
        return []