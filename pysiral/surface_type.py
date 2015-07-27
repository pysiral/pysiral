# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 11:25:04 2015

@author: Stefan
"""

class SurfaceType(object):
    """
    Container for surface type information.

    Possible classifications (Adapted from CryoSat-2 conventions)
        - open ocean
        - closed sea/lakes

        - lead
        - large lead/polynya
        - continental ice
        - land
    """
    _SURFACE_TYPE_DICT = {
        "open_ocean": 0,
        "closed_sea": 1,
        "continental_ice": 2,
        "sea_ice": 3,
        "lead": 4,
        "polynya": 5}

    def __init__(self):
        self._surface_type_flag = None
        self._received_attributes = []


class IceType(object):
    """
    Container for ice type information

    Possible classifications
        - young thin ice
        - first year ice
        - multi year ice
    """
    _ICE_TYPE_DICT = {
        "thin_ice": 0,
        "fyi": 1,
        "myi": 2}

    def __init__(self):
        self._ice_type_flag = None
