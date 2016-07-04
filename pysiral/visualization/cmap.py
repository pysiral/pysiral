# -*- coding: utf-8 -*-
"""
Created on Mon Jul 04 15:20:29 2016

@author: shendric
"""

from matplotlib.colors import ListedColormap


def surface_type_cmap():

    from pysiral.surface_type import SurfaceType
    surface_type = SurfaceType()

    # Color definitions for different surface types
    rgb_list = [
        "#bcbdbf",  # Unkown
        "#003e6e",  # Ocean
        "#00ace5",  # Lead
        "#DDA0DD",  # Polynya
        "#76FF7A",  # Sea Ice
        "#4b4b4d",  # Lakes
        "#FFCC00",  # Land Ice
        "#826644",  # Land
        "#FF1700"   # invalid
        ]
    cmap = ListedColormap(rgb_list, name='surface_type')
    labels = [surface_type.name(i) for i in range(9)]
    return cmap, labels


def radar_mode_cmap():

    from pysiral.config import RadarModes
    radar_modes = RadarModes()

    # Color definitions for different surface types
    rgb_list = [
        (0.050383, 0.029803, 0.527975, 1.0),  # lrm
        (0.798216, 0.280197, 0.469538, 1.0),  # sar
        (0.940015, 0.975158, 0.131326, 1.0)  # sin
        ]
    cmap = ListedColormap(rgb_list, name='radar_modes')
    labels = [radar_modes.name(i) for i in range(3)]
    return cmap, labels


def is_valid_cmap():

    # Color definitions for different surface types
    rgb_list = ["#FF1700", "#76FF7A"]
    cmap = ListedColormap(rgb_list, name='is_valid')
    labels = ["invalid", "valid"]
    return cmap, labels
