# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 12:29:13 2016

@author: shendric
"""

import scipy.ndimage as ndimage
from pyproj import Proj


def get_along_track_gridval(gridval, gridlon, gridlat, gridproj_dict,
                            griddim, orbitlon, orbitlat):

    # Establish grid projection and get extent in projection coordinates
    p = Proj(**gridproj_dict)
    x, y = p(gridlon, gridlat)

    # Get orbit position in grid projection coordinates
    l2x, l2y = p(orbitlon, orbitlat)

    # Convert track projection coordinates to image coordinates
    # x: 0 < n_lines; y: 0 < n_cols
    dim = griddim
    x_min = x[dim.n_lines-1, 0]
    y_min = y[dim.n_lines-1, 0]
    ix, iy = (l2x-x_min)/dim.dx, (l2y-y_min)/dim.dy

    # Extract along track data from grid
    orbitval = ndimage.map_coordinates(gridval, [iy, ix], order=0)

    return orbitval
