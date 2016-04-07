import numpy as np
from pyproj import Proj

from matplotlib.collections import LineCollection


class GeoPcolorGrid():

    def __init__(self, lon, lat):

        self.n, self.m = np.shape(lon)
        self.glon = lon
        self.glat = lat
        self.longitude = np.ndarray((self.n+1, self.m+1), dtype=np.float32)
        self.latitude = np.ndarray((self.n+1, self.m+1), dtype=np.float32)

    def calc_from_proj(self, **kwargs):

        p = Proj(**kwargs)
        x, y = p(self.glon, self.glat)

        xp = np.ndarray((self.n+1, self.m+1), dtype=np.float32)
        yp = np.ndarray((self.n+1, self.m+1), dtype=np.float32)

        dx = -0.5 * (x[0, 1] - x[0, 0])
        dy = -0.5 * (y[1, 0] - y[0, 0])

        # Main part
        xp[0:self.n, 0:self.m] = x[0:self.n, 0:self.m] + dx
        yp[0:self.n, 0:self.m] = y[0:self.n, 0:self.m] + dy

        # lower boundary
        xp[self.n, 0:self.m] = x[self.n-1, 0:self.m] + dx
        yp[self.n, 0:self.m] = y[self.n-1, 0:self.m] - dy

        # Right boundary
        xp[0:self.n, self.m] = x[0:self.n, self.m-1] - dx
        yp[0:self.n, self.m] = y[0:self.n, self.m-1] + dy

        # Last tiny piece
        xp[self.n, self.m] = x[self.n-1, self.m-1] - dx
        yp[self.n, self.m] = y[self.n-1, self.m-1] - dy

        self.longitude, self.latitude = p(xp, yp, inverse=True)


def get_landcoastlines(basemap, color="0.0", linewidth=1):
    landpolygons = np.where(np.array(basemap.coastpolygontypes) == 1)[0]
    landsegs = []
    for index in landpolygons:
        landsegs.append(basemap.coastsegs[index])
    landcoastlines = LineCollection(landsegs, antialiaseds=(1,))
    landcoastlines.set_color(color)
    landcoastlines.set_linewidth(linewidth)
    return landcoastlines
