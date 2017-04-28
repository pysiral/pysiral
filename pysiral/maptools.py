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

    def calc_from_proj(self, xres=None, yres=None, **kwargs):

        p = Proj(**kwargs)
        x, y = p(self.glon, self.glat)

        xp = np.ndarray((self.n+1, self.m+1), dtype=np.float32)
        yp = np.ndarray((self.n+1, self.m+1), dtype=np.float32)

        if xres is None:
            dx = -0.5 * (x[0, 1] - x[0, 0])
        else:
            dx = -0.5*xres

        if yres is None:
            dy = -0.5 * (y[1, 0] - y[0, 0])
        else:
            dy = -0.5*yres

        # Main part
        xp[0:self.n, 0:self.m] = x[0:self.n, 0:self.m] + dx
        yp[0:self.n, 0:self.m] = y[0:self.n, 0:self.m] - dy

        # lower boundary
        xp[self.n, 0:self.m] = x[self.n-1, 0:self.m] + dx
        yp[self.n, 0:self.m] = y[self.n-1, 0:self.m] + dy

        # Right boundary
        xp[0:self.n, self.m] = x[0:self.n, self.m-1] - dx
        yp[0:self.n, self.m] = y[0:self.n, self.m-1] - dy

        # Last tiny piece
        xp[self.n, self.m] = x[self.n-1, self.m-1] - dx
        yp[self.n, self.m] = y[self.n-1, self.m-1] + dy

#        import matplotlib.pyplot as plt
#
#        ii, jj = 310, 310
#        plt.figure()
#        plt.scatter([x[ii, jj]], [y[ii, jj]], s=120, color="red")
#        plt.scatter([xp[ii, jj], xp[ii+1, jj], xp[ii, jj+1], xp[ii+1, jj+1]],
#                    [yp[ii, jj], yp[ii+1, jj], yp[ii, jj+1], yp[ii+1, jj+1]],
#                    color="black")
#        plt.show()

#        plt.figure()
#        plt.scatter(x, y, color="red")
#        plt.scatter(xp, yp, color="black")
#        plt.xlim(np.nanmin(x), np.nanmax(x))
#        plt.ylim(np.nanmin(y), np.nanmax(y))
#        plt.show()
#
#        stop

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
