import numpy as np
from pyproj import Proj

from matplotlib.collections import LineCollection


class StereoMapAutoProjExtent(object):

    def __init__(self, lons, lats, aspect=1, scale=1.1, res='h'):
        self.aspect = aspect
        self.scale = scale
        self.res = res
        self.lons = lons
        self.lats = lats

    def _compute_proj_parameters_and_extent(self, **kwargs):

        # 1. use approximate center to find actual center
        approx_lon_0 = np.mean([np.nanmin(self.lons), np.nanmax(self.lons)])
        approx_lat_0 = np.mean([np.nanmin(self.lats), np.nanmax(self.lats)])
        p = Proj(proj='stere', lon_0=approx_lon_0, lat_0=approx_lat_0,
                 lat_ts=approx_lat_0, ellps='WGS84')

        # 2. find a circle that just encompasses all points
        xx, yy = p(self.lons, self.lats)
        points = [[xx[i], yy[i]] for i in np.arange(len(xx))]
        c = make_circle(points)
        self.lon_0, self.lat_0 = p(c[0], c[1], inverse=True)
        size = self.scale * 2.0 * c[2]

        # Apply scale and aspect
        self.width = size
        self.height = size
        if self.aspect > 1:
            self.width *= self.aspect
        if self.aspect < 1:
            self.height *= self.aspect

    def get_basemap_args(self, **kwargs):

        self._compute_proj_parameters_and_extent(**kwargs)
        basemap_kwargs = {
                'projection': 'stere',
                'width': np.ceil(self.width),
                'height': np.ceil(self.height),
                'lon_0': self.lon_0,
                'lat_0': self.lat_0,
                'lat_ts': self.lat_0,
                'resolution': self.res}
        return basemap_kwargs


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



#
# Smallest enclosing circle - Library (Python)
#
# Copyright (c) 2017 Project Nayuki
# https://www.nayuki.io/page/smallest-enclosing-circle
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program (see COPYING.txt and COPYING.LESSER.txt).
# If not, see <http://www.gnu.org/licenses/>.
#

import math, random


# Data conventions: A point is a pair of floats (x, y). A circle is a triple of floats (center x, center y, radius).

#
# Returns the smallest circle that encloses all the given points. Runs in expected O(n) time, randomized.
# Input: A sequence of pairs of floats or ints, e.g. [(0,5), (3.1,-2.7)].
# Output: A triple of floats representing a circle.
# Note: If 0 points are given, None is returned. If 1 point is given, a circle of radius 0 is returned.
#
# Initially: No boundary points known
def make_circle(points):
	# Convert to float and randomize order
	shuffled = [(float(p[0]), float(p[1])) for p in points]
	random.shuffle(shuffled)

	# Progressively add points to circle or recompute circle
	c = None
	for (i, p) in enumerate(shuffled):
		if c is None or not is_in_circle(c, p):
			c = _make_circle_one_point(shuffled[0 : i + 1], p)
	return c


# One boundary point known
def _make_circle_one_point(points, p):
	c = (p[0], p[1], 0.0)
	for (i, q) in enumerate(points):
		if not is_in_circle(c, q):
			if c[2] == 0.0:
				c = make_diameter(p, q)
			else:
				c = _make_circle_two_points(points[0 : i + 1], p, q)
	return c


# Two boundary points known
def _make_circle_two_points(points, p, q):
	circ = make_diameter(p, q)
	left = None
	right = None

	# For each point not in the two-point circle
	for r in points:
		if is_in_circle(circ, r):
			continue

		# Form a circumcircle and classify it on left or right side
		cross = _cross_product(p[0], p[1], q[0], q[1], r[0], r[1])
		c = make_circumcircle(p, q, r)
		if c is None:
			continue
		elif cross > 0.0 and (left is None or _cross_product(p[0], p[1], q[0], q[1], c[0], c[1]) > _cross_product(p[0], p[1], q[0], q[1], left[0], left[1])):
			left = c
		elif cross < 0.0 and (right is None or _cross_product(p[0], p[1], q[0], q[1], c[0], c[1]) < _cross_product(p[0], p[1], q[0], q[1], right[0], right[1])):
			right = c

	# Select which circle to return
	if left is None and right is None:
		return circ
	elif left is None:
		return right
	elif right is None:
		return left
	else:
		return left if (left[2] <= right[2]) else right


def make_circumcircle(p0, p1, p2):
	# Mathematical algorithm from Wikipedia: Circumscribed circle
	ax = p0[0]; ay = p0[1]
	bx = p1[0]; by = p1[1]
	cx = p2[0]; cy = p2[1]
	ox = (min(ax, bx, cx) + max(ax, bx, cx)) / 2.0
	oy = (min(ay, by, cy) + max(ay, by, cy)) / 2.0
	ax -= ox; ay -= oy
	bx -= ox; by -= oy
	cx -= ox; cy -= oy
	d = (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by)) * 2.0
	if d == 0.0:
		return None
	x = ox + ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by)) / d
	y = oy + ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax)) / d
	ra = math.hypot(x - p0[0], y - p0[1])
	rb = math.hypot(x - p1[0], y - p1[1])
	rc = math.hypot(x - p2[0], y - p2[1])
	return (x, y, max(ra, rb, rc))


def make_diameter(p0, p1):
	cx = (p0[0] + p1[0]) / 2.0
	cy = (p0[1] + p1[1]) / 2.0
	r0 = math.hypot(cx - p0[0], cy - p0[1])
	r1 = math.hypot(cx - p1[0], cy - p1[1])
	return (cx, cy, max(r0, r1))


_MULTIPLICATIVE_EPSILON = 1 + 1e-14

def is_in_circle(c, p):
	return c is not None and math.hypot(p[0] - c[0], p[1] - c[1]) <= c[2] * _MULTIPLICATIVE_EPSILON


# Returns twice the signed area of the triangle defined by (x0, y0), (x1, y1), (x2, y2).
def _cross_product(x0, y0, x1, y1, x2, y2):
	return (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0)