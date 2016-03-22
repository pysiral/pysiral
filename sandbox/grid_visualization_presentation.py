# -*- coding: utf-8 -*-
"""
Sandbox script to create appealing visualization of grid data
(e.g. http://ozonewatch.gsfc.nasa.gov/ozone_maps/images/)

Created on Mon Mar 21 20:20:08 2016

@author: Stefan
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.basemap import Basemap

import os
from PIL import Image

from plib.files.scientific_ff import ReadNC
from plib.geo.grid import GeoPcolorGrid

from plib.map.basemap import cmBlueMarble
from matplotlib.image import pil_to_array


# AWI eisblau #00ace5
# AWI tiefblau #003e6e
# AWI grau 1 #4b4b4d
# AWI grau 2 #bcbdbf


def grid_visualization_presentation():

    # create_background_image()
    grid = get_grid()
    create_map(grid)


def get_grid():
    projection = {
        "proj": "laea",
        "lat_0": -90.0,
        "lon_0": 0.0,
        "ellps": "WGS84",
        "datum": "WGS84",
        "units": "m"}
    nc_file = os.path.join(r"D:\awi\product\optimal_interp\v1",
                           r"cs2_smos_ice_thickness_20150302_20150308.nc")
    grid = ReadNC(nc_file)
    grid.pcolor = GeoPcolorGrid(grid.longitude, grid.latitude)
    grid.pcolor.calc_from_proj(**projection)
    return grid


def create_background_image():

    plt.figure("Background Image", figsize=(12, 6))
    plt.gca().set_position([0, 0, 1, 1])
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='i')
    m.drawmapboundary(color="none", fill_color='#003e6e')
    m.fillcontinents(color='#4b4b4d', lake_color='#4b4b4d')
    coastlines = get_landcoastlines(m, color="#bcbdbf", linewidth=0.05)
    plt.gca().add_collection(coastlines)
    plt.savefig("temp.png", dpi=300)
    plt.close()

    background = Image.open("temp.png")
    background_array = np.array(background)
    foreground_array = get_diffuse_highlighted_image(
        -45, -60, background.size, n=1)
    foreground_array_wet = get_diffuse_highlighted_image(
        -45, -60, background.size, n=5)
    # XXX: This is working but really really slow
    for i in range(background.size[0]):
        for j in range(background.size[1]):
            if np.array_equal(background_array[j, i, :], [0,  62, 110, 255]):
                foreground_array[j, i, :] = foreground_array_wet[j, i, :]

    foreground = Image.fromarray(foreground_array)
    background.paste(foreground, (0, 0), foreground)
    background.save('temp2.png', 'PNG')


def get_foreground_array(size, cutoff=0.5, color=[0, 0, 0], alpha_ref=150):
    foreground_template = Image.new("RGBA", size, (0, 0, 0, 200))
    foreground_array = np.array(foreground_template)
    latimin, latimax = 0, int(cutoff*size[1])
    scale_factor = 256./np.float(latimax)
    for i in np.arange(latimin, latimax).astype(int):
        alpha = int(i*scale_factor)
        if alpha > 200:
            alpha = 200
        if alpha < 0:
            alpha = 0
        foreground_array[i, :, :] = [0, 0, 0, alpha]
    return foreground_array


def create_map(grid):
    lat_0 = 80
    lon_0 = 0
    h = 6000.
    grid_keyw = {"dashes": (None, None), "color": "#bcbdbf",
                 "linewidth": 0.1, "latmax": 85, "zorder": 120}
    plt.figure("Grid Data Visualization",  figsize=(8, 8),
               dpi=150, facecolor="#4b4b4b")
    m = Basemap(projection='nsper', lon_0=lon_0, lat_0=lat_0,
                satellite_height=h*1000., resolution='l')
#    m = Basemap(projection='ortho', lon_0=lon_0, lat_0=lat_0,
#                resolution='l')
    m.warpimage("temp2.png", scale=1.0, zorder=100)
    x, y = m(grid.pcolor.longitude, grid.pcolor.latitude)
    data = getattr(grid, "background_thickness")
    dataMask = np.isnan(data)
    data = np.ma.masked_where(dataMask, data)
    m.pcolor(x, y, data, cmap=plt.get_cmap("plasma"), vmin=0, vmax=5,
             zorder=110)

    m.drawparallels(np.arange(-90., 120., 15.), **grid_keyw)
    m.drawmeridians(np.arange(0., 420., 30.), **grid_keyw)
    plt.show()


def get_diffuse_highlighted_image(lon_0, lat_0, size,
                                  color=[0, 0, 0], alpha_ref=255, **kwargs):
    intensity = np.ndarray(shape=(size[0], size[1]))
    # angle = np.ndarray(shape=(size[0], size[1]))
    r0 = [1, np.deg2rad(lat_0), np.deg2rad(lon_0)]
    ref_rgba = (color[0], color[1], color[2], alpha_ref)
    diffuse_template = Image.new("RGBA", size, ref_rgba)
    diffuse_array = np.array(diffuse_template)
    lons = np.deg2rad(np.linspace(-180, 180, size[0]))
    lats = np.deg2rad(np.linspace(-90, 90, size[1]))
    for i in range(size[0]):
        for j in range(size[1]):
            r1 = [1, lats[j], lons[i]]
            intensity[i, j] = get_diffuse_intensity(r0, r1,**kwargs)
            diffuse_array[j, i, 3] = 255-int(intensity[i, j]*255)
#    plt.figure("intensity")
#    plt.imshow(intensity.T)
#    plt.colorbar()
#    plt.figure("angle")
#    plt.imshow(np.rad2deg(angle.T))
#    plt.colorbar()
#    plt.show(block=True)

    return diffuse_array


def get_diffuse_intensity(r0, r1, k=1.0, b=1, n=4):
    v0 = llarToWorld(r0[1], r0[2], 1, 0)
    v1 = llarToWorld(r1[1], r1[2], 1, 0)
    angle = angle_between(v0, v1)
    if angle > np.pi/2.0:
        angle = np.pi/2.0
    intensity = k*np.cos(b*angle)**n
    if intensity < 0:
        intensity = 0
    return intensity# , angle


def get_landcoastlines(basemap, color="0.0", linewidth=1):
    landpolygons = np.where(np.array(basemap.coastpolygontypes) == 1)[0]
    landsegs = []
    for index in landpolygons:
        landsegs.append(basemap.coastsegs[index])
    landcoastlines = LineCollection(landsegs, antialiaseds=(1,))
    landcoastlines.set_color(color)
    landcoastlines.set_linewidth(linewidth)
    return landcoastlines


def spherical_to_cartesian(spherical_vect):
    """
    Convert the spherical coordinate vector [r, theta, phi] to the Cartesian
    vector [x, y, z].

    The parameter r is the radial distance, theta is the polar angle, and phi
    is the azimuth.

    @param spherical_vect:  The spherical coordinate vector [r, theta, phi].
    @type spherical_vect:   3D array or list
    @param cart_vect:       The Cartesian vector [x, y, z].
    @type cart_vect:        3D array or list
    """
    # Trig alias.
    cart_vect = np.copy(spherical_vect)
    sin_theta = np.sin(spherical_vect[1])
    # The vector.
    cart_vect[0] = spherical_vect[0] * np.cos(spherical_vect[2]) * sin_theta
    cart_vect[1] = spherical_vect[0] * np.sin(spherical_vect[2]) * sin_theta
    cart_vect[2] = spherical_vect[0] * np.cos(spherical_vect[1])
    return cart_vect


def llarToWorld(lat, lon, alt, rad):
    # see: http://www.mathworks.de/help/toolbox/aeroblks/llatoecefposition.html
    f = 0                              # flattening
    ls = np.arctan((1 - f)**2 * np.tan(lat))    # lambda
    x = rad * np.cos(ls) * np.cos(lon) + alt * np.cos(lat) * np.cos(lon)
    y = rad * np.cos(ls) * np.sin(lon) + alt * np.cos(lat) * np.sin(lon)
    z = rad * np.sin(ls) + alt * np.sin(lat)
    return [x, y, z]


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = np.arccos(np.dot(v1_u, v2_u))
    if np.isnan(angle):
        if (v1_u == v2_u).all():
            return 0.0
        else:
            return np.pi
    return angle


if __name__ == "__main__":
    grid_visualization_presentation()
