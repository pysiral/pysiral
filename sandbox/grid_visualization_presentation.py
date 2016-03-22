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


from plib.map.basemap import cmBlueMarble
from matplotlib.image import pil_to_array


# AWI eisblau #00ace5
# AWI tiefblau #003e6e
# AWI grau 1 #4b4b4d
# AWI grau 2 #bcbdbf


def grid_visualization_presentation():

    create_background_image()
    create_map()


def create_background_image():

    plt.figure("Background Image", figsize=(12, 6))
    plt.gca().set_position([0, 0, 1, 1])
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='i')
    m.drawmapboundary(color="none", fill_color='#003e6e')
    m.fillcontinents(color='#4b4b4d', lake_color='#4b4b4d')
    coastlines = get_landcoastlines(m, color="#bcbdbf", linewidth=0.05)
    plt.gca().add_collection(coastlines)
    plt.savefig("temp.png", dpi=1200)
    plt.close()

    background = Image.open("temp.png")
    foreground_array = get_foreground_array(background.size)
    foreground_array_wet = get_foreground_array(background.size, cutoff=0.8)
#    mask = np.array(background)
#    wet_pixel_list = np.where(mask == [0, 62, 110, 255])
#    foreground_array[wet_pixel_list] = foreground_array_wet[wet_pixel_list]
#    foreground_array[wet_pixel_list] = [100, 62, 110, 255]
    foreground = Image.fromarray(foreground_array)
    background.paste(foreground, (0, 0), foreground)
    background.save('temp2.png', 'PNG')


def get_foreground_array(size, cutoff=0.5, color=[0, 0, 0], alpha_ref=200):
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


def create_map():
    lat_0 = 80
    lon_0 = 0
    h = 5000.
    grid_keyw = {"dashes": (None, None), "color": "#ffffff",
                 "linewidth": 0.05, "latmax": 85, "zorder": 120}
    plt.figure("Grid Data Visualization",  dpi=300, facecolor="#4b4b4d")
    m = Basemap(projection='nsper', lon_0=lon_0, lat_0=lat_0,
                satellite_height=h*1000., resolution='l')
    m.warpimage("temp2.png", scale=1.0, zorder=100)
    m.drawparallels(np.arange(-90., 120., 15.), **grid_keyw)
    m.drawmeridians(np.arange(0., 420., 30.), **grid_keyw)
    plt.show()


def fig2data(fig):
    """
    @brief Convert a Matplotlib figure to a 4D numpy array with RGBA channels
           and return it
    @param fig a matplotlib figure
    @return a numpy 3D array of RGBA values
    """
    # draw the renderer
    fig.canvas.draw()
    # Get the RGBA buffer from the figure
    w, h = fig.canvas.get_width_height()
    buf = np.fromstring(fig.canvas.tostring_argb(), dtype=np.uint8)
    buf.shape = (w, h, 4)
    # canvas.tostring_argb give pixmap in ARGB mode.
    # Roll the ALPHA channel to have it in RGBA mode
    buf = np.roll(buf, 3, axis=2)
    return buf


def fig2img(fig):
    """
    @brief Convert a Matplotlib figure to a PIL Image in RGBA format
           and return it
    @param fig a matplotlib figure
    @return a Python Imaging Library ( PIL ) image
    """
    # put the figure pixmap into a numpy array
    buf = fig2data(fig)
    w, h, d = buf.shape
    return Image.fromstring("RGBA", (w, h), buf.tostring())


def get_landcoastlines(basemap, color="0.0", linewidth=1):
    landpolygons = np.where(np.array(basemap.coastpolygontypes) == 1)[0]
    landsegs = []
    for index in landpolygons:
        landsegs.append(basemap.coastsegs[index])
    landcoastlines = LineCollection(landsegs, antialiaseds=(1,))
    landcoastlines.set_color(color)
    landcoastlines.set_linewidth(linewidth)
    return landcoastlines


def draw_shading(basemap):
    import tempfile
    import os
    from PIL import image
    tempdir = tempfile.gettempdir()

    pass


if __name__ == "__main__":
    grid_visualization_presentation()
