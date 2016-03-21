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

from plib.map.basemap import cmBlueMarble
from matplotlib.image import pil_to_array


# AWI eisblau #00ace5
# AWI tiefblau #003e6e
# AWI grau 1 #4b4b4d
# AWI grau 2 #bcbdbf


def grid_visualization_presentation():

    from PIL import Image

#    folder = r"D:\awi\workspace\python-tools\plib\map\bluemarble\8100x4050"
#    filename = r"world.topo.200401.3x8100x4050.jpg"
#    image_test_file = os.path.join(folder, filename)
#    im = Image.open(image_test_file)
#    im = im.convert('RGBA')
#    imarr = np.array(im)
#    print imarr.shape
#    print imarr.dtype
#    stop

    # template_image = Image.new("RGBA", (1535, 768), (0, 0, 0, 255))
    template_image = Image.new("RGBA", (1024, 512), (0, 0, 0, 255))
    shading_image_array = np.array(template_image)
    for i in range(256):
        shading_image_array[i, :, :] = [255-i, 255-i, 255-i, 255]
    shading_image = Image.fromarray(shading_image_array)
#    plt.figure()
#    plt.imshow(shading_image)
#    plt.show(block=True)
    png_info = shading_image.info
    shading_image.save('temp.png', 'PNG', **png_info)


    grid_keyw = {"dashes": (None, None),
                 "color": "#bcbdbf",
                 "linewidth": 0.5,
                 "latmax": 85,
                 "zorder": 120}

    lat_0 = 80
    lon_0 = 0

    plt.figure("Grid Data Visualization", facecolor="white")
    h = 6500.
#    m = Basemap(projection='nsper', lon_0=lon_0, lat_0=lat_0,
#                satellite_height=h*1000., resolution='l')
    m = Basemap(projection='ortho', lon_0=lon_0, lat_0=lat_0,
                resolution='l')
    m.drawmapboundary(fill_color='#003e6e')
    m.fillcontinents(color='#4b4b4d', lake_color='#4b4b4d')
    coastlines = get_landcoastlines(m, color="#bcbdbf", linewidth=1.0)
    plt.gca().add_collection(coastlines)
    # draw parallels and meridians.
    im = m.warpimage("temp.png", scale=1.0, zorder=100, alpha=0.55)
    # cmBlueMarble(m, basesize="10800x5400", scale=1.0)
    m.drawparallels(np.arange(-90., 120., 15.), **grid_keyw)
    m.drawmeridians(np.arange(0., 420., 30.), **grid_keyw)
    plt.show()


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
