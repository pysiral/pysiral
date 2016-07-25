# -*- coding: utf-8 -*-
"""
Created on Tue Jul 05 18:56:49 2016

@author: shendric
"""

from pysiral.visualization.dataplot import DataPlot
from pysiral.iotools import get_temp_png_filename

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap
from matplotlib.ticker import MultipleLocator

import numpy as np
import os


class PlotL2IdataOrbit(DataPlot):

    def __init__(self):

        super(PlotL2IdataOrbit, self).__init__()

        # Specific location of ticks (e.g. to locate regions in map)
        self.xtick_locations = None

        # Figure settings
        self.fig_args = {"facecolor": "#ffffff", "figsize": (12, 12)}
        self.grid_keyw = {"dashes": (None, None),
                          "color": "#bcbdbf",
                          "linewidth": 0.5,
                          "latmax": 88, "zorder": 10}

    def create_plot(self, l2i):

        self._create_shading(l2i)
        self._create_map(l2i)

    def _create_shading(self, l2i):

        # Check the hemisphere
        hemisphere = "north"
        if l2i.latitude[0] < 0:
            hemisphere = "south"

        xsize, ysize = 4000, 2000

        lats = np.linspace(-0.5*np.pi, 0.5*np.pi, num=ysize)

        if hemisphere == "south":
            shading = np.zeros(shape=(ysize, xsize))
            indices = np.where(lats > np.deg2rad(-90))[0]
            colormap_name = "Greys"
        elif hemisphere == "north":
            shading = np.ones(shape=(ysize, xsize)) * 2.0
            indices = np.where(lats < np.deg2rad(90))[0]
            colormap_name = "Greys_r"
        for i in indices:
            shading[i, :] = (np.sin(lats[i])+1.0)**3.

        self.temp_filename = get_temp_png_filename()
        fig = plt.figure(figsize=(10, 5))
        ax = plt.axes([0, 0, 1, 1])
        plt.axis('off')
        ax.imshow(np.flipud(shading), cmap=plt.get_cmap(colormap_name),
                  vmin=0, vmax=1.1*np.amax(shading))
        plt.savefig(self.temp_filename, dpi=600)
        plt.close(fig)

    def _create_map(self, l2i):

        lat_0 = np.median(l2i.latitude)
        lon_0 = np.median(l2i.longitude)

        self.fig = plt.figure(**self.fig_args)
        m = Basemap(projection='ortho', lon_0=lon_0, lat_0=lat_0,
                    resolution='i')

        # Fill the continents
        m.fillcontinents(color='#4b4b4d', lake_color='#4b4b4d', zorder=20)

        # draw parallels and meridians.
        m.drawparallels(np.arange(-90., 120., 10.), **self.grid_keyw)
        m.drawmeridians(np.arange(0., 420., 15.), **self.grid_keyw)

        # plot the track
        x, y = m(l2i.longitude, l2i.latitude)

        # plot the track start marker
        m.scatter(x, y, s=10, color="#FF1700", edgecolors="none", zorder=202)
        m.scatter(x[0], y[0], marker="D", s=80, edgecolors="#FF1700",
                  color="#FF1700", zorder=201)

        plt.title("")

        m.warpimage(self.temp_filename, zorder=200, alpha=0.80)
        os.remove(self.temp_filename)
        plt.tight_layout()
