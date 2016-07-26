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


class PlotL2iDataOrbit(DataPlot):

    def __init__(self):

        super(PlotL2iDataOrbit, self).__init__()

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


class PlotL2iAuxiliaryData(DataPlot):

    def __init__(self):

        super(PlotL2iAuxiliaryData, self).__init__()

        # Specific location of ticks (e.g. to locate regions in map)
        self.xtick_locations = None

        # Number of plots
        self.n_plots = 3

        # Figure settings
        self.fig_args = {"facecolor": "white", "figsize": (10, 14)}
        self.axis_bg_color = '0.98'
        self.line_settings = {"lw": 2, "color": "#00ace5"}

    def create_plot(self, l2i, xtick_locations=None):

        # Save the pointer
        self.l2i = l2i

        # Read to create the plot
        self._create_plot()

    def _create_plot(self):

        # Number of sub panels
        n = self.n_plots

        # Create the figure
        gs = gridspec.GridSpec(n, 1)
        self.fig = plt.figure(n, **self.fig_args)
        plt.subplots_adjust(bottom=0.075, top=0.95)

        # Plot sea ice concentration
        ax = plt.subplot(gs[0])
        ax.plot(self.l2i.sea_ice_concentration, **self.line_settings)
        self._set_plot_style(ax, "Sea Ice Concentration", [-5, 105])

        # Plot myi fraction
        ax = plt.subplot(gs[1])
        ax.plot(self.l2i.sea_ice_type, **self.line_settings)
        self._set_plot_style(ax, "Sea Ice Type (MYI Fraction)", [-0.05, 1.05])

        # Plot myi fraction
        ax = plt.subplot(gs[2])
        ax.plot(self.l2i.snow_depth, **self.line_settings)
        self._set_plot_style(ax, "Snow Depth", [0, 0.8])

    def _set_plot_style(self, ax, name, ylim):

        # Labels and axis style
        ax.set_axis_bgcolor(self.axis_bg_color)
        ax.set_title(name)
        ax.yaxis.grid(True, which='major', dashes=(None, None), lw=0.1)
        ax.yaxis.set_tick_params(direction='out')
        ax.set_ylim(ylim)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks([])

        spines_to_remove = ["top", "bottom"]
        for spine in spines_to_remove:
            ax.spines[spine].set_visible(False)


class PlotL2iSurfaceTypeClassification(DataPlot):

    def __init__(self):

        super(PlotL2iSurfaceTypeClassification, self).__init__()

        # Number of records that shall be used for alongtrack histogram
        self.section_hist_size = 100

        # Figure settings
        self.fig_args = {"facecolor": "white", "figsize": (8, 10)}

    def create_plot(self, l2i):

        from pysiral.visualization.cmap import surface_type_cmap

        # Get the colorbar for surface types
        cmap, labels = surface_type_cmap()

        # Get surface type fractions
        fractions = self._get_surface_type_fractions(l2i)

        # Create the figure (with sub grids)
        self.fig = plt.figure(**self.fig_args)
        gridshape = (7, 15)

        # Create the pie chart
        plt.subplot2grid(gridshape, (0, 0), colspan=14, rowspan=4)

        non_zero = np.nonzero(fractions)[0]
        colors = np.array(cmap.colors)

        plt.gca().set_aspect(1.0)
        wedges, texts, autotexts = plt.pie(
            fractions[non_zero], labels=np.array(labels)[non_zero],
            colors=colors[non_zero], autopct='%1.1f%%', startangle=90,
            labeldistance=1.1, pctdistance=0.75)

        for wedge in wedges:
            wedge.set_width(0.5)
            wedge.set_ec("none")


        # Create the along track flag indicator plot
        plt.subplot2grid(gridshape, (4, 0), colspan=15)
        ax = plt.gca()

        aspect = l2i.surface_type.shape[0]/20.
        im = plt.imshow([l2i.surface_type], interpolation=None, aspect=aspect,
                        vmin=-0.5, vmax=8.5, cmap=cmap)

        ax.xaxis.set_tick_params(direction='out')
        ax.xaxis.set_ticks_position('bottom')
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        # ax.spines["bottom"].set_position(("data", 0.75))

        spines_to_remove = ["top", "right", "left", "bottom"]

        for spine in spines_to_remove:
            ax.spines[spine].set_visible(False)

        # Create the label colorbar
        plt.subplot2grid(gridshape, (0, 14), rowspan=4)
        cbar = plt.colorbar(im, orientation="vertical", cax=plt.gca())
        cbar_ticks = np.arange(9)
        cbar.set_ticks(cbar_ticks)
        cbar.set_ticklabels(labels)
        cbar.solids.set_edgecolor("1.0")
        cbar.outline.set_linewidth(5)
        cbar.outline.set_alpha(0.0)
        cbar.ax.tick_params('both', length=0.1, which='major', pad=5)

        # Create the along-track histogram
        sect_fractions = self._get_section_fractions(l2i)
        baseline = np.zeros(shape=(sect_fractions.shape[0]))
        index = np.arange(sect_fractions.shape[0])

        plt.subplot2grid(gridshape, (5, 0), colspan=15)
        ax = plt.gca()

        for i in np.arange(9):
            fractions = np.ravel(sect_fractions[:, i])
            ax.bar(index, fractions, 1., color=colors[i], bottom=baseline,
                   edgecolor="none")
            baseline += fractions

        ax.set_xlim([0, sect_fractions.shape[0]])
        ax.set_ylim([0, 1])
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        spines_to_remove = ["top", "right", "left", "bottom"]

        for spine in spines_to_remove:
            ax.spines[spine].set_visible(False)

        plt.tight_layout()

    def _get_surface_type_fractions(self, l2i):

        # Calculate Fractions
        n = float(l2i.n_records)
        fractions = []
        for i in np.arange(9):
            num = len(np.where(l2i.surface_type == i)[0])
            fractions.append(float(num)/n)

        return np.array(fractions)

    def _get_section_fractions(self, l2i):

        n_sections = int(np.ceil(l2i.n_records/self.section_hist_size))
        n_surface_types = 9

        result = np.ndarray(shape=(n_sections, n_surface_types))

        for section in np.arange(n_sections):
            i0 = section * self.section_hist_size
            i1 = i0 + self.section_hist_size
            surface_type_flag = l2i.surface_type[i0:i1]
            for i in np.arange(n_surface_types):
                num = len(np.where(surface_type_flag == i)[0])
                result[section, i] = float(num)/float(self.section_hist_size)

        return result


