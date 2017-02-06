# -*- coding: utf-8 -*-
"""
Created on Tue Jul 05 18:56:49 2016

@author: shendric
"""

from pysiral.iotools import get_temp_png_filename
from pysiral.visualization.dataplot import DataPlot

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap
from matplotlib.ticker import MultipleLocator

import numpy as np
import os


class PlotL1bdataOrbit(DataPlot):

    def __init__(self):

        super(PlotL1bdataOrbit, self).__init__()

        # Specific location of ticks (e.g. to locate regions in map)
        self.xtick_locations = None

        # Figure settings
        self.fig_args = {"facecolor": "#ffffff", "figsize": (12, 12)}
        self.grid_keyw = {"dashes": (None, None),
                          "color": "#bcbdbf",
                          "linewidth": 0.5,
                          "latmax": 88, "zorder": 10}

    def create_plot(self, l1b):

        self._create_shading(l1b)
        self._create_map(l1b)

    def _create_shading(self, l1b):

        # Check the hemisphere
        hemisphere = "north"
        if l1b.time_orbit.latitude[0] < 0:
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

    def _create_map(self, l1b):

        lat_0 = np.median(l1b.time_orbit.latitude)
        lon_0 = np.median(l1b.time_orbit.longitude)

        self.fig = plt.figure(**self.fig_args)
        m = Basemap(projection='ortho', lon_0=lon_0, lat_0=lat_0,
                    resolution='i')

        # Fill the continents
        m.fillcontinents(color='#4b4b4d', lake_color='#4b4b4d', zorder=20)

        # draw parallels and meridians.
        m.drawparallels(np.arange(-90., 120., 10.), **self.grid_keyw)
        m.drawmeridians(np.arange(0., 420., 15.), **self.grid_keyw)

        # plot the track
        x, y = m(l1b.time_orbit.longitude, l1b.time_orbit.latitude)

        # plot the track start marker
        m.scatter(x, y, s=10, color="#FF1700", edgecolors="none", zorder=202)
        m.scatter(x[0], y[0], marker="D", s=80, edgecolors="#FF1700",
                  color="#FF1700", zorder=201)

        plt.title("")

        m.warpimage(self.temp_filename, zorder=200, alpha=0.80)
        os.remove(self.temp_filename)
        plt.tight_layout()


class PlotL1bdataWaveform(DataPlot):

    def __init__(self):

        super(PlotL1bdataWaveform, self).__init__()

        # Specific location of ticks (e.g. to locate regions in map)
        self.xtick_locations = None

        self.elevation_limit = [-75, 50]

        # Figure settings
        self.fig_args = {"facecolor": "white", "figsize": (12, 6)}
        self.axis_bg_color = '0.98'
#        self.line_settings = {"lw": 2, "color": "#00ace5"}

    def create_plot(self, l1b):

        # Set up the figure and axis location
        self.fig = plt.figure(**self.fig_args)
        ax = plt.gca()
        ax.set_axis_bgcolor(self.axis_bg_color)
        ax.set_position([0.1, 0.05, 0.7, 0.9])

        # Align the waveform based on the window delay
        image, elevation_range = align_echo_power(
            l1b.waveform.power, l1b.waveform.range, l1b.time_orbit.altitude,
            elevation_limit=self.elevation_limit)

        # Compute the extent if the image
        image_extent = (0, len(l1b.time_orbit.altitude),
                        np.amin(elevation_range), np.amax(elevation_range))

        # Plot the image as background image
        im = ax.imshow(
            np.log10(image).transpose(), cmap=plt.get_cmap("magma"),
            interpolation='none', origin='lower', extent=image_extent,
            aspect='auto')

        # Axes style
        ax.yaxis.set_tick_params(direction='out')
        ax.yaxis.set_ticks_position('left')
        ax.set_ylabel("Elevation (m)")
        ax.xaxis.set_ticks([])
        ax.spines["left"].set_position(("data", -20))

        spines_to_remove = ["top", "right", "bottom", "left"]
        for spine in spines_to_remove:
            ax.spines[spine].set_color("#4b4b4b")
            ax.spines[spine].set_linewidth(0.25)

        # Add a colobar
        cax = plt.gcf().add_axes([0.85, 0.05, 0.02, 0.9])
        cb = plt.colorbar(im, cax=cax)
        cb.set_label("Echo Power log10(Watt)")


class PlotL1bdataWaveformFlags(DataPlot):

    def __init__(self):

        super(PlotL1bdataWaveformFlags, self).__init__()

        # Specific location of ticks (e.g. to locate regions in map)
        self.xtick_locations = None

        self.elevation_limit = [-75, 50]

        # Figure settings
        self.fig_args = {"facecolor": "white", "figsize": (12, 6)}

    def create_plot(self, l1b):

        from pysiral.visualization.cmap import (
            is_valid_cmap, surface_type_cmap, radar_mode_cmap)

        # Create the plot
        self.fig = plt.figure(**self.fig_args)

        # set up the axis
        # position of flag images
        axgrid = gridspec.GridSpec(8, 1)

        # position of axis is tied to waveform plot
        axgrid.update(left=0.1, right=0.8)

        # position of pie fraction charts
        y0 = 0.05
        pie_height = 0.4
        pie_width = 0.2
        pie_position = [
            (0.05, y0, pie_width, pie_height),
            (0.35, y0, pie_width, pie_height),
            (0.65, y0, pie_width, pie_height)]

        # position of label colobars
        cb_height = 0.4
        cb_width = 0.02
        cb_position = [
            (0.27, y0, cb_width, cb_height),
            (0.57, y0, cb_width, cb_height),
            (0.87, y0, cb_width, cb_height)]

        # prepare the information for looping over the flags
        flags = [l1b.waveform.is_valid,
                 l1b.waveform.radar_mode,
                 l1b.surface_type.flag]

        cmap_info = [is_valid_cmap(), radar_mode_cmap(), surface_type_cmap()]

        title = ["waveform is_valid flag", "waveform radar mode",
                 "l1b surface type flag"]

        # loop over the three flags
        for i in np.arange(3):

            # Plot the flag as time line
            ax = plt.subplot(axgrid[i])

            # Plot the flag images
            flag = flags[i]
            cmap, labels = cmap_info[i]
            im = plt.imshow([flag], interpolation=None, aspect='auto',
                            vmin=-0.5, vmax=len(labels)-0.5, cmap=cmap)

            # Axis style
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            spines_to_remove = ["top", "right", "bottom", "left"]
            for spine in spines_to_remove:
                ax.spines[spine].set_visible(False)

            # Plot the colorbar with labels
            cb_x0, cb_y0, cb_width, pie_height = cb_position[i]
            cax = plt.gcf().add_axes([cb_x0, cb_y0, cb_width, cb_height])

            cbar = plt.colorbar(im, orientation="vertical", cax=cax)
            cbar_ticks = np.arange(len(labels))
            cbar.set_ticks(cbar_ticks)
            cbar.set_ticklabels(labels)
            cbar.solids.set_edgecolor("1.0")
            cbar.outline.set_alpha(0.0)
            cbar.ax.tick_params('both', length=0.1, which='major', pad=5)

            # Plot a pie diagram with fractions
            pie_x0, pie_y0, pie_width, pie_height = pie_position[i]
            pax = plt.gcf().add_axes([pie_x0, pie_y0, pie_width, pie_height])

            pax.set_title(title[i])
            n = float(len(flag))
            fractions = []
            for i in np.arange(9):
                num = len(np.where(flag == i)[0])
                fractions.append(float(num)/n)

            pie_fractions = np.array(fractions)
            pie_colors = np.array(cmap.colors)

            non_zero = np.nonzero(pie_fractions)[0]
            pie_fractions = pie_fractions[non_zero]
            pie_colors = pie_colors[non_zero]

            wedges, texts, autotexts = pax.pie(
                pie_fractions, colors=pie_colors,
                autopct='%1.1f%%', startangle=90, labeldistance=1.1,
                pctdistance=0.75)

            for wedge in wedges:
                wedge.set_width(0.5)
                wedge.set_ec("none")


class PlotL1bdataRangeCorrections(DataPlot):

    def __init__(self):

        super(PlotL1bdataRangeCorrections, self).__init__()

        # Specific location of ticks (e.g. to locate regions in map)
        self.xtick_locations = None

        # Indices of l1b geophysical range corrections to plot
        self.grc_indices = None

        # Number of panels
        self.max_number_plots = None

        # Horizontal line at fixed interval to illustrate grc magnitude
        self.range_indicator = 0.05

        # Figure settings
        self.fig_args = {"facecolor": "white", "figsize": (10, 14)}
        self.axis_bg_color = '0.98'
        self.line_settings = {"lw": 2, "color": "#00ace5"}

    def create_plot(self, l1b, grc_indices=None, max_number_plots=None,
                    xtick_locations=None):

        # Save the pointer
        self.l1b = l1b

        # Check which range correction shall be plotted
        # (this option is required for l1bdata_report)
        if grc_indices is None:
            self.grc_indices = range(0, len(l1b.correction.parameter_list))
        else:
            self.grc_indices = grc_indices

        # This option is related to index_range, but determines the space
        # reserved for each plot (max_number_plots >= len(grc_indices))
        if max_number_plots is None:
            self.max_number_plots = len(l1b.correction.parameter_list)
        else:
            self.max_number_plots = max_number_plots
        if self.max_number_plots < len(self.grc_indices):
            self.max_number_plot = len(self.grc_indices)

        # Read to create the plot
        self._create_plot()

    def _create_plot(self):

        # Number of sub panels
        n = self.max_number_plots

        # Create the figure
        gs = gridspec.GridSpec(n, 1)
        self.fig = plt.figure(n, **self.fig_args)
        plt.subplots_adjust(bottom=0.075, top=0.95)

        # Loop over parameters
        for i, grc_index in enumerate(self.grc_indices):

            ax = plt.subplot(gs[i])

            # Retrieve and plot the range correction
            grc, name = self.l1b.correction.get_parameter_by_index(grc_index)
            ax.plot(grc, **self.line_settings)

            # Labels and axis style
            ax.set_axis_bgcolor(self.axis_bg_color)
            ax.set_title(name)
            # ax.yaxis.set_minor_locator(MultipleLocator(self.range_indicator))
            ax.yaxis.grid(True, which='minor')
            ax.yaxis.set_tick_params(direction='out')
            ax.yaxis.set_ticks_position('left')
            ax.xaxis.set_ticks([])

            spines_to_remove = ["top", "right", "bottom"]
            for spine in spines_to_remove:
                ax.spines[spine].set_visible(False)


class PlotL1bdataWaveformClassifier(DataPlot):

    def __init__(self):

        super(PlotL1bdataWaveformClassifier, self).__init__()

        # Specific location of ticks (e.g. to locate regions in map)
        self.xtick_locations = None

        # Indices of l1b geophysical range corrections to plot
        self.clf_indices = None

        # Number of panels
        self.max_number_plots = None

        # Figure settings
        self.fig_args = {"facecolor": "white", "figsize": (10, 14)}
        self.axis_bg_color = '0.98'
        self.marker_args = {"s": 1, "edgecolors": "none", "color": "#00ace5"}
        self.barh_args = {"align": 'center', "edgecolor": "none",
                          "color": "#00ace5"}

    def create_plot(self, l1b, clf_indices=None, max_number_plots=None,
                    xtick_locations=None):

        # Save the pointer
        self.l1b = l1b

        # Check which range correction shall be plotted
        # (this option is required for l1bdata_report)
        if clf_indices is None:
            self.grc_indices = range(0, len(l1b.classifier.parameter_list))
        else:
            self.clf_indices = clf_indices

        # This option is related to index_range, but determines the space
        # reserved for each plot (max_number_plots >= len(grc_indices))
        if max_number_plots is None:
            self.max_number_plots = len(l1b.classifier.parameter_list)
        else:
            self.max_number_plots = max_number_plots
        if self.max_number_plots < len(self.clf_indices):
            self.max_number_plot = len(self.clf_indices)

        # Read to create the plot
        self._create_plot()

    def _create_plot(self):

        import matplotlib.gridspec as gridspec

        # Number of sub panels
        n = self.max_number_plots
        classifier_names = sorted(self.l1b.classifier.parameter_list)

        # Create the figure
        gs = gridspec.GridSpec(n, 5)
        self.fig = plt.figure(**self.fig_args)
        plt.subplots_adjust(bottom=0.075, top=0.95)

        x = np.arange(self.l1b.n_records)

        # Loop over parameters
        for i, clf_index in enumerate(self.clf_indices):

            # Scatterplot axes
            ax0 = plt.subplot(gs[i, 0:-1])

            # Histogram axes
            ax1 = plt.subplot(gs[i, -1], sharey=ax0)

            parameter_name = classifier_names[clf_index]

            # Retrieve and plot the classifier
            classifier = self.l1b.classifier.get_parameter(parameter_name)
            ax0.scatter(x, classifier, **self.marker_args)

            # Axis Style
            ax0.set_axis_bgcolor(self.axis_bg_color)
            ax0.set_title(parameter_name)
            ax0.yaxis.set_tick_params(direction='out')
            ax0.yaxis.set_ticks_position('left')
            ax0.xaxis.set_ticks([])
            ax0.set_xlim(0, len(classifier))

            spines_to_remove = ["top", "right", "bottom"]
            for spine in spines_to_remove:
                ax0.spines[spine].set_visible(False)

            # Compute histogram
            valid = np.where(np.isfinite(classifier))[0]
            hist, bin_edges = np.histogram(
                classifier[valid], bins=500, density=True)
            bin_width = bin_edges[1]-bin_edges[0]
            bin_center = bin_edges[0:-1] + 0.5*bin_width
            ax1.barh(bin_center, hist, height=bin_width, **self.barh_args)

            ax1.set_axis_bgcolor(self.axis_bg_color)
            ax1.yaxis.set_tick_params(direction='out')
            ax1.yaxis.set_ticks_position('right')
            ax1.xaxis.set_ticks([])
            spines_to_remove = ["top", "left", "bottom"]
            for spine in spines_to_remove:
                ax1.spines[spine].set_visible(False)


def align_echo_power(power, range, altitude, elevation_limit=None):

    from scipy.interpolate import interp1d

    n_range_bins = len(power[0, :])
    n_records = len(power[:, 0])

    elevation = np.repeat(altitude, n_range_bins)
    elevation = elevation.reshape(n_records, n_range_bins)
    elevation -= range

    range_step = range[0, 1] - range[0, 0]

    elevation_range = np.arange(np.amin(elevation)-range_step*0.5,
                                np.amax(elevation)+range_step*0.5,
                                range_step)

    aligned_power = np.ndarray(shape=(n_records, len(elevation_range)))

    for i in np.arange(n_records):
        f = interp1d(elevation[i, :].flatten(), power[i, :].flatten(),
                     bounds_error=False, fill_value=np.nan)
        aligned_power[i, :] = f(elevation_range)

    if elevation_limit is not None:
        in_limit_flag = np.logical_and(elevation_range >= elevation_limit[0],
                                       elevation_range <= elevation_limit[1])
        in_limit = np.where(in_limit_flag)[0]
        elevation_range = elevation_range[in_limit]
        aligned_power = aligned_power[:, in_limit]

    return aligned_power, elevation_range
