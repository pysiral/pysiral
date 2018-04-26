# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 15:41:30 2016

@author: shendric
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

import numpy as np

from pysiral.visualization.mapstyle import custom_font


class DataPlot(object):

    def __init__(self):
        self.filename = None
        self.dpi = 300
        self.l1b = None
        self.fig = None
        self.fig_args = {"figsize": (10, 16), "facecolor": "white"}
        set_mpl_default_style()

    def savefig(self, filename, close_figure=True):
        self.filename = filename
        plt.savefig(self.filename, dpi=self.dpi,
                    facecolor=self.fig.get_facecolor())
        if close_figure:
            plt.close(self.fig)

    @property
    def canvas_aspect(self):
        figsize = self.fig_args["figsize"]
        return float(figsize[1])/float(figsize[0])


def set_mpl_default_style():
    mpl.rcParams['font.sans-serif'] = "arial"
    for target in ["xtick.color", "ytick.color", "axes.edgecolor",
                   "axes.labelcolor"]:
        mpl.rcParams[target] = "#4b4b4d"


class ParameterScatterPlot(object):

    def __init__(self, x, y, boxprop_dict, title=None):
        
        # Save parameters
        self.x = x
        self.y = y
        self.boxprop_dict = boxprop_dict
        self.title = title or ""

        # Create figures
        self._init_figure()
        self._add_content()
        self._format_figure()

    def savefig(self, filename, dpi=300, close=True):
        plt.savefig(filename, dpi=dpi)
        if close:
            plt.close(self.figure)

    def _init_figure(self):
        plt.ioff()
        self.figure = plt.figure(figsize=(6, 6))
        self.ax = plt.gca()
        self.ax.set_position([0.15, 0.12, 0.8, 0.8])

    def _add_content(self):
 
        # Create hexbin plot of point density in background
        hb = self.ax.hexbin(self.x, self.y, gridsize=200, cmap=plt.cm.GnBu, edgecolors="white", linewidths=0.2, zorder=200)
        bin_max_val = np.max(hb.get_array())
        hb.set(norm=Normalize(vmin=0, vmax=1.5*bin_max_val))
        
        # compute boxplot dataset
        binsize = 0.2
        positions = np.arange(0.0, 4.0+binsize/2.0, binsize)
        dataset = []
        n_points = []
        for pos in positions:
            ll, ul = pos-binsize/2., pos+binsize/2.
            indices = np.where(np.logical_and(self.x>=ll, self.x<ul))
            dataset.append(self.y[indices])
            n_points.append(len(indices[0]))

        n_points = list(np.cumsum(n_points))
        min_val, max_val = np.percentile(n_points, [5, 95])
        perc5, perc95 = positions[n_points.index(min_val)], positions[n_points.index(max_val)]

        bp = self.ax.boxplot(dataset, positions=positions, widths=0.6*binsize, zorder=205, showfliers=False,
                             whis=[5, 95], 
                             manage_xticks=False, boxprops=dict(linewidth=0.75),
                             whiskerprops=dict(linewidth=0.75), capprops=dict(linewidth=0.75),
                             meanprops=dict(color="0.0", linestyle="-"), showmeans=True, meanline=True,
                             medianprops=dict(linewidth=0.0))


        for target in ["boxes", "means"]:
            for i, obj in enumerate(bp[target]):
                if positions[i] <= perc5 or positions[i] >= perc95:
                    obj.set_color("#B0C4DE")
        for target in ["whiskers", "caps"]:
            for i, obj in enumerate(bp[target]):
                if positions[i/2] <= perc5 or positions[i/2] >= perc95:
                    obj.set_color("#B0C4DE")

        self.ax.plot([-1, 6], [-1, 6], linewidth=0.5, linestyle="--", dashes=(5, 5), color="#4b4b4d", zorder=201)
        # self.ax.plot([-0.9, 1.1], [-1, 1], linewidth=0.75, linestyle=":", color="#4b4b4d")
        # self.ax.plot([-1.1, .9], [-1, 1], linewidth=0.75,linestyle=":", color="#4b4b4d")
        self.ax.axvline(linewidth=0.5, color="#4b4b4d", zorder=201)
        self.ax.axhline(linewidth=0.5, color="#4b4b4d", zorder=201)

    def _format_figure(self):

        # Set period as title
        self.ax.set_title(self.title, **custom_font(18))
        self.ax.set_xlim(-0.2, 5.2)
        self.ax.set_ylim(-0.2, 5.2)
        ticks = np.linspace(0.0, 5, 6)
        self.ax.xaxis.set_ticks(ticks)
        self.ax.xaxis.set_ticklabels(["%.1f" % tick for tick in ticks])
        self.ax.yaxis.set_ticks(ticks)
        self.ax.set_xlabel("Envisat Sea Ice Thickness (m)", **custom_font(12))
        self.ax.set_ylabel("CryoSat-2 Sea Ice Thickness (m)", **custom_font(12))

        # Remove mirror axis (Exception first 2d hist plot)
        x_off = -0.02
        y_off = -0.02
        spine_props = dict(length=8, which='major', pad=8)

        spines_to_remove = ["top", "right"]
        for spine in spines_to_remove:
            self.ax.spines[spine].set_visible(False)

        self.ax.spines["bottom"].set_position(("axes", x_off))
        self.ax.tick_params('bottom', **spine_props)
        self.ax.spines["left"].set_position(("axes", y_off))
        self.ax.tick_params('bottom', **spine_props)

        cl = plt.getp(self.ax, 'xmajorticklabels')
        plt.setp(cl, **custom_font(12))

        cl = plt.getp(self.ax, 'ymajorticklabels')
        plt.setp(cl, **custom_font(12))



    @property
    def data_bins(self):
        return get_nphistbins([-0.2, 0.5], 0.005)
