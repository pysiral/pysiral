# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 15:41:30 2016

@author: shendric
"""

import matplotlib as mpl
import matplotlib.pyplot as plt


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
