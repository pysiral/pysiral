# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 12:42:13 2016

@author: shendric
"""

import argparse
import glob
import os
import numpy as np
from datetime import datetime
from calendar import monthrange

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import FixedLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from pysiral.config import ConfigInfo
from pysiral.logging import DefaultLoggingClass
from pysiral.path import filename_from_path, validate_directory
from pysiral.output import PysiralOutputFilenaming


def l1bdata_database_summary_graph():

    # get argument parser
    parser = get_argparser()
    args = parser.parse_args()
    job = L1bdataDatabaseSummary(args)
    job.execute()


class L1bdataDatabaseSummary(DefaultLoggingClass):

    def __init__(self, options):

        super(L1bdataDatabaseSummary, self).__init__(self.__class__.__name__)
        self.options = options
        self.l1brepo = L1bdataRepo()
        self.timestamp = datetime.now()

    def execute(self):
        """
        Loop over available l1bdata folders and create yearly
        summary plots
        """

        # Get a summary of all files
        self.l1brepo.set_l1bdata_folder(self.options.folder)
        self.l1brepo.set_l1bdata_version(self.options.version)
        self.l1brepo.get_l1bdata_inventory()

        # Create output plots
        plot = L1bdataSummaryPlot()
        plot.set_inventory(self.l1brepo.inventory)
        plot.add_metadata(mission_id=self.l1brepo.mission_id,
                          timestamp=self.timestamp,
                          version=self.options.version)
        for year in self.l1brepo.years:
            plot.create_yearly_graph(year)
            plot.save_graph(self.get_graph_output_directory(year))

    def get_graph_output_directory(self, year):
        output_directory = os.path.join(
            self.options.folder,
            "database_inventory",
            self.options.version)
        validate_directory(output_directory)
        return output_directory


class L1bdataRepo(DefaultLoggingClass):

    def __init__(self):
        super(L1bdataRepo, self).__init__(self.__class__.__name__)
        self.l1bdata_folder = None
        self.l1bdata_version = "*"
        self.mission_id = None
        self.years = []
        self.inventory = None
        self.has_hemisphere = {"north": False, "south": False}

    def set_l1bdata_folder(self, folder):
        self.l1bdata_folder = folder
        self.log.info("Check Folder: %s" % folder)
        # Verify that north & south exists
        for hemisphere in self.hemispheres:
            self.has_hemisphere[hemisphere] = os.path.isdir(
                os.path.join(folder, hemisphere))
            self.log.info("%s subfolder exists: %s" % (
                hemisphere, str(self.has_hemisphere[hemisphere])))

    def set_l1bdata_version(self, version):
        self.l1bdata_version = version

    def get_l1bdata_inventory(self):
        """ Count l1bdata files for each day in l1bdata repository """
        self.get_years()
        self.init_inventory()
        for hemisphere in self.existing_hemispheres:
            self.count_l1bdata_files(hemisphere)

    def get_years(self):
        """ Get the list of north/south years """
        dirs = []
        for hemisphere in self.existing_hemispheres:
            dirs.extend(
                os.listdir(os.path.join(self.l1bdata_folder, hemisphere)))
        dirs = np.unique(dirs)
        self.years = [int(year) for year in dirs if self.is_year_folder(year)]

    def is_year_folder(self, folder_name):
        """ Check if folder name is a valid year """
        is_year = False
        try:
            int_year = int(folder_name)
        except:
            return is_year
        if int_year > 1991 and int_year <= datetime.now().year:
            is_year = True
        return is_year

    def init_inventory(self):
        """ Counter for each day in both hemisphere """
        ref = np.full((31), 0, dtype=np.int16)
        self.inventory = dict()
        for hemisphere in self.hemispheres:
            self.inventory[hemisphere] = dict()
            for year in self.years:
                self.inventory[hemisphere][year] = dict()
                for month in np.arange(12):
                    self.inventory[hemisphere][year][month+1] = np.array(ref)

    def count_l1bdata_files(self, hemisphere):
        """ Loop over all l1bdata files and count segments per day """
        basedir = os.path.join(self.l1bdata_folder, hemisphere)
        l1binfo = PysiralOutputFilenaming()
        # pysiral organizes data in the structury yyyy/mm
        for year in self.years:
            for month in np.arange(1, 13):
                # Get full diretory
                l1bdir = os.path.join(basedir, "%04g" % year, "%02g" % month)
                # Folder might not exist
                if not os.path.isdir(l1bdir):
                    continue
                # Get list of l1bdata files
                version = self.l1bdata_version
                search = os.path.join(l1bdir, "l1bdata_%s_*.nc" % version)
                l1bdata_files = glob.glob(search)
                self.log.info("%s %04g-%02g [%g files]" % (
                    hemisphere, year, month, len(l1bdata_files)))
                # Parse date from filename and increase counter
                for l1bdata_file in l1bdata_files:
                    filename = filename_from_path(l1bdata_file)
                    l1binfo.parse_filename(filename)
                    day = l1binfo.start.day
                    self.mission_id = l1binfo.mission_id
                    self.inventory[hemisphere][year][month][day-1] += 1

    @property
    def hemispheres(self):
        return ["north", "south"]

    @property
    def existing_hemispheres(self):
        return [hemisphere for hemisphere in self.hemispheres
                if self.has_hemisphere[hemisphere]]


class L1bdataSummaryPlot(DefaultLoggingClass):

    def __init__(self):
        super(L1bdataSummaryPlot, self).__init__(self.__class__.__name__)
        self.l1bdata_inventory = None
        self.figure = None
        self.output_filename = None
        self.axarr = []
        self.timestamp = None
        self.version = None
        self.mission_id = None
        self.cmap_props = {"vmin": 1, "vmax": 50}

    def set_inventory(self, l1bdata_inventory):
        self.l1bdata_inventory = l1bdata_inventory

    def add_metadata(self, **props):
        for key in props.keys():
            setattr(self, key, props[key])

    def create_yearly_graph(self, year):
        self._setup_figure()
        self._create_plot(year)
        self._add_colorbar()
        self._annotate_graph(year)
        self._set_output_filename(year)

    def save_graph(self, output_directory):
        output = os.path.join(output_directory, self.output_filename)
        plt.savefig(output, dpi=600, facecolor=self.figure.get_facecolor())
        plt.close(self.figure)

    def _setup_figure(self):
        mpl.rcParams['font.sans-serif'] = "arial"
        targets = ["xtick.color", "ytick.color", "axes.edgecolor",
                   "axes.labelcolor"]
        for target in targets:
            mpl.rcParams[target] = "#4b4b4d"
        subplot_props = {"figsize": (18, 5), "sharex": True,
                         "facecolor": "1.0"}
        self.figure, self.axarr = plt.subplots(1, 2, **subplot_props)
        plt.subplots_adjust(top=0.98, bottom=0.15, wspace=0.1, left=0.05,
                            right=0.95)

    def _create_plot(self, year):
        ax_indices = {"north":  0, "south": 1}
        yaxis_position = {"north": "left", "south": "right"}
        title = {"north": "Arctic", "south": "Antarctic"}

        for hemisphere in ["north", "south"]:
            counts = self._get_pcolor_array(hemisphere, year)
            ax = self.axarr[ax_indices[hemisphere]]
            ax.pcolor(counts, color="1.0", lw=2,
                      cmap=self.cmap, **self.cmap_props)
            ax.set_xlim(0, 31)
            ax.set_xticks(self.xticks)
            ax.set_xticklabels(self.xticklabels)
            ax.set_yticks(self.yticks)
            ax.set_yticklabels(self.yticklabels)
            ax.set_title(title[hemisphere])
            ax.yaxis.set_ticks_position(yaxis_position[hemisphere])
            spines_to_remove = ["top", "bottom", "left", "right"]
            ax.set(adjustable='box-forced', aspect='equal')
            ax.tick_params(axis=u'both', which=u'both', length=0)
            for spine in spines_to_remove:
                ax.spines[spine].set_visible(False)

    def _add_colorbar(self):
        sm = plt.cm.ScalarMappable(cmap=self.cmap,
                                   norm=plt.Normalize(**self.cmap_props))
        ax = self.axarr[0]
        sm._A = []
        cb_ax_kwargs = {
            'loc': 3, 'bbox_to_anchor': (0.0, -0.2, 1, 1),
            'width': "30%", 'height': "5%", 'bbox_transform': ax.transAxes,
            'borderpad': 0}
        ticks = FixedLocator([self.cmap_props["vmin"],
                              self.cmap_props["vmax"]])
        axins = inset_axes(ax, **cb_ax_kwargs)
        cb = plt.colorbar(sm, cax=axins, ticks=ticks,
                          orientation="horizontal", extend="both")
        cl = plt.getp(cb.ax, 'xmajorticklabels')
        plt.setp(cl, fontsize="small")
        cb.set_label("# L1b Segments / day", fontsize="small")
        cb.outline.set_linewidth(0.2)
        cb.outline.set_alpha(0.0)
        for t in cb.ax.get_yticklines():
            t.set_color("1.0")
        cb.ax.tick_params('both', length=0.1, which='major', pad=10)
        plt.sca(ax)

    def _annotate_graph(self, year):

        plt.annotate("%04g" % year, (0.01, 0.92), xycoords="figure fraction",
                     fontsize="xx-large")

        timestamp_str = "{dt:%Y-%m-%d %H:%M}".format(dt=self.timestamp)
        created_str = "Created on: %s" % timestamp_str
        plt.annotate(created_str, (0.995, 0.02), xycoords="figure fraction",
                     fontsize="x-small", ha="right")

        info = ConfigInfo()
        mission_name = info.get_mission_info(self.mission_id).long_name
        version_str = "%s (version:%s)" % (mission_name, self.version)
        plt.annotate(version_str, (0.2, 0.145), xycoords="figure fraction",
                     fontsize="larger")

    def _get_pcolor_array(self, hemisphere, year):
        inventory = self.l1bdata_inventory
        counts = np.ma.array(np.zeros((12, 31)))
        mask = np.ones((12, 31), dtype=np.bool)
        for month in np.arange(1, 13):
            counts[12-month, :] = inventory[hemisphere][year][month][:]
            days = monthrange(year, month)
            days_list = np.arange(days[1])
            mask[12-month, days_list] = False
        counts.mask = mask
        return counts

    @property
    def cmap(self):
        cdict1 = {'red':  ((0.0, 1.0, 1.0),
                           (0.25, 1.0, 1.0),
                           (1.0, 0.5, 0.5)),
                  'green': ((0.0, 0.38, 0.38),
                            (0.25, 0.84, 0.84),
                            (1.0, 1.0, 1.0)),
                  'blue':  ((0.0, 0.28, 0.28),
                            (0.25, 0.0, 0.0),
                            (1.0, 0.0, 0.0))}
        cmap = LinearSegmentedColormap('quality', cdict1)
        cmap.set_under('0.9')
        cmap.set_over('#e066ff')
        return cmap

    @property
    def xticks(self):
        return np.arange(31)+0.5

    @property
    def xticklabels(self):
        return ["%02g" % (intday+1) for intday in np.arange(31)]

    @property
    def yticks(self):
        return np.arange(12)+0.5

    @property
    def yticklabels(self):
        return ["Dec", "Nov", "Oct", "Sep", "Aug", "Jul", "Jun", "May",
                "Apr", "Mar", "Feb", "Jan"]

    def _set_output_filename(self, year):
        self.output_filename = "l1bdata_database_inventory_%04g.png" % year


def get_argparser():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-version",
        required=True,
        action='store',
        dest='version',
        type=str,
        help='l1bdata data version')

    # Main product folder
    parser.add_argument(
        action='store',
        dest='folder',
        help='l1bdata product folder')

    return parser


if __name__ == "__main__":
    l1bdata_database_summary_graph()
