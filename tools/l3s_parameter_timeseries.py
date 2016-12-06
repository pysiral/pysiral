# -*- coding: utf-8 -*-
"""
Created on Tue Dec 06 09:50:05 2016

@author: shendric
"""

from pysiral.config import ConfigInfo
from pysiral.logging import DefaultLoggingClass
from pysiral.iotools import NCMaskedGridData
from pysiral.path import validate_directory
from pysiral.visualization.mapstyle import get_custom_font

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import date2num
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

from datetime import datetime
import argparse
import glob
import os


class L3sParameterTimeSeries(DefaultLoggingClass):

    def __init__(self):
        super(L3sParameterTimeSeries, self).__init__(self.__class__.__name__)
        self.latlonbox = LatLonBox(area="world")
        self.output_folder = "auto"
        self.parameter_list = []
        self.l3s_stack = L3sSourceList()
        self.l4plotdef = Level4ParameterPlotSettings()

    def init_from_args(self):
        """
        Initialize the job from command line arguments
        XXX: There should be other ways of doing this
        """

        # Parse command line args
        self.args = self.arg_parser.parse_args()

        # Update latlonbox
        if self.args.latlonbox is not None:
            self.latlonbox.set_from_values(*self.args.latlonbox)
        if self.args.roi_name is not None:
            self.latlonbox.set_region_name(self.args.roi_name)

        # Update plot settings
        self.l4plotdef.update(labels=self.args.plot_label,
                              annotation=self.args.plot_annotation,
                              export_map=self.args.export_map)

        # Set parameter_list
        self.parameter_list = self.args.parameter_list

        # Get output folder
        # XXX: Automatic folder generation assumes that lowest level folder
        #      is 'l3s'
        if self.args.output_folder != "auto":
            self.output_folder = self.args.output_folder
        else:
            dirs = os.path.split(self.args.l3s_folders[0])
            self.output_folder = os.path.join(os.sep.join(dirs[0:-1]), "l4")
        self.log.info("Output Folder: %s" % self.output_folder)

    def execute(self):

        # Get all files
        self.l3s_stack.get_files(self.args.l3s_folders)

        # Loop over all Parameter
        for parameter_name in self.parameter_list:

            # Extract data (from ROI)
            self.log.info("Extracting parameter:%s in region:%s" % (
                          parameter_name, self.latlonbox.full_name))
            l4par = Level4MultiMissionParameter()
            l4par.set_roi(self.latlonbox)
            l4par.get_from_l3s(self.l3s_stack.list, parameter_name)

            # Create output plots
            self.l4plotdef.set_parameter_name(parameter_name)
            self.l4plotdef.set_region_name(l4par.region_name)
            l4par_plot = Level4ParameterPlot(self.l4plotdef)
            l4par_plot.set_parameter_data(l4par)
            l4par_plot.export(self.output_folder)

    @property
    def arg_parser(self):
        """ Return the Argument parser for l3s_time_series """

        parser = argparse.ArgumentParser()

        parser.add_argument(action='store', dest='l3s_folders',
                            nargs="+", help='l3s product folder')

        parser.add_argument("-label", action='store', dest='plot_label',
                            type=str, required=False,
                            help='Plot Label')

        parser.add_argument("-annotation", action='store',
                            dest='plot_annotation', type=str, required=False,
                            help='Plot Label')

        parser.add_argument("-latlonbox", action='store', dest='latlonbox',
                            type=float, nargs=4, required=False,
                            help='lon_min, lat_min, lon_max, lat_max')

        parser.add_argument("-roi-name", action='store', dest='roi_name',
                            type=str, required=False, default=None,
                            help='label of latlonbox')

        parser.add_argument("-output-folder", action='store',
                            dest='output_folder', type=str, default="auto",
                            required=False, help='folder for writing results')

        parser.add_argument("--export-map", action='store_true',
                            dest='export_map', default=False,
                            help='Export map with lat/lon Box ')

        parser.add_argument('-parameter', action='store',
                            dest='parameter_list', nargs='+', required=True,
                            help='list of parameter')

        return parser


class L3sSourceList(DefaultLoggingClass):

    def __init__(self):
        super(L3sSourceList, self).__init__(self.__class__.__name__)
        self._l3s_file_list = []

    def get_files(self, folders):

        for folder in folders:
            self._l3s_recursive_search(folder)
        self._l3s_file_list = sorted(self._l3s_file_list)

    def _l3s_recursive_search(self, folder):
        self.log.info("Search L3s netCDF files in %s" % folder)
        for root, dirs, files in os.walk(folder):
            for directory in dirs:
                l3s_query = os.path.join(root, directory, self.search_pattern)
                l3s_files = glob.glob(l3s_query)
                n_l3s_files = len(l3s_files)
                if n_l3s_files > 0:
                    self._l3s_file_list.extend(l3s_files)
                    self.log.info("- Subfolder %s: found %g" % (
                                  directory, n_l3s_files))

    @property
    def search_pattern(self):
        return "L3S*.nc"

    @property
    def list(self):
        return self._l3s_file_list


class Level4MultiMissionParameter(DefaultLoggingClass):

    def __init__(self):
        class_name = self.__class__.__name__
        super(Level4MultiMissionParameter, self).__init__(class_name)
        self.roi = None
        self.period = None
        self.region_name = None
        self.l4_points = []

    def set_roi(self, roi):
        self.roi = roi
        self.region_name = roi.full_name

    def get_from_l3s(self, l3s_filenames, parameter_name):
        """ Extract parameter in roi for each l3s file """

        # loop over files
        for l3s_filename in l3s_filenames:

            # Read data
            l3s = NCMaskedGridData(l3s_filename)

            # Compute the mask (may be different for each l3s file in stack)
            # If no roi is specified, use entire grid
            if self.roi is not None:
                roi_mask = self.roi.get_grid_mask(l3s.longitude, l3s.latitude)
            else:
                roi_mask = np.full(l3s.longitude.shape, True)
                self.region_name = l3s.hemisphere

            # Extract the parameter
            parameter_grid = l3s.get_by_name(parameter_name)
            roi_indices = np.where(roi_mask)
            parameter = parameter_grid[roi_indices]

            # Create a L4 data point
            l4dat = Level4DataPoint()
            # XXX: Hard coded to monthly data
            period = (l3s.start_time.year, l3s.start_time.month)
            l4dat.set_metadata(mission=l3s.mission_ids, period=period,
                               parameter_name=parameter_name)
            l4dat.set_parameter_data(parameter)
            self.l4_points.append(l4dat)

    def get_data_range(self, paddays=0):
        dates = self.dates
        return([np.amin(dates)-paddays, np.amax(dates)+paddays])

    def get_boxplot_data(self, mission):
        data = [l4.array for l4 in self.l4_points if l4.mission == mission]
        dates = [l4.datenum for l4 in self.l4_points if l4.mission == mission]
        if self.n_missions > 1:
            mission_index = list(self.missions).index(mission)
            offset = -5 + mission_index*10
        else:
            offset = 0
        return data, np.array(dates)+offset

    def get_all_periods(self):
        return self.dates, self.labels

    @property
    def dates(self):
        return [l4.datenum for l4 in self.l4_points]

    @property
    def labels(self):
        return [l4.label for l4 in self.l4_points]

    @property
    def years(self):
        years_list = [l4.period[0] for l4 in self.l4_points]
        return np.unique(sorted(years_list))

    @property
    def missions(self):
        sort_order = ["ers1", "ers2", "envisat", "cryosat2"]
        mission_list = [l4.mission for l4 in self.l4_points]
        all_missions = np.unique(sorted(mission_list))
        missions = []
        for mission in sort_order:
            if mission in all_missions:
                missions.append(mission)
        return missions

    @property
    def n_missions(self):
        return len(self.missions)


class Level4DataPoint(object):

    def __init__(self):
        self.parameter_name = None
        self.array = None
        self.period = None
        self.mean = None
        self.sdev = None
        self.median = None
        self.min = None
        self.max = None

    def set_metadata(self, **keyws):
        for key in keyws.keys():
            setattr(self, key, keyws[key])

    def set_parameter_data(self, data):
        no_nans = np.where(np.isfinite(data))
        self.array = data[no_nans].ravel()

    @property
    def datenum(self):
        return date2num(self.date)

    @property
    def label(self):
        return self.date.strftime("%b")

    @property
    def date(self):
        return datetime(self.period[0], self.period[1], 15)


class LatLonBox(DefaultLoggingClass):

    def __init__(self, area=None):
        super(LatLonBox, self).__init__(self.__class__.__name__)
        self.lonbounds = [None, None]
        self.latbounds = [None, None]
        self._area_name = None
        if area is not None:
            self.set_from_area_name(area)

    def set_from_values(self, lon_min, lat_min, lon_max, lat_max):
        self.lonbounds = [lon_min, lon_max]
        self.latbounds = [lat_min, lat_max]

    def set_from_area_name(self, area):
        if area == "world":
            self.lonbounds = [-180, 180]
            self.latbounds = [-90.0, 90.0]
            self._area_name = "Global"
        elif area == "north":
            self.lonbounds = [-180, 180]
            self.latbounds = [50.0, 90.0]
            self._area_name = "Arctic"
        elif area == "south":
            self.lonbounds = [0.0, 360]
            self.latbounds = [-90.0, -50.0]
            self._area_name = "Antarctic"
        else:
            self.log.warn("Unknown Area: %s" % str(area))

    def set_region_name(self, name):
        self._area_name = name

    def get_grid_mask(self, lon, lat):

        # Sanity check (currently longitude range -180 - 180 assumed)
        if len(np.where(lon.ravel() > 180.0)[0]) > 0:
            raise ValueError("Longitude not in -180 to 180 range")

        # Simple longitude/latitude mask
        in_lon_range = np.logical_and(lon >= self.lonbounds[0],
                                      lon <= self.lonbounds[1])
        in_lat_range = np.logical_and(lat >= self.latbounds[0],
                                      lat <= self.latbounds[1])
        return np.logical_and(in_lon_range, in_lat_range)

    @property
    def full_name(self):
        return self._area_name


class Level4ParameterPlotSettings(DefaultLoggingClass):

    def __init__(self):
        class_name = self.__class__.__name__
        super(Level4ParameterPlotSettings, self).__init__(class_name)
        self.label = None
        self.annotation = None
        self.plot_range = None
        self.export_map = False
        self.region_name = None

    def update(self, **kwargs):
        pass

    def set_region_name(self, region_name):
        self.region_name = region_name

    def set_parameter_name(self, parameter_name):
        self.parameter_name = parameter_name


class Level4ParameterPlot(DefaultLoggingClass):

    def __init__(self, plotdef):
        super(Level4ParameterPlot, self).__init__(self.__class__.__name__)
        self.plotdef = plotdef
        self.l4dat = None

    def set_parameter_data(self, l4dat):
        self.l4 = l4dat

    def export(self, output_folder):
        self._init_figure()
        self._create_plot()
        self._export(output_folder)

    def _init_figure(self):
        self.figure = plt.figure(figsize=(16, 6))
        self.ax = plt.gca()
        self.ax.set_position([0.05, 0.15, 0.92, 0.80])
        self.ax.set_axis_bgcolor('0.98')

    def _create_plot(self):
        pysiral_config = ConfigInfo()
        colors = ["#00ace5", "#003e6e"]
        for i, mission in enumerate(self.l4.missions):

            flierprops = dict(marker=',', color=colors[i], markersize=1)
            boxprops = dict(color=colors[i], linewidth=1)
            whiskerprops = dict(linestyle="-", color=colors[i], linewidth=0.75)
            capprops = whiskerprops
            meanprops = dict(zorder=101, mec="1.0", markerfacecolor=colors[i])
            medianprops = dict(linewidth=0)

            data, dates = self.l4.get_boxplot_data(mission)
            bp = plt.boxplot(data, positions=dates,
                             flierprops=flierprops, whiskerprops=whiskerprops,
                             capprops=capprops, boxprops=boxprops,
                             meanprops=meanprops, medianprops=medianprops,
                             showmeans=True, widths=5,
                             showfliers=False)

            # Draw boxes as filled objects
            for box in bp["boxes"]:
                xdata, ydata = box.get_xdata(), box.get_ydata()
                xmin, xmax = np.amin(xdata), np.amax(xdata)
                ymin, ymax = np.amin(ydata), np.amax(ydata)
                rect = Rectangle([xmin, ymin], xmax-xmin, ymax-ymin,
                                 ec="none", fc=colors[i], zorder=100)
                self.ax.add_patch(rect)

            fontprops = {
                "color": colors[i],
                "fontproperties": get_custom_font(fontsize=16, awi_font=True)}

            mission_info = pysiral_config.get_mission_info(mission)
            labelx = 0.05 + i*0.075
            plt.annotate(mission_info.long_name, (labelx, 0.05),
                         xycoords="figure fraction",
                         **fontprops)

        # Force correct xticks
        xticks, xticklabels = self.l4.get_all_periods()
        self.ax.set_xticks(xticks)
        self.ax.set_xticklabels(xticklabels)

        fontprops = {
            "color": "#4b4b4d",
            "fontproperties": get_custom_font(fontsize=12, awi_font=True)}

        patches = []
        for year in self.l4.years:
            labelx = date2num(datetime(year, 7, 1))
            plt.annotate(str(year), (labelx, 3.8), xycoords="data",
                         ha="center", **fontprops)
            if np.mod(year, 2) == 0:
                continue
            xmin = date2num(datetime(year, 1, 1))
            ymin = -0.5
            width = 365
            height = 4.5
            rect = Rectangle([xmin, ymin], width, height, ec="none", fc="0.94")
            self.ax.add_patch(rect)

        plt.xlim(self.l4.get_data_range(paddays=30))
        plt.ylim(-0.5, 4)

        spines_to_remove = ["top", "right"]
        for spine in spines_to_remove:
            self.ax.spines[spine].set_visible(False)

        self.ax.yaxis.set_tick_params(direction='out')
        self.ax.yaxis.set_ticks_position('left')
        self.ax.spines["left"].set_position(("axes", -0.01*6./16.))

        self.ax.xaxis.set_tick_params(direction='out')
        self.ax.xaxis.set_ticks_position('bottom')
        self.ax.spines["bottom"].set_position(("axes", -0.01))

        fontprops = {
            "color": "#4b4b4d",
            "fontproperties": get_custom_font(fontsize=5, awi_font=True)}
        cl = plt.getp(self.ax, 'xmajorticklabels')
        plt.setp(cl, **fontprops)

        fontprops = {
            "color": "#4b4b4d",
            "fontproperties": get_custom_font(fontsize=12, awi_font=True)}

        cl = plt.getp(self.ax, 'ymajorticklabels')
        plt.setp(cl, **fontprops)
        self.ax.set_ylabel("Sea Ice Thickness (meter)", **fontprops)

    def _export(self, output_folder):
        validate_directory(output_folder)
        filename = os.path.join(output_folder, "test.png")
        plt.savefig(filename, dpi=600)
        plt.close(self.figure)


if __name__ == "__main__":
    job = L3sParameterTimeSeries()
    job.init_from_args()
    job.execute()
