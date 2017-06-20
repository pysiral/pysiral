# -*- coding: utf-8 -*-
"""
Created on Tue Dec 06 09:50:05 2016

@author: shendric
"""

from pysiral.config import ConfigInfo
from pysiral.logging import DefaultLoggingClass
from pysiral.iotools import NCMaskedGridData
from pysiral.path import validate_directory, file_basename
from pysiral.visualization.mapstyle import get_custom_font
from pysiral.visualization.mapstyle import GridMapPaperStyle
from pysiral.maptools import get_landcoastlines

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import date2num
from matplotlib.patches import Rectangle
from mpl_toolkits.basemap import Basemap

from datetime import datetime
import argparse
import glob
import os


class L3ParameterTimeSeries(DefaultLoggingClass):

    def __init__(self):
        super(L3ParameterTimeSeries, self).__init__(self.__class__.__name__)
        self.latlonbox = LatLonBox(area="world")
        self.output_folder = "auto"
        self.parameter_list = []
        self.l4roimap = Level4ROIMap()
        self.l3_stack = L3SourceList()
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
                              custom_plot_range=self.args.plot_range,
                              export_map=self.args.export_map)

        # Update map plot
        self.l4roimap.is_active = self.args.export_map
        self.l4roimap.annotation = self.args.plot_annotation
        if self.args.roi_name is not None:
            self.l4roimap.region_name = self.args.roi_name
        else:
            self.l4roimap.region_name = "global"

        # Set parameter_list
        self.parameter_list = self.args.parameter_list

        # Get output folder
        # XXX: Automatic folder generation assumes that lowest level folder
        #      is 'l3'
        if self.args.output_folder != "auto":
            self.output_folder = self.args.output_folder
        else:
            dirs = os.path.split(self.args.l3_folders[0])
            self.output_folder = os.path.join(os.sep.join(dirs[0:-1]), "l4")
        self.log.info("Output Folder: %s" % self.output_folder)

    def execute(self):

        # Get all files
        self.l3_stack.get_files(self.args.l3_folders)

        # Loop over all Parameter
        for parameter_name in self.parameter_list:

            # Extract data (from ROI)
            self.log.info("Extracting parameter:%s in region:%s" % (
                          parameter_name, self.latlonbox.full_name))
            l4par = Level4MultiMissionParameter()
            l4par.set_roi(self.latlonbox)
            l4par.get_from_l3(self.l3_stack.list, parameter_name)

            # Create output plots
            self.l4plotdef.set_parameter_name(parameter_name)
            self.l4plotdef.set_region_name(l4par.region_name)
            l4par_plot = Level4ParameterPlot(self.l4plotdef)
            l4par_plot.set_parameter_data(l4par)
            l4par_plot.export(self.output_folder)

        # Export one map
        if self.l4roimap.is_active:
            self.l4roimap.set_roi(self.latlonbox)
            self.l4roimap.export(self.output_folder)

    @property
    def arg_parser(self):
        """ Return the Argument parser for l3_time_series """

        parser = argparse.ArgumentParser()

        parser.add_argument(action='store', dest='l3_folders',
                            nargs="+", help='l3 product folder')

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

        parser.add_argument('-plot-range', action='store',
                            dest="plot_range", nargs=4, type=int,
                            help='plot range: yyyy mm till yyyy mm')
        return parser


class L3SourceList(DefaultLoggingClass):

    def __init__(self):
        super(L3SourceList, self).__init__(self.__class__.__name__)
        self._l3_file_list = []

    def get_files(self, folders):

        for folder in folders:
            self._l3_recursive_search(folder)
        self._l3_file_list = sorted(self._l3_file_list)

    def _l3_recursive_search(self, folder):
        self.log.info("Search L3 netCDF files in %s" % folder)
        for root, dirs, files in os.walk(folder):
            for directory in dirs:
                l3_query = os.path.join(root, directory, self.search_pattern)
                l3_files = glob.glob(l3_query)
                n_l3_files = len(l3_files)
                if n_l3_files > 0:
                    self._l3_file_list.extend(l3_files)
                self.log.info("- Subfolder %s: found %g" % (
                              directory, n_l3_files))

    @property
    def search_pattern(self):
        return "*L3C*.nc"

    @property
    def list(self):
        return self._l3_file_list


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

    def get_from_l3(self, l3_filenames, parameter_name):
        """ Extract parameter in roi for each l3 file """

        # loop over files
        for l3_filename in l3_filenames:

            # Read data
            l3 = NCMaskedGridData(l3_filename)

            # Compute the mask (may be different for each l3 file in stack)
            # If no roi is specified, use entire grid
            if self.roi is not None:
                roi_mask = self.roi.get_grid_mask(l3.longitude, l3.latitude)
            else:
                roi_mask = np.full(l3.longitude.shape, True)
                self.region_name = l3.hemisphere

            # Extract the parameter
            parameter_grid = l3.get_by_name(parameter_name)
            roi_indices = np.where(roi_mask)
            parameter = parameter_grid[roi_indices]

            # Create a L4 data point
            l4dat = Level4DataPoint()
            # XXX: Hard coded to monthly data
            # start_time = dateutil.parser.parse(l3.time_coverage_start)
            # period = (start_time.year, start_time.month)
            basename = file_basename(l3_filename)
            strarr = basename.split("-")
            year = int(strarr[7][0:4])
            month = int(strarr[7][4:6])
            period = (year, month)
            # XXX: Hotfix for SICCI2 L3C v0.9
            sensor_catalog = {"RA2": "envisat", "SIRAL": "cryosat2"}
            mission_id = sensor_catalog[l3.sensor]
            l4dat.set_metadata(mission=mission_id, period=period,
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

    def get_boundary_line(self, n_points_per_line=100):
        lons = np.array([])
        lats = np.array([])
        lon_indices = [(0, 1,), (1, 1), (1, 0), (0, 0)]
        lat_indices = [(0, 0), (0, 1), (1, 1), (1, 0)]
        for i in np.arange(4):
            lon_index, lat_index = lon_indices[i], lat_indices[i]
            lon_line = np.linspace(self.lonbounds[lon_index[0]],
                                   self.lonbounds[lon_index[1]],
                                   n_points_per_line)
            lat_line = np.linspace(self.latbounds[lat_index[0]],
                                   self.latbounds[lat_index[1]],
                                   n_points_per_line)
            lons = np.append(lons, lon_line)
            lats = np.append(lats, lat_line)
        return lons, lats

    @property
    def full_name(self):
        return self._area_name


class Level4ParameterPlotSettings(DefaultLoggingClass):
    """ Contains all settings for the plot """

    def __init__(self):
        class_name = self.__class__.__name__
        super(Level4ParameterPlotSettings, self).__init__(class_name)
        self.label = None
        self.annotation = None
        self.plot_range = None
        self.export_map = False
        self.region_name = None
        self.custom_plot_range = None

    def update(self, **kwargs):
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])

    def set_region_name(self, region_name):
        self.region_name = region_name

    def set_parameter_name(self, parameter_name):
        self.parameter_name = parameter_name

    def get_boxplotprops(self, mission_index):
        color = self.mission_colors[mission_index]
        boxplot_props = dict(
            flierprops=dict(marker=',', color=color, markersize=1),
            boxprops=dict(color=color, linewidth=1),
            whiskerprops=dict(linestyle="-", color=color, linewidth=0.75),
            capprops=dict(linestyle="-", color=color, linewidth=0.75),
            meanprops=dict(zorder=101, mec="1.0", markerfacecolor=color,
                           linewidth=1),
            medianprops=dict(linewidth=0), showmeans=True, widths=5,
            showfliers=False)
        return boxplot_props

    def get_boxfillprops(self, mission):
        return dict(ec="none", fc=self.mission_colors[mission],
                    zorder=100)

    def get_labelfontprops(self, mission):
        return dict(color=self.mission_colors[mission],
                    fontproperties=get_custom_font(fontsize=16, awi_font=True))

    @property
    def output_filename(self):
        # Simplify region name
        regstr = self.region_name.lower()
        regstr = regstr.replace(" ", "_")
        # Simplify annotation
        if self.annotation is None:
            annstr = ""
        else:
            annstr = "_"+self.annotation.lower()
            annstr = annstr.replace(" ", "_")
        output_filename = "l4_%s_%s%s.png" % (self.parameter_name, regstr,
                                              annstr)
        return output_filename

    @property
    def savefigprops(self):
        return dict(dpi=600)

    @property
    def yearfontprops(self):
        return dict(color="#4b4b4d",
                    fontproperties=get_custom_font(fontsize=12, awi_font=True))

    def get_monthfontprops(self, xlim):
        # Get font size
        years = np.floor((xlim[1]-xlim[0])/365.)
        if years <= 2:
            fs = 11
        if years >= 5:
            fs = 5
        else:
            fs = 8
        return dict(color="#4b4b4d",
                    fontproperties=get_custom_font(fontsize=fs, awi_font=True))

    @property
    def yearboxprops(self):
        return dict(ec="#4b4b4d", fc="0.94", lw=0.05)

    @property
    def parlabelfontprops(self):
        return dict(color="#4b4b4d",
                    fontproperties=get_custom_font(fontsize=12, awi_font=True))

    @property
    def mission_colors(self):
        return dict(envisat="#00ace5", cryosat2="#003e6e")

    @property
    def custom_plot_xlim(self):
        if self.custom_plot_range is None:
            xlim = None
        else:
            rng = self.custom_plot_range
            startdt = datetime(rng[0], rng[1], 1)
            stopdt = datetime(rng[2], rng[3], 1)
            xlim = ([date2num(startdt), date2num(stopdt)])
        return xlim


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
        self._set_plot_ranges()
        self._set_plot_style()
        self._export(output_folder)

    def _init_figure(self):
        self.figure = plt.figure(figsize=(16, 6))
        self.ax = plt.gca()
        self.ax.set_position([0.05, 0.15, 0.92, 0.80])
        self.ax.set_axis_bgcolor('0.98')

    def _create_plot(self):
        """ Boxplots and Mission labels """

        # Get mission info (for full mission names)
        pysiral_config = ConfigInfo()

        # Loop over missions
        for i, mission in enumerate(self.l4.missions):

            # First draw the boxplots
            boxplotprops = self.plotdef.get_boxplotprops(mission)
            data, dates = self.l4.get_boxplot_data(mission)
            bp = plt.boxplot(data, positions=dates, **boxplotprops)

            # Fill boxes (for better visibility)
            for box in bp["boxes"]:
                xdata, ydata = box.get_xdata(), box.get_ydata()
                xmin, xmax = np.amin(xdata), np.amax(xdata)
                ymin, ymax = np.amin(ydata), np.amax(ydata)
                boxfillprops = self.plotdef.get_boxfillprops(mission)
                rect = Rectangle([xmin, ymin], xmax-xmin, ymax-ymin,
                                 **boxfillprops)
                self.ax.add_patch(rect)

            # Add the mission label
            mission_info = pysiral_config.get_mission_info(mission)
            labelx = 0.05 + i*0.075
            plt.annotate(mission_info.long_name, (labelx, 0.03),
                         xycoords="figure fraction",
                         **self.plotdef.get_labelfontprops(mission))

    def _set_plot_ranges(self):

        xticks, xticklabels = self.l4.get_all_periods()
        self.ax.set_xticks(xticks)
        self.ax.set_xticklabels(xticklabels)

        # Set plot ranges
#        if self.plotdef.custom_plot_xlim is None:
#            plt.xlim(self.l4.get_data_range(paddays=30))
#        else:
#            plt.xlim(self.plotdef.custom_plot_xlim)

        plt.ylim(-0.5, 4)

    def _set_plot_style(self):
        """ Axis Style """

        # Add year rectangles for better visibility
        for year in self.l4.years:
            labelx = date2num(datetime(year, 7, 1))
            plt.annotate(str(year), (labelx, 3.8), xycoords="data",
                         ha="center", **self.plotdef.yearfontprops)
            if np.mod(year, 2) == 0:
                continue
            xmin = date2num(datetime(year, 1, 1))
            ymin, width, height = -0.5, 365, 4.5
            rect = Rectangle([xmin, ymin], width, height,
                             **self.plotdef.yearboxprops)
            self.ax.add_patch(rect)

        # Remove missor axis
        spines_to_remove = ["top", "right"]
        for spine in spines_to_remove:
            self.ax.spines[spine].set_visible(False)

        # x-axis (x-ticks are forced to remove shifts for multiple missions)
        self.ax.xaxis.set_tick_params(direction='out')
        self.ax.xaxis.set_ticks_position('bottom')
        self.ax.spines["bottom"].set_position(("axes", -0.01))
        cl = plt.getp(self.ax, 'xmajorticklabels')
        plt.setp(cl, **self.plotdef.get_monthfontprops(self.ax.get_xlim()))

        # y-axis
        self.ax.set_ylabel("Sea Ice Thickness (meter)",
                           **self.plotdef.parlabelfontprops)
        self.ax.yaxis.set_tick_params(direction='out')
        self.ax.yaxis.set_ticks_position('left')
        self.ax.spines["left"].set_position(("axes", -0.01*6./16.))
        cl = plt.getp(self.ax, 'ymajorticklabels')
        plt.setp(cl, **self.plotdef.parlabelfontprops)

    def _export(self, output_folder):
        """ png export """
        validate_directory(output_folder)
        filename = os.path.join(output_folder, self.plotdef.output_filename)
        plt.savefig(filename, **self.plotdef.savefigprops)
        plt.close(self.figure)


class Level4ROIMap(object):

    def __init__(self):
        self.is_active = False
        self.roi = None
        self.annotation = None

    def set_roi(self, roi):
        self.roi = roi

    def export(self, output_folder):

        # get style from pysiral definition
        style = GridMapPaperStyle()

        # switch off interactive plotting
        plt.ioff()
        figure = plt.figure(**style.figure.keyw)

        # Basemap settings
        m = Basemap(**self.basemap_args)
        m.drawmapboundary(**style.mapboundary.keyw)
        m.fillcontinents(zorder=120, **style.continents.keyw)

#        parallels = np.arange(0., 88., 2.)
#        m.drawparallels(parallels, labels=[0, 0, 0, 0], fontsize=20, latmax=88,
#                        linewidth=0.25, dashes=(None, None), color="#4b4b4b")
#        meridians = np.arange(0., 360., 10.)
#        m.drawmeridians(meridians, labels=[0, 0, 0, 0], fontsize=20, latmax=88,
#                        linewidth=0.25, dashes=(None, None), color="#4b4b4b")

        if style.coastlines.is_active:
            coastlines = get_landcoastlines(m, **style.coastlines.keyw)
            coastlines.set_zorder(120)
            plt.gca().add_collection(coastlines)

        lons, lats = self.roi.get_boundary_line()
        x, y = m(lons, lats)
        plt.plot(x, y, linewidth=4, color="#FF2400", zorder=400)
#        plt.plot(x, y, linewidth=4, color="#4b4b4b", zorder=400)
        # "#76FF7A"
#        x = config.basemap_args["height"] / 2.
#        plt.plot(np.array([-x, x, x, -x, -x])+size/2.,
#                 np.array([-x, -x, x, x, -x])+size/2.,
#                 linewidth=4, color="#4b4b4b", zorder=400)

        # save the file
#        mapfile = config.output_basename + "_overview_map.png"
        filename = os.path.join(output_folder, self.output_filename)
        plt.savefig(filename, dpi=150, facecolor=figure.get_facecolor(),
                    bbox_inches="tight")
        plt.close(figure)

    @property
    def basemap_args(self):
        return dict(projection="laea", resolution="i",  width=6000000,
                    height=6000000, lon_0=-45, lat_0=90)

    @property
    def output_filename(self):
        # Simplify region name
        regstr = self.region_name.lower()
        regstr = regstr.replace(" ", "_")
        # Simplify annotation
        if self.annotation is None:
            annstr = ""
        else:
            annstr = "_"+self.annotation.lower()
            annstr = annstr.replace(" ", "_")
        output_filename = "l4map_%s%s.png" % (regstr, annstr)
        return output_filename


if __name__ == "__main__":
    job = L3ParameterTimeSeries()
    job.init_from_args()
    job.execute()
