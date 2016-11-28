# -*- coding: utf-8 -*-

from pysiral.logging import DefaultLoggingClass
from pysiral.config import ConfigInfo, TimeRangeRequest
from pysiral.iotools import get_l1bdata_files
from pysiral.l1bdata import L1bdataNCFile
from pysiral.path import validate_directory
from pysiral.l3proc import L3MetaData
from pysiral.output import L3SDataNC

from scipy.stats import skew, kurtosis
from pyproj import Proj
import numpy as np
import sys
import os


def create_classifier_grids():

    # get the pysiral configuration info
    config = ConfigInfo()

    # parse command line arguments
    jobdef = ClassifierGridSettings(config)
    jobdef.parse_command_line_arguments()
    jobdef.collect_job_parameter()

    # start the job
    job = ClassifierGrid(jobdef)
    job.run()


class ClassifierGridSettings(DefaultLoggingClass):

    def __init__(self, config):
        super(ClassifierGridSettings, self).__init__(self.__class__.__name__)
        self.pysiral_config = config
        self.args = None
        self.griddef = None
        self.output_folder = None
        self.l1b_files = []

    def parse_command_line_arguments(self):
        self.args = self.parser.parse_args()
        self._validate_args()

    def collect_job_parameter(self):
        # get list of l1b data files
        self.get_l1b_file_list()
        # get grid definition
        self.griddef = self.pysiral_config.area[self.args.grid_id]
        # Create output folder
        self.create_output_folder()

    def create_output_folder(self):
        path = self.args.output_directory
        month = self.args.month
        path = os.path.join(path, "%04g" % month[0], "%02g" % month[1])
        validate_directory(path)
        self.log.info("Export directory: %s" % path)
        self.output_folder = path

    def get_l1b_file_list(self):
        year, month = int(self.args.month[0]), int(self.args.month[1])
        self.l1b_files = get_l1bdata_files(
            self.args.mission_id, self.args.hemisphere,
            year, month, config=self.pysiral_config)
        self.log.info("Found %g l1bdata files" % len(self.l1b_files))

    def _validate_args(self):

        # Check if Valid Mission
        args = self.args
        if args.mission_id not in self.pysiral_config.mission_ids:
            error_message = "Error: Invalid mission id (%s)" % args.mission_id
            error_message += "\nuse: " + " | ".join(
                self.pysiral_config.mission_ids)
            self.log.error(error_message)
            sys.exit(1)
        else:
            info_message = "valid mission id \"%s\"" % self.args.mission_id
            self.log.info(info_message)

        # Simple Check if month data is ok
        try:
            int(args.month[0])
            int(args.month[1])
            info_message = "valid month \"%s\"" % str(self.args.month)
            self.log.info(info_message)
        except:
            error_message = "Wrong month definition (yyyy mm): %s"
            error_message = error_message % str(args.month)
            self.log.error(error_message)
            sys.exit(1)

        # check if region id exists
        from pysiral.config import td_branches
        pysiral_grid_ids, objects = td_branches(self.pysiral_config.area)
        if self.args.grid_id not in pysiral_grid_ids:
            error_message = "grid id \"%s\" not recognized" % self.args.grid_id
            self.log.error(error_message)
            sys.exit(1)
        else:
            info_message = "valid grid id \"%s\"" % self.args.grid_id
            self.log.info(info_message)

        # Test if output path is ok
        path = self.args.output_directory
        if not os.access(os.path.dirname(path), os.W_OK):
            error_message = "Invalid output folder: %s" % path
            self.log.error(error_message)
            sys.exit(1)
        else:
            info_message = "valid output folder: %s" % path
            self.log.info(info_message)

    @property
    def parser(self):

        import argparse

        parser = argparse.ArgumentParser()

        parser.add_argument(
            '-mission', action='store', dest='mission_id',
            help='mission_id', required=True)

        parser.add_argument(
            '-output', action='store', dest='output_directory',
            help='Output Directory')

        parser.add_argument(
            '-grid', action='store', dest='grid_id',
            help='grid id (see config/area_def.yaml')

        parser.add_argument(
            '-hemisphere', action='store', dest='hemisphere',
            help='hemisphere (north or south)', type=str, default="north")

        parser.add_argument(
            '-month', action='store', dest='month',
            nargs='+', type=int, required=True,
            help='year and month (-start yyyy mm)')

        parser.add_argument(
            '--export-stack', action='store_true', dest='export_stack',
            default=False, help='export stack as *.npz')

        return parser

    @property
    def metadata_dict(self):
        period = TimeRangeRequest()
        start_date = [int(self.args.month[0]), int(self.args.month[1])]
        stop_date = start_date
        period.set_range(start_date, stop_date)
        metadata_dict = {}
        metadata_dict["start_time"] = period.start_dt
        metadata_dict["stop_time"] = period.stop_dt
        metadata_dict["period_label"] = period.start_dt.strftime("%B %Y")
        metadata_dict["grid_tag"] = self.griddef.grid_tag
        metadata_dict["hemisphere"] = self.griddef.hemisphere
        metadata_dict["resolution_tag"] = self.griddef.resolution_tag
        metadata_dict["mission_ids"] = self.args.mission_id
        return metadata_dict

    @property
    def nc_dimdict(self):
        from collections import OrderedDict
        dimdict = OrderedDict([("numx", self.griddef.extent.numx),
                               ("numy", self.griddef.extent.numy)])
        return dimdict


class ClassifierGrid(DefaultLoggingClass):

    def __init__(self, jobdef):
        super(ClassifierGrid, self).__init__(self.__class__.__name__)
        self.jobdef = jobdef
        self.classifier_list = []
        self.parameter_stack = {}
        self.parameter_grid = {}
        self._progress_percent = 0

    def run(self):
        self.init_parameter_stacks()
        self.init_parameter_grids()
        self.init_projection()
        self.init_longitude_latitude_grids()
        self.create_orbit_stacks()
        self.grid_orbit_stacks()
        if self.jobdef.args.export_stack:
            self.export_stacks()
        self.export_l3_grid()

#        import matplotlib.pyplot as plt
#
#        for parameter_name in self.classifier_list:
#            data = self.parameter_grid[parameter_name]
#            plt.figure(parameter_name, figsize=(12, 12),
#                       facecolor="white")
#            plt.imshow(np.flipud(data), interpolation="none",
#                       cmap=plt.get_cmap("plasma"))
#            plt.colorbar()
#        plt.show()

    def init_parameter_stacks(self):
        # get the list of clasiifiers from the first file in the list
        l1b_file = self.jobdef.l1b_files[0]
        l1b = L1bdataNCFile(l1b_file)
        l1b.parse()
        self.classifier_list = l1b.classifier.parameter_list
        # Initialize the grid for the classifier stacke
        griddef = self.jobdef.griddef
        dimx, dimy = griddef.extent.numx, griddef.extent.numy
        self.surface_type_stack = \
            [[[] for _ in range(dimx)] for _ in range(dimy)]
        for parameter_name in self.classifier_list:
            self.parameter_stack[parameter_name] = \
                [[[] for _ in range(dimx)] for _ in range(dimy)]

    def init_parameter_grids(self):
        # Initialize the grid for the classifier stacke
        griddef = self.jobdef.griddef
        shape = (griddef.extent.numx, griddef.extent.numy)
        for parameter_name in self.classifier_list:
            self.parameter_grid[parameter_name] = \
                np.ndarray(shape=shape, dtype='f4')*np.nan
            for statmoment in self.statmoments:
                statmoment_name = parameter_name+"_"+statmoment
                self.parameter_grid[statmoment_name] = \
                    np.ndarray(shape=shape, dtype='f4')*np.nan

    def init_projection(self):
        self.proj = Proj(**self.jobdef.griddef.projection)

    def init_longitude_latitude_grids(self):

        griddef = self.jobdef.griddef
        x0 = griddef.extent.xoff
        xsize = griddef.extent.xsize
        numx = griddef.extent.numx
        xmin, xmax = x0-(xsize/2.), x0+(xsize/2.)

        y0 = griddef.extent.yoff
        ysize = griddef.extent.ysize
        numy = griddef.extent.numy
        ymin, ymax = y0-ysize/2., y0+ysize/2.

        x = np.linspace(xmin, xmax, num=numx)
        y = np.linspace(ymin, ymax, num=numy)
        xx, yy = np.meshgrid(x, y)
        lon, lat = self.proj(xx, yy, inverse=True)

        self.longitude = lon
        self.latitude = lat
        self.parameter_grid["longitude"] = lon
        self.parameter_grid["latitude"] = lat

    def create_orbit_stacks(self):
        """ Loop over l1b data files """
        n_l1b_files = len(self.jobdef.l1b_files)
        for i, l1b_file in enumerate(self.jobdef.l1b_files):
            self._log_progress(i, n_l1b_files)
            l1b = L1bdataNCFile(l1b_file)
            l1b.parse()
            self.append_to_stack(l1b)

    def append_to_stack(self, l1b):
        """ Append classifier data of single l1b object to stack """
        xi, yj = self.get_l1b_track_grid_indices(l1b)
        for i in np.arange(l1b.n_records):
            x, y = int(xi[i]), int(yj[i])
            self.surface_type_stack[y][x].append(l1b.surface_type.flag[i])
            for parameter_name in self.classifier_list:
                data = l1b.classifier.get_parameter(parameter_name)
                self.parameter_stack[parameter_name][y][x].append(data[i])

    def grid_orbit_stacks(self):
        for name in self.classifier_list:
            self.log.info("Gridding classifier: %s" % name)
            for xi in self.grid_xi_range:
                for yj in self.grid_yj_range:
                    surface_type = np.array(self.surface_type_stack[yj][xi])
                    is_land = len(np.where(surface_type == 7)[0]) > 0
                    if is_land:
                        continue
                    self.create_gridcell_statistics(name, xi, yj)

    def create_gridcell_statistics(self, name, xi, yj):
        scipy_nanprops = {"nan_policy": "omit"}
        grid = self.parameter_grid
        data = np.array(self.parameter_stack[name][yj][xi])
        if len(np.where(np.isfinite(data))[0]) > 1:
            grid[name][yj, xi] = np.nanmean(data)
            namevar = name+"_variance"
            grid[namevar][yj, xi] = np.nanvar(data)
            nameskew = name+"_skewness"
            grid[nameskew][yj, xi] = skew(data, **scipy_nanprops)
            namekurt = name+"_kurtosis"
            grid[namekurt][yj, xi] = kurtosis(data, **scipy_nanprops)

    def get_l1b_track_grid_indices(self, l1b):
        projx, projy = self.proj(l1b.time_orbit.longitude,
                                 l1b.time_orbit.latitude)
        # Convert projection coordinates to grid indices
        extent = self.jobdef.griddef.extent
        xi = np.floor((projx + extent.xsize/2.0)/extent.dx)
        yj = np.floor((projy + extent.ysize/2.0)/extent.dy)
        return xi, yj

    def export_stacks(self):
        for parameter_name in self.classifier_list:
            filename = self._get_stack_output_filename(parameter_name)
            path = os.path.join(self.jobdef.output_folder, filename)
            self.log.info("Exporting %s stack to file %s" % (
                parameter_name, path))
            array = np.array(self.parameter_stack[parameter_name])
            np.savez_compressed(path, array)

    def export_l3_grid(self):
        l3_metadata = L3MetaData()
        metadata_dict = self.jobdef.metadata_dict
        for key in metadata_dict.keys():
            l3_metadata.set_attribute(key, metadata_dict[key])

        # Write grid to L3S netcdf file:
        self.log.info("Exporting file")
        output = L3SDataNC()
        output.set_metadata(l3_metadata)
        output.set_export_folder(self.jobdef.output_folder)
        output.export_parameter_dict(self.parameter_grid,
                                     dimdict=self.jobdef.nc_dimdict)

    def _get_stack_output_filename(self, parameter_name):
        mission_id = self.jobdef.args.mission_id
        grid_tag = self.jobdef.griddef.grid_tag.lower()
        hemisphere_tag = self.jobdef.griddef.hemisphere
        resolution_tag = self.jobdef.griddef.resolution_tag
        year = self.jobdef.args.month[0]
        month = self.jobdef.args.month[1]
        filename = "%s_classifier_%04g%02g_%s_%s_%s_%s.npz" % (
            mission_id, year, month, grid_tag, hemisphere_tag,
            resolution_tag, parameter_name)
        return filename

    def _log_progress(self, i, n):
        progress_percent = float(i)/float(n-1)*100.
        current_reminder = np.mod(progress_percent, 10)
        last_reminder = np.mod(self._progress_percent, 10)
        if last_reminder > current_reminder:
            self.log.info("Parsing l1b orbit stack: %3g%% (%g of %g)" % (
                          progress_percent-current_reminder, i+1, n))
        self._progress_percent = progress_percent

    @property
    def grid_xi_range(self):
        return np.arange(self.jobdef.griddef.extent.numx)

    @property
    def grid_yj_range(self):
        return np.arange(self.jobdef.griddef.extent.numy)

    @property
    def statmoments(self):
        return ["variance", "skewness", "kurtosis"]

if __name__ == "__main__":
    create_classifier_grids()
