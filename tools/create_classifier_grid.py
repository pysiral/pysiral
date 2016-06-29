# -*- coding: utf-8 -*-

from pysiral.logging import DefaultLoggingClass
from pysiral.config import ConfigInfo
from pysiral.iotools import get_l1bdata_files
from pysiral.l1bdata import L1bdataNCFile
from pysiral.path import validate_directory

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

        return parser


class ClassifierGrid(DefaultLoggingClass):

    def __init__(self, jobdef):
        super(ClassifierGrid, self).__init__(self.__class__.__name__)
        self.jobdef = jobdef
        self.classifier_list = []
        self.parameter_stack = {}
        self.parameter_grid = {}

    def run(self):
        self.init_parameter_stacks()
        self.init_parameter_grids()
        self.init_projection()
        self.init_longitude_latitude_grids()
        self.create_orbit_stacks()
        self.grid_orbit_stacks()
        self.export_stacks()

        import matplotlib.pyplot as plt

        for parameter_name in self.classifier_list:
            data = self.parameter_grid[parameter_name]
            plt.figure(parameter_name, figsize=(12, 12),
                       facecolor="white")
            plt.imshow(np.flipud(data), interpolation="none",
                       cmap=plt.get_cmap("plasma"))
            plt.colorbar()
        plt.show()

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

    def create_orbit_stacks(self):
        """ Loop over l1b data files """
        n_l1b_files = len(self.jobdef.l1b_files)
        for i, l1b_file in enumerate(self.jobdef.l1b_files):
            self.log.info("+ Parsing file %g of %g: %s" % (
                i+1, n_l1b_files, l1b_file))
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
        for xi in self.grid_xi_range:
            for yj in self.grid_yj_range:
                surface_type = np.array(self.surface_type_stack[yj][xi])
                is_land = len(np.where(surface_type == 7)[0]) > 0
                if is_land:
                    continue
                for name in self.classifier_list:
                    data = np.array(self.parameter_stack[name][yj][xi])
                    if len(np.where(np.isfinite(data))[0]) > 1:
                        self.parameter_grid[name][yj, xi] = np.nanmean(data)

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

    def _get_stack_output_filename(self, parameter_name):
        mission_id = self.jobdef.args.mission_id
        grid_tag = self.jobdef.griddef.grid_tag.lower()
        hemisphere_tag = self.jobdef.griddef.hemisphere
        resolution_tag = self.jobdef.griddef.resolution_tag
        year = self.jobdef.args.month[0]
        month = self.jobdef.args.month[1]
        filename = "%s_classifier_%04g%02g_%s_%s_%s_%s.npz" % (
            mission_id, year, month, parameter_name, grid_tag, hemisphere_tag,
            resolution_tag)
        return filename

    @property
    def grid_xi_range(self):
        return np.arange(self.jobdef.griddef.extent.numx)

    @property
    def grid_yj_range(self):
        return np.arange(self.jobdef.griddef.extent.numy)

if __name__ == "__main__":
    create_classifier_grids()
