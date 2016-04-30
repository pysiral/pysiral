# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 14:04:27 2015

@author: Stefan
"""

from pysiral.logging import DefaultLoggingClass
from pysiral.l2data import L2iNCFile
import numpy as np


class Level3Processor(DefaultLoggingClass):

    def __init__(self, job):
        super(Level3Processor, self).__init__("Level3Processor")
        self.job = job
        self._l2_files = []

    def set_l2_files(self, l2_files):
        self._l2_files = l2_files

    def run(self):

        # Initialize the stack for the l2i orbit files
        stack = L2DataStack()
        stack.set_grid_definition(self.job.grid)
        stack.set_l2_parameter(self.job.l2_parameter)
        stack.initialize()

        # Parse all orbit files and add to the stack
        for i, l2_file in enumerate(self._l2_files):
            self.log.info("+ Parsing file %g of %g: %s" % (
                i, len(self._l2_files), l2_file))
            l2data = L2iNCFile(l2_file)
            l2data.project(self.job.grid)
            stack.append(l2data)

        # Initialize the data grid
        grid = L3DataGrid(self.job.grid)
        grid.init_l2_parameter(self.job.l2_parameter)
        grid.init_l3_parameter(self.job.l3_parameter)

        # Average level-2 parameter for each grid cell
        grid.set_l2_stack(stack)
        grid.average_l2_parameter()

        # Compute level-2 parameter (lead_fraction ....)


class L2DataStack(object):

    def __init__(self):
        self.griddef = None
        self.l2_parameter = None
        self.n_records = 0
        self.stack = {}

    def set_grid_definition(self, griddef):
        self.griddef = griddef

    def set_l2_parameter(self, l2_parameter):
        self.l2_parameter = l2_parameter

    def initialize(self):
        dimx, dimy = self.griddef.extent.numx, self.griddef.extent.numy
        for parameter_name in self.l2_parameter:
            self.stack[parameter_name] = \
                [[[] for _ in range(dimx)] for _ in range(dimy)]

    def append(self, orbit):
        self.n_records += orbit.n_records
        for i in np.arange(orbit.n_records):
            x, y = int(orbit.xi[i]), int(orbit.yj[i])
            for parameter_name in self.l2_parameter:
                data = getattr(orbit, parameter_name)
                self.stack[parameter_name][y][x].append(data[i])


class L3DataGrid(object):
    """ Contains gridded orbit data parameters plus derived parameter  """

    def __init__(self, griddef):
        self.griddef = griddef
        self._l2_parameter = None
        self._l2 = None
        self._l3 = {}

    def init_l2_parameter(self, l2_parameter):
        shape = (self.griddef.extent.numx, self.griddef.extent.numy)
        # Taken from gridded SICCI data
        self.longitude = np.ndarray(shape=shape, dtype='f4')*np.nan
        self.latitude = np.ndarray(shape=shape, dtype='f4')*np.nan
        for parameter_name in l2_parameter:
            self._l3[parameter_name] = np.ndarray(
                shape=shape, dtype='f4')*np.nan

    def init_l3_parameter(self, l3_parameter):
        shape = (self.griddef.extent.numx, self.griddef.extent.numy)
        # Taken from gridded SICCI data
        pass

    def set_l2_stack(self, stack):
        self._l2 = stack

    def average_l2_parameter(self):
        n_grid_cells = self.griddef.extent.numx * self.griddef.extent.numy

        # Loop over all grid cells
        # XXX: Is there a better way?
        for xi in self.grid_xi_range:
            for yj in self.grid_yj_range:
                sit = np.array(self._l2.stack["sea_ice_thickness"][yj][xi])
                self._l3["sea_ice_thickness"][yj, xi] = np.nanmean(sit)
        import matplotlib.pyplot as plt
        plt.figure()
        plt.imshow(
            np.flipud(self._l3["sea_ice_thickness"]), vmin=0, vmax=5,
            cmap=plt.get_cmap("plasma"), interpolation="none")
        plt.colorbar()
        plt.show()
        stop

    def flipud(self):
        for parameter in self.parameters:
            setattr(self, parameter, np.flipud(getattr(self, parameter)))

    def sync_nans(self):
        """
        Transfer NaN mask from freeboard fields to parameters derived
        from orbit data. (There is additional filtering in the freeboard
        grids and there should not be a lead fraction where there is no
        freeboard)
        """
        nan_list = np.where(np.isnan(self.radar_freeboard))
        except_list = ['n_valid_waveforms', 'lead_fraction', 'ocean_fraction',
                       'ice_fraction', 'radar_freeboard']
        for parameter in self.parameters:
            if parameter in except_list:
                continue
            data = getattr(self, parameter)
            data[nan_list] = np.nan
            setattr(self, parameter, data)

    @property
    def grid_xi_range(self):
        return np.arange(self.griddef.extent.numx)

    @property
    def grid_yj_range(self):
        return np.arange(self.griddef.extent.numy)