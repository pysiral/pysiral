# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 14:04:27 2015

@author: Stefan
"""

from pysiral.logging import DefaultLoggingClass
from pysiral.l2data import L2iNCFileImport
from pysiral.output import L3SDataNC
import numpy as np


# %% Level 3 Processor

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
        stack.initialize(len(self._l2_files))

        # Initialize the data grid
        grid = L3DataGrid(self.job.grid)
        grid.init_l2_parameter(self.job.l2_parameter)
        grid.init_l3_parameter(self.job.l3_parameter)
        grid.calculate_longitude_latitude_fields()

        # Parse all orbit files and add to the stack
        for i, l2_file in enumerate(self._l2_files):
            self.log.info("+ Parsing file %g of %g: %s" % (
                i+1, len(self._l2_files), l2_file))
            l2data = L2iNCFileImport(l2_file)
            l2data.project(self.job.grid)
            stack.append(l2data)

        # Average level-2 parameter for each grid cell
        grid.set_l2_stack(stack)
        grid.average_l2_parameter()
        grid.set_freeboard_nan_mask(self.job.l2_freeboard_nan_mask_targets)

        # TODO: Compute level-3 parameter (lead_fraction ....)

        # Get the metadata information from the L2 stack
        l3_metadata = L3MetaData()
        l3_metadata.get_missions_from_stack(stack)
        l3_metadata.get_data_period_from_stack(stack)
        l3_metadata.get_projection_parameter(self.job.grid)

        # Write grid to L3S netcdf file:
        output = L3SDataNC()
        output.set_metadata(l3_metadata)
        output.set_export_folder(self.job.export_folder)
        output.export(grid)


# %% Helper Classes

class L2DataStack(object):

    def __init__(self):
        self.griddef = None
        self.l2_parameter = None
        self.n_records = 0
        self.stack = {}
        self.l2_count = 0

    def set_grid_definition(self, griddef):
        self.griddef = griddef

    def set_l2_parameter(self, l2_parameter):
        self.l2_parameter = l2_parameter

    def initialize(self, n_files):
        # attributes for grid metadata
        self.start_time = np.ndarray(shape=(n_files), dtype=object)
        self.stop_time = np.ndarray(shape=(n_files), dtype=object)
        self.mission = np.ndarray(shape=(n_files), dtype=object)
        # gridded parameter
        dimx, dimy = self.griddef.extent.numx, self.griddef.extent.numy
        for parameter_name in self.l2_parameter:
            self.stack[parameter_name] = \
                [[[] for _ in range(dimx)] for _ in range(dimy)]

    def append(self, orbit):
        # Save the metadata from the orbit data
        self.start_time[self.l2_count] = orbit.timestamp[0]
        self.stop_time[self.l2_count] = orbit.timestamp[-1]
        self.mission[self.l2_count] = orbit.mission
        self.l2_count += 1
        # Stack the l2 parameter in the corresponding grid cells
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
        self._l2_parameter = l2_parameter
        shape = (self.griddef.extent.numx, self.griddef.extent.numy)
        # Taken from gridded SICCI data
        self.longitude = np.ndarray(shape=shape, dtype='f4')*np.nan
        self.latitude = np.ndarray(shape=shape, dtype='f4')*np.nan
        for parameter_name in l2_parameter:
            self._l3[parameter_name] = np.ndarray(
                shape=shape, dtype='f4')*np.nan

    def init_l3_parameter(self, l3_parameter):
        # shape = (self.griddef.extent.numx, self.griddef.extent.numy)
        # derived parameters from the stack
        pass

    def calculate_longitude_latitude_fields(self):

        # TODO: Move all this into a grid container
        from pyproj import Proj
        proj = Proj(**self.griddef.projection)

        x0 = self.griddef.extent.xoff
        xsize = self.griddef.extent.xsize
        numx = self.griddef.extent.numx
        xmin, xmax = x0-(xsize/2.), x0+(xsize/2.)

        y0 = self.griddef.extent.yoff
        ysize = self.griddef.extent.ysize
        numy = self.griddef.extent.numy
        ymin, ymax = y0-ysize/2., y0+ysize/2.

        x = np.linspace(xmin, xmax, num=numx)
        y = np.linspace(ymin, ymax, num=numy)
        xx, yy = np.meshgrid(x, y)
        lon, lat = proj(xx, yy, inverse=True)

        self._l3["longitude"] = lon
        self._l3["latitude"] = lat

#        import matplotlib.pyplot as plt
#
#        plt.figure("xx")
#        plt.imshow(xx, cmap=plt.get_cmap("rainbow"))
#        plt.colorbar()
#
#        plt.figure("yy")
#        plt.imshow(yy, cmap=plt.get_cmap("rainbow"))
#        plt.colorbar()
#
#        plt.figure("longitude")
#        plt.imshow(lon, cmap=plt.get_cmap("rainbow"))
#        plt.colorbar()
#
#        plt.figure("lat")
#        plt.imshow(lat, cmap=plt.get_cmap("rainbow"))
#        plt.colorbar()
#
#        plt.show()
#        stop

    def set_l2_stack(self, stack):
        self._l2 = stack

    def average_l2_parameter(self):
        # Loop over all grid cells
        # XXX: Is there a better way?
        for xi in self.grid_xi_range:
            for yj in self.grid_yj_range:
                for name in self._l2_parameter:
                    data = np.array(self._l2.stack[name][yj][xi])
                    self._l3[name][yj, xi] = np.nanmean(data)

    def set_freeboard_nan_mask(self, targets):
        """
        Apply the freeboard nan mask to a selected number of parameters
        see setting/l3/*.yaml for details
        """
        frb_is_nan = np.where(np.isnan(self._l3["freeboard"]))
        for target in targets:
            self._l3[target][frb_is_nan] = np.nan

    def get_parameter_by_name(self, name):
        return self._l3[name]


#        import matplotlib.pyplot as plt
#
#        plt.figure("sit")
#        plt.imshow(
#            np.flipud(self._l3["sea_ice_thickness"]), vmin=0, vmax=5,
#            cmap=plt.get_cmap("plasma"), interpolation="none")
#        plt.colorbar()
#
#        plt.figure("frb")
#        plt.imshow(
#            np.flipud(self._l3["freeboard"]), vmin=0, vmax=0.5,
#            cmap=plt.get_cmap("viridis"), interpolation="none")
#        plt.colorbar()
#
#        plt.figure("snow")
#        plt.imshow(
#            np.flipud(self._l3["snow_depth"]), vmin=0, vmax=0.6,
#            cmap=plt.get_cmap("magma"), interpolation="none")
#        plt.colorbar()
#
#        plt.figure("ice_density")
#        plt.imshow(
#            np.flipud(self._l3["ice_density"]), vmin=800, vmax=1000,
#            cmap=plt.get_cmap("magma"), interpolation="none")
#        plt.colorbar()
#
#        plt.figure("sea_surface_anomaly")
#        plt.imshow(
#            np.flipud(self._l3["sea_surface_anomaly"]), vmin=-1, vmax=1,
#            cmap=plt.get_cmap("bwr"), interpolation="none")
#        plt.colorbar()
#
#        plt.figure("mean_sea_surface")
#        plt.imshow(
#            np.flipud(self._l3["mean_sea_surface"]), vmin=-10, vmax=30,
#            cmap=plt.get_cmap("gist_earth"), interpolation="none")
#        plt.colorbar()
#
#        plt.figure("sea_ice_type")
#        plt.imshow(
#            np.flipud(self._l3["sea_ice_type"]), vmin=0, vmax=1,
#            cmap=plt.get_cmap("inferno"), interpolation="none")
#        plt.colorbar()
#
#        plt.show()
#        stop

    def flipud(self):
        for parameter in self.parameters:
            setattr(self, parameter, np.flipud(getattr(self, parameter)))

    @property
    def grid_xi_range(self):
        return np.arange(self.griddef.extent.numx)

    @property
    def grid_yj_range(self):
        return np.arange(self.griddef.extent.numy)

    @property
    def parameter_list(self):
        # TODO: Only L2 parameter for now
        parameter_list = self._l2_parameter
        parameter_list.append("longitude")
        parameter_list.append("latitude")
        return parameter_list

    @property
    def dimdict(self):
        from collections import OrderedDict
        dimdict = OrderedDict([("numx", self.griddef.extent.numx),
                               ("numy", self.griddef.extent.numy)])
        return dimdict


class L3MetaData(object):

    """
    Container for L3S Metadata information
    (see property attribute_list for a list of attributes)
    """

    _attribute_list = [
        "mission_ids", "start_time", "stop_time", "grid_name", "period_label",
        "pysiral_version", "projection_str"]

    def __init__(self):
        # Init all fields
        for field in self.attribute_list:
            setattr(self, field, None)

    def get_missions_from_stack(self, stack):
        """
        Get a list of missions that went into the stack
        (must be a list, since multi-mission grids are supported)
        """
        self.set_attribute("mission_ids", ",".join(np.unique(stack.mission)))

    def get_data_period_from_stack(self, stack):
        """ Get the first and last timestamp """
        self.set_attribute("start_time", np.amin(stack.start_time))
        self.set_attribute("stop_time", np.amax(stack.stop_time))
        # XXX: Only monthly periods are currently supported
        self.set_attribute("period_label", self.start_time.strftime("%B %Y"))

    def get_projection_parameter(self, griddef):
        pass

    def __repr__(self):
        output = "pysiral.L3S Metadata Container:\n"
        for field in self._attribute_list:
            output += "%22s: %s" % (field, getattr(self, field))
            output += "\n"
        return output

    @property
    def attribute_list(self):
        return self._attribute_list

    @property
    def attdict(self):
        """ Return attributes as dictionary (e.g. for netCDF export) """
        attdict = {}
        for field in self.attribute_list:
            attdict[field] = getattr(self, field)
        return attdict

    @property
    def mission(self):
        mission_ids = self.mission_ids.split(",")
        return "_".join(mission_ids)

    @property
    def start_period(self):
        return self.start_time.strftime("%Y%m%d")

    @property
    def stop_period(self):
        return self.stop_time.strftime("%Y%m%d")

    def set_attribute(self, tag, value):
        if tag not in self.attribute_list:
            raise ValueError("Unknown attribute: ", tag)
        setattr(self, tag, value)
