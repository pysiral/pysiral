# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 14:04:27 2015

@author: Stefan
"""

from pysiral.logging import DefaultLoggingClass
from pysiral.l2data import L2iNCFileImport
from pysiral.output import L3SDataNC
from pysiral.flag import ORCondition
from pysiral.surface_type import SurfaceType
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
        grid.init_parameter_fields(self.job.l2_parameter, "l2")
        grid.init_parameter_fields(self.job.l3_parameter, "l3")
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


        # Get the level-3 parameter
        grid.compute_l3_parameter()

        # Set parameters nan if freeboard is nan
        # (list in output definition file)
        grid.set_freeboard_nan_mask(self.job.freeboard_nan_mask_targets)

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

class L2DataStack(DefaultLoggingClass):

    def __init__(self):
        super(L2DataStack, self).__init__(self.__class__.__name__)
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
        self.stack["surface_type"] = \
            [[[] for _ in range(dimx)] for _ in range(dimy)]
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
#        import matplotlib.pyplot as plt
#        plt.figure()
#        plt.plot(orbit.surface_type)
#        plt.show()
        for i in np.arange(orbit.n_records):
            # Add the surface type per default
            # (will not be gridded, therefore not in list of l2 parameter)
            x, y = int(orbit.xi[i]), int(orbit.yj[i])
            self.stack["surface_type"][y][x].append(orbit.surface_type[i])
            for parameter_name in self.l2_parameter:
                data = getattr(orbit, parameter_name)
                self.stack[parameter_name][y][x].append(data[i])


class L3DataGrid(DefaultLoggingClass):
    """ Contains gridded orbit data parameters plus derived parameter  """

    def __init__(self, griddef):
        super(L3DataGrid, self).__init__(self.__class__.__name__)
        self.griddef = griddef
        self._surface_type_dict = SurfaceType.SURFACE_TYPE_DICT
        self._l2_parameter = None
        self._l3_parameter = None
        self._l2 = None
        self._l3 = {}

    def init_parameter_fields(self, parameter_names, level):
        setattr(self, "_"+level+"_parameter", parameter_names)
        shape = (self.griddef.extent.numx, self.griddef.extent.numy)
        # Taken from gridded SICCI data
        # XXX: There needs to be a better handling of data types
        #       (requires better definition of l3 output file format)
        self.longitude = np.ndarray(shape=shape, dtype='f4')*np.nan
        self.latitude = np.ndarray(shape=shape, dtype='f4')*np.nan
        self.is_land = np.ndarray(shape=shape, dtype='f4')*np.nan
        for parameter_name in parameter_names:
            msg = "Adding %s parameter: %s" % (level, parameter_name)
            self.log.info(msg)
            self._l3[parameter_name] = np.ndarray(
                shape=shape, dtype='f4')*np.nan

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

                # Filter land values
                surface_type = np.array(self._l2.stack["surface_type"][yj][xi])
                # TODO: this needs to be done more formalized
                # surface_type (land = 7)
                self.is_land[xi, yj] = len(np.where(surface_type == 7)[0] > 0)
                if self.is_land[xi, yj]:
                    continue
                for name in self._l2_parameter:
                    data = np.array(self._l2.stack[name][yj][xi])
                    if len(np.where(np.isfinite(data))[0]) > 1:
                        self._l3[name][yj, xi] = np.nanmean(data)

    def set_freeboard_nan_mask(self, targets):
        """
        Apply the freeboard nan mask to a selected number of parameters
        see setting/l3/*.yaml for details
        """
        frb_is_nan = np.where(np.isnan(self._l3["freeboard"]))
        for target in targets:
            self._l3[target][frb_is_nan] = np.nan

    def compute_l3_parameter(self):

        # Loop over grid items
        for xi in self.grid_xi_range:
            for yj in self.grid_yj_range:
                for l3_parameter_name in self._l3_parameter:
                    result = self.get_l3_parameter(l3_parameter_name, xi, yj)
                    self._l3[l3_parameter_name][yj, xi] = result

    def get_l3_parameter(self, l3_parameter_name, xi, yj):
        """
        XXX: Currently only works for a defined set of l3 parameter
             Better implementation needed for future flexibility
        """
        surface_type = np.array(self._l2.stack["surface_type"][yj][xi])
        stflags = self._surface_type_dict

        # Number of all Waveforms (including unknown, invalid etc)
        if l3_parameter_name == "n_total_waveforms":
            return len(surface_type)

        # Only positively identified waveforms (either lead or ice)
        # XXX: what about polynay and ocean?
        elif l3_parameter_name == "n_valid_waveforms":
            valid_waveform = ORCondition()
            valid_waveform.add(surface_type == stflags["lead"])
            valid_waveform.add(surface_type == stflags["sea_ice"])
            return valid_waveform.num

        # Fractions of leads on valid_waveforms
        elif l3_parameter_name == "valid_fraction":
            n_valid = self.get_l3_parameter("n_valid_waveforms", xi, yj)
            n_total = self.get_l3_parameter("n_total_waveforms", xi, yj)
            try:
                valid_fraction = float(n_valid)/float(n_total)
            except:
                valid_fraction = 0.0
            return valid_fraction

        # Fractions of leads on valid_waveforms
        elif l3_parameter_name == "lead_fraction":
            n_leads = len(np.where(surface_type == stflags["lead"]))
            n_valid = self.get_l3_parameter("n_valid_waveforms", xi, yj)
            try:
                lead_fraction = float(n_leads)/float(n_valid)
            except:
                lead_fraction = 0.0
            return lead_fraction

        # Fractions of leads on valid_waveforms
        elif l3_parameter_name == "ice_fraction":
            n_valid = self.get_l3_parameter("n_valid_waveforms", xi, yj)
            n_ice = np.where(surface_type == stflags["sea_ice"])
            try:
                ice_fraction = float(n_ice)/float(n_valid)
            except:
                ice_fraction = 0.0
            return ice_fraction

        # Raise if unknown l3 parameter
        else:
            msg = "Unknown L3 parameter name: %s" % l3_parameter_name
            raise ValueError(msg)

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
        parameter_list.extend(self._l3_parameter)
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
        "pysiral_version", "projection_str", "grid_tag", "resolution_tag"]

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
        self.set_attribute("grid_tag", griddef.grid_tag)
        self.set_attribute("resolution_tag", griddef.resolution_tag)

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
