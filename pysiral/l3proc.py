# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 14:04:27 2015

@author: Stefan
"""

from pysiral.config import ConfigInfo, get_yaml_config
from pysiral.errorhandler import ErrorStatus
from pysiral.grid import GridDefinition
from pysiral.logging import DefaultLoggingClass
from pysiral.l2data import L2iNCFileImport
from pysiral.output import L3SDataNC, OutputHandlerBase
from pysiral.flag import ORCondition
from pysiral.surface_type import SurfaceType
import numpy as np
import os


# %% Level 3 Processor

class Level3Processor(DefaultLoggingClass):

    def __init__(self, product_def):
        super(Level3Processor, self).__init__(self.__class__.__name__)
        self.error = ErrorStatus(caller_id=self.__class__.__name__)
        self._job = product_def
        self._l3_progress_percent = 0.0

    def process_l2i_files(self, l2i_files):

        # Store l2i_files
        self._l2i_files = l2i_files

        # Initialize the stack for the l2i orbit files
        self.log.info("Initialize l2i data stack")
        stack = L2DataStack(self._job.grid, self._job.l2_parameter)
        stack.initialize(len(l2i_files))

        # Parse all orbit files and add to the stack
        for i, l2i_file in enumerate(l2i_files):
            self._log_progress(i)
            l2data = L2iNCFileImport(l2i_file)
            stack.append(l2data)

        # Initialize the data grid
        self.log.info("Initialize l3 data grid")
        grid = L3DataGrid(self._job.grid)
        grid.init_parameter_fields(self._job.l2_parameter, "l2")
        grid.init_parameter_fields(self._job.l3_parameter, "l3")
        grid.init_mandatory_parameter_fields()
        grid.calculate_longitude_latitude_fields()

        # Average level-2 parameter for each grid cell
        grid.set_l2i_stack(stack)
        self.log.info("Compute masks and mandatory grid statistics")
        grid.compute_l3_mandatory_parameter()

        self.log.info("Compute l2i parameter averages")
        grid.compute_l2_parameter_averages()

        # Get the level-3 parameter
        self.log.info("Compute level-3 ouput parameter")
        grid.compute_l3_output_parameter()

        # Set parameters nan if freeboard is nan
        # (list in output definition file)
        self.log.info("Apply data masks")
        grid.set_freeboard_nan_mask(self._job.freeboard_nan_mask_targets)
        grid.set_sic_mask(self._job.sea_ice_concentration_mask_targets)

        # Get the metadata information from the L2 stack
        self.log.info("Compile metadata")
        l3_metadata = L3MetaData()
        l3_metadata.get_missions_from_stack(stack)
        l3_metadata.get_data_period_from_stack(stack)
        l3_metadata.get_projection_parameter(self._job.grid)

        # Write grid to L3S netcdf file:
#        self.log.info("Exporting file")
#        output = L3SDataNC()
#        output.set_metadata(l3_metadata)
#        output.set_export_folder(self._job.export_folder)
#        output.export(grid)

    def _log_progress(self, i):
        """ Concise logging on the progress of l2i stack creation """
        n = len(self._l2_files)
        progress_percent = float(i)/float(n-1)*100.
        current_reminder = np.mod(progress_percent, 10)
        last_reminder = np.mod(self._l3_progress_percent, 10)
        if last_reminder > current_reminder:
            self.log.info("Creating l2i orbit stack: %3g%% (%g of %g)" % (
                          progress_percent-current_reminder, i+1, n))
        self._l3_progress_percent = progress_percent


# %% Data Containers

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
                try:
                    data = getattr(orbit, parameter_name)
                    self.stack[parameter_name][y][x].append(data[i])
                except:
                    pass


class L3DataGrid(DefaultLoggingClass):
    """
    Container for computing gridded data sets based on a l2i data stack
    (averaged l2i parameter, grid cell statistics)
    """

    def __init__(self, griddef):

        super(L3DataGrid, self).__init__(self.__class__.__name__)

        # Grid size definition
        self.griddef = griddef

        # Shortcut to the surface type flag dictionalry
        self._surface_type_dict = SurfaceType.SURFACE_TYPE_DICT

        # Name and data type of mandatory surface type statistics
        self._surface_type_l3par = {
            "n_total_waveforms": "i4",
            "n_valid_waveforms": "i4",
            "valid_fraction": "f4",
            "lead_fraction": "f4",
            "ice_fraction": "f4",
            "is_land": "i2"}

        # List of level-2 parameter
        # (gridded parameter that are already in l2i)
        self._l2_parameter = None

        # list of level-3 parameter
        # (grid specific parameter, e.g. surface type statistics)
        self._l3_parameter = None

        # list of stacked l2 parameters for each grid cell
        self._l2 = None

        # container for gridded parameters
        self._l3 = {}

    def init_parameter_fields(self, parameter_names, level):
        """ Initialize output parameter fields """
        setattr(self, "_"+level+"_parameter", parameter_names)
        shape = self.grid_shape
        for parameter_name in parameter_names:
            msg = "Adding %s parameter: %s" % (level, parameter_name)
            self.log.info(msg)
            self._l3[parameter_name] = np.ndarray(
                shape=shape, dtype='f4')*np.nan

    def init_mandatory_parameter_fields(self):
        """
        Initialize mandatory parameter field that may be needed for masking
        and thus will be calculated whether they are included in the
        output files or not (computational cost is small)

        Note: These parameter fields will be computed separately from the
        loops over the l2/l3 parameter loops. But as the fields can
        be also specified in the output format definition, these will be
        than ignored when looping over all output parameters
        """

        # grid output dimensions
        shape = self.grid_shape

        # XXX: There needs to be a better handling of data types
        #       (requires better definition of l3 output file format)
        self.lon = np.ndarray(shape=shape, dtype='f4')*np.nan
        self.lat = np.ndarray(shape=shape, dtype='f4')*np.nan

        self.log.info("Adding parameter: lon")
        self.log.info("Adding parameter: lat")

        # Surface type statistics
        for surface_type_statistics_par in self._surface_type_l3par.keys():
            # Check if already created, which will be the case if
            # the parameter is in the l3 output definition
            if surface_type_statistics_par in self._l3_parameter:
                continue
            self.log.info("Adding parameter: %s" % surface_type_statistics_par)
            dtype = self._surface_type_l3par[surface_type_statistics_par]
            self._l3[surface_type_statistics_par] = np.ndarray(
                shape=shape, dtype=dtype)

        # Sea Ice Concentration (for potential masking)
        if "sea_ice_concentration" not in self._l2_parameter:
            self.log.info("Adding parameter: sea_ice_concentration")
            self._l3["sea_ice_concentration"] = np.ndarray(
                shape=shape, dtype='f4')*np.nan

    def calculate_longitude_latitude_fields(self):

        # TODO: Move all this into a grid definition container
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

        self._l3["lon"] = lon
        self._l3["lat"] = lat

    def set_l2i_stack(self, l2i_stack):
        """
        Set the l2i data stack
        (list of all individual measurements per grid cell
         grouped by parameter)
        """
        # Input validation
        if not isinstance(l2i_stack, L2DataStack):
            msg = "Input must be of type pysiral.l3proc.L2DataStack, was %s"
            msg = msg % type(l2i_stack)
            raise ValueError(msg)
        self._l2 = l2i_stack

    def compute_l2_parameter_averages(self):
        """
        Compute averages of all l2i parameter for each grid cell.
        The list of l2i parameter is from the output format definition

        No averages are computed for grid cells that are tagged with
        a land flag.
        """
        # Loop over all grid cells
        # XXX: Is there a better way?
        for xi in self.grid_xi_range:
            for yj in self.grid_yj_range:

                # Exclude land (or near land grid cells)
                if self._l3["is_land"][xi, yj]:
                    continue

                for name in self._l2_parameter:
                    data = np.array(self._l2.stack[name][yj][xi])

                    # nanmean needs at least 2 valid items
                    if len(np.where(np.isfinite(data))[0]) > 1:
                        self._l3[name][yj, xi] = np.nanmean(data)

    def compute_l3_mandatory_parameter(self):
        """
        Wrapper method for computing the surface type statistics for
        each grid cell
        """
        for xi in self.grid_xi_range:
            for yj in self.grid_yj_range:
                for l3_parameter_name in self._l3_parameter:
                    self._compute_surface_type_grid_statistics(xi, yj)

    def compute_l3_output_parameter(self):
        """
        Compute level-3 parameter for each grid cell. A parameter is
        classified as level-3 if it only exists on the grid cell level
        (e.g. total number of waveforms). Parameters that are averaged
        from  the l2i orbits, are called level-2 in the terminology
        of pysiral.Level2Processor
        """
        # Loop over grid items
        for xi in self.grid_xi_range:
            for yj in self.grid_yj_range:
                for l3_parameter_name in self._l3_parameter:
                    # level-3 parameter can to be computed
                    result = self.get_l3_parameter(l3_parameter_name, xi, yj)
                    if result is not None:
                        self._l3[l3_parameter_name][yj, xi] = result

    def get_l3_parameter(self, l3_parameter_name, xi, yj):
        """
        Compution of all level-3 parameter for a given grid cell
        Since level-3 parameter may be computed in very different ways,
        this method is supplied with the grid index and the name of the
        parameter and then redirects to the corresponding computation
        method.

        It returns the parameter value, which will be ignored if None. This
        is the case for the mandatory level-3 parameter (surface type
        statistics), which are computed in a seperate method and can safely
        be ignored here.
        """
        # Surface type based parameter are computed anyway, skip
        if l3_parameter_name in self._surface_type_l3par:
            return None
        # No other l3 parameter computations at the moment
        else:
            raise ValueError("Unknown l3 parameter name: %s" %
                             l3_parameter_name)

    def set_freeboard_nan_mask(self, nan_masks_targets):
        """
        Apply the freeboard nan mask to a selected number of parameters
        see setting/l3/*.yaml for details
        """
        frb_is_nan = np.where(np.isnan(self._l3["freeboard"]))
        for target in nan_masks_targets:
            self.log.info("freeboard nan mask to %s" % target)
            self._l3[target][frb_is_nan] = np.nan

    def set_sic_mask(self, nan_masks_targets, sic_threshold=5.0):
        """
        Apply the sea ice concentration mask (0 or nan) to a selected number
        of parameters see setting/l3/*.yaml for details
        """
        sic_is_nan = np.isnan(self._l3["sea_ice_concentration"])
        sic_is_zero = self._l3["sea_ice_concentration"] <= sic_threshold
        sic_mask = np.logical_or(sic_is_nan, sic_is_zero)
        for target in nan_masks_targets:
            self.log.info("sea ice concentration mask to %s" % target)
            self._l3[target][sic_mask] = np.nan

    def _compute_surface_type_grid_statistics(self, xi, yj):
        """
        Computes the mandatory surface type statistics for a given
        grid index based on the surface type stack flag

        The current list
          - is_land (land flag exists in l2i stack)
          - n_total_waveforms (size of l2i stack)
          - n_valid_waveforms (tagged as either lead or sea ice )
          - valid_fraction (n_valid/n_total)
          - lead_fraction (n_leads/n_valid)
          - ice_fraction (n_ice/n_valid)
        """
        surface_type = np.array(self._l2.stack["surface_type"][yj][xi])

        # Stack can be empty
        if len(surface_type) == 0:
            return

        stflags = self._surface_type_dict

        # Create a land flag
        is_land = len(np.where(surface_type == stflags["land"])[0] > 0)
        self._l3["is_land"][xi, yj] = is_land

        # Compute total waveforms in grid cells
        n_total_waveforms = len(surface_type)
        self._l3["n_total_waveforms"][yj, xi] = n_total_waveforms

        # Compute valid wavefords
        # Only positively identified waveforms (either lead or ice)
        # XXX: what about polynay and ocean?
        valid_waveform = ORCondition()
        valid_waveform.add(surface_type == stflags["lead"])
        valid_waveform.add(surface_type == stflags["sea_ice"])
        n_valid_waveforms = valid_waveform.num
        self._l3["n_valid_waveforms"][yj, xi] = n_valid_waveforms

        # Fractions of leads on valid_waveforms
        try:
            valid_fraction = float(n_valid_waveforms)/float(n_total_waveforms)
        except:
            valid_fraction = np.nan
        self._l3["valid_fraction"][yj, xi] = valid_fraction

        # Fractions of leads on valid_waveforms
        n_leads = len(np.where(surface_type == stflags["lead"])[0])
        try:
            lead_fraction = float(n_leads)/float(n_valid_waveforms)
        except:
            lead_fraction = np.nan
        self._l3["lead_fraction"][yj, xi] = lead_fraction

        # Fractions of leads on valid_waveforms
        n_ice = len(np.where(surface_type == stflags["sea_ice"])[0])
        try:
            ice_fraction = float(n_ice)/float(n_valid_waveforms)
        except:
            ice_fraction = np.nan
        self._l3["ice_fraction"][yj, xi] = ice_fraction

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
        parameter_list = list(self._l2_parameter)
        parameter_list.extend(self._l3_parameter)
        parameter_list.append("lon")
        parameter_list.append("lat")
        return parameter_list

    @property
    def dimdict(self):
        from collections import OrderedDict
        dimdict = OrderedDict([("time", 1),
                               ("lat", self.griddef.extent.numx),
                               ("lon", self.griddef.extent.numy)])
        return dimdict

    @property
    def grid_shape(self):
        return (self.griddef.extent.numx, self.griddef.extent.numy)


class L3MetaData(object):

    """
    Container for L3S Metadata information
    (see property attribute_list for a list of attributes)
    """

    _attribute_list = [
        "mission_ids", "start_time", "stop_time", "grid_name", "period_label",
        "pysiral_version", "projection_str", "grid_tag", "resolution_tag",
        "hemisphere"]

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
        self.set_attribute("hemisphere", griddef.hemisphere)
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
        return self.start_time

    @property
    def stop_period(self):
        return self.stop_time

    def set_attribute(self, tag, value):
        if tag not in self.attribute_list:
            raise ValueError("Unknown attribute: ", tag)
        setattr(self, tag, value)


class Level3OutputHandler(OutputHandlerBase):

    # Some fixed parameters for this class
    default_file_location = ["settings", "outputdef", "l3_default.yaml"]
    subfolder_tags = ["year"]
    applicable_data_level = 3

    def __init__(self, output_def="default", base_directory="l3proc_default",
                 overwrite_protection=True):

        if output_def == "default":
            output_def = self.default_output_def_filename

        super(Level3OutputHandler, self).__init__(output_def)
        self.error.caller_id = self.__class__.__name__
        self.log.name = self.__class__.__name__
        self.overwrite_protection = overwrite_protection
        self._init_product_directory(base_directory)

    def _init_product_directory(self, base_directory_or_id):
        """ Initializes the product directory. If `base_directory` is already
        a directory, it is used as is. Else, it is treated as subfolder of
        the default pysiral product directory. Product level id and
        overwrite protection subfolders are added for both options"""
        # argument is directory
        if os.path.isdir(base_directory_or_id):
            basedir = base_directory_or_id
        # argument is id
        else:
            basedir = self.pysiral_config.local_machine.product_repository
            basedir = os.path.join(basedir, base_directory_or_id)
        # add product level subfolder
        basedir = os.path.join(basedir, self.product_level_subfolder)
        # optional (subfolder with current time)
        if self.overwrite_protection:
            basedir = os.path.join(basedir, self.now_directory)
        # set the directory
        self._set_basedir(basedir)

    @property
    def default_output_def_filename(self):
        pysiral_config = ConfigInfo()
        local_settings_path = pysiral_config.pysiral_local_path
        return os.path.join(local_settings_path, *self.default_file_location)


class Level3GridDefinition(GridDefinition):
    """ This is a variation of GridDefinition with a mandatory link to
    a griddef yaml file"""

    def __init__(self, l3_settings_file):
        super(Level3GridDefinition, self).__init__(self)
        self.set_from_griddef_file(l3_settings_file)


class Level3ProductDefinition(DefaultLoggingClass):

    def __init__(self, l3_settings_file, grid, output):
        """ Container for the Level3Processor settings

        Arguments:
            l3_settings_file (str): Full filename to l3 settings file
            grid (pysiral.grid.GridDefinition): Output grid class
            output (Level-3 compliant output handler from pysiral.output)
        """
        super(Level3ProductDefinition, self).__init__(self.__class__.__name__)
        self.error = ErrorStatus(caller_id=self.__class__.__name__)
        self._l3_settings_file = l3_settings_file
        self._output = [output]
        self._grid = grid
        self._parse_l3_settings()

        # Report settings to log handler
        self.log.info("Output grid id: %s" % str(self._grid.grid_id))
        for output in self._output:
            msg = "L3 product directory (%s): %s"
            msg = msg % (str(output.id), str(output.basedir))
            self.log.info(msg)

    def _parse_l3_settings(self):
        self.log.info("Parsing settings: %s" % str(self._l3_settings_file))
        try:
            self._l3 = get_yaml_config(self._l3_settings_file)
        except Exception, msg:
            self.error.add_error("l3settings-parser-error", msg)
            self.error.raise_on_error()

    def validate(self):
        pass

    @property
    def grid(self):
        return self._grid

    @property
    def outputs(self):
        return self._output

    @property
    def n_outputs(self):
        return len(self._output)

    @property
    def l3def(self):
        return self._l3
