# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 14:04:27 2015

@author: Stefan
"""

from pysiral.config import (ConfigInfo, get_yaml_config, SENSOR_NAME_DICT,
                            MISSION_NAME_DICT)
from pysiral.errorhandler import ErrorStatus
from pysiral.grid import GridDefinition
from pysiral.logging import DefaultLoggingClass
from pysiral.l2data import L2iNCFileImport
from pysiral.output import OutputHandlerBase, Level3Output
from pysiral.flag import ORCondition
from pysiral.surface_type import SurfaceType

from datetime import datetime
import numpy as np
import sys
import os
import re


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
        stack = L2iDataStack(self._job.grid, self._job.l2_parameter)

        self.log.info("Parsing products (prefilter active: %s)" % (
                str(self._job.l3def.l2i_prefilter.active)))

        # Parse all orbit files and add to the stack
        for i, l2i_file in enumerate(l2i_files):

            self._log_progress(i)

            # Parse l2i source file
            l2i = L2iNCFileImport(l2i_file)

            # Prefilter l2i product
            # Note: In the l2i product only the minimum set of nan are used
            #       for different parameters (e.g. the radar freeboard mask
            #       does not equal the thickness mask). This leads to
            #       inconsistent results during gridding and therefore it is
            #       highly recommended to harmonize the mask for thickness
            #       and the different freeboard levels
            prefilter = self._job.l3def.l2i_prefilter
            if prefilter.active:
                l2i.transfer_nan_mask(prefilter.nan_source,
                                      prefilter.nan_targets)
            # Add to stack
            stack.add(l2i)

        # Initialize the data grid
        self.log.info("Initialize l3 data grid")
        l3 = L3DataGrid(self._job, stack)

        # Write output(s)
        for output_handler in self._job.outputs:
            output = Level3Output(l3, output_handler)
            self.log.info("Write %s product: %s" % (output_handler.id,
                                                    output.filename))

    def _log_progress(self, i):
        """ Concise logging on the progress of l2i stack creation """
        n = len(self._l2i_files)
        progress_percent = float(i)/float(n-1)*100.
        current_reminder = np.mod(progress_percent, 10)
        last_reminder = np.mod(self._l3_progress_percent, 10)
        if last_reminder > current_reminder:
            self.log.info("Creating l2i orbit stack: %3g%% (%g of %g)" % (
                          progress_percent-current_reminder, i+1, n))
        self._l3_progress_percent = progress_percent


# %% Data Containers

class L2iDataStack(DefaultLoggingClass):

    def __init__(self, griddef, l2_parameter):
        """ A container for stacking l2i variables (geophysical paramters
        at sensor resolution) in L3 grid cells. For each parameters
        a (numx, numy) array is created, with an list containing all
        l2i data points that fall into the grid cell area. This list
        can be averaged or used for histogram computation in later stages
        of the Level-3 processor.

        Args:
            griddef (obj): pysiral.grid.GridDefinition or inheritated objects
            l2_parameter (str list): list of l2i parameter names

        Returns:
            class instance

        """
        super(L2iDataStack, self).__init__(self.__class__.__name__)

        # Grid Definition Type
        self.griddef = griddef

        # A list of level-2 parameters to be stacked
        self.l2_parameter = l2_parameter

        # Statistics
        self._n_records = 0
        self._l2i_count = 0
        self.start_time = []
        self.stop_time = []
        self.mission = []

        # Flags
        self._has_surface_type = False

        # Save global attricbutes from l2i (will be overwritten for each
        # l2i file, but it is assumed that general information, e.g.
        # on auxdata remains the same)
        self._l2i_info = None

        # Create parameter stacks
        self._initialize_stacks()

    def _initialize_stacks(self):
        """ Create all data stacks, content will be added sequentially
        with `add` method """

        # Stack dictionary that will hold the data
        self.stack = {}

        # init stack arrays
        # surface_type is mandatory for level-3 parameters
        # (e.g. n_total_wave_forms, lead_fraction, ...)
        self.stack["surface_type"] = self.parameter_stack

        # create a stack for each l2 parameter
        for parameter_name in self.l2_parameter:
            self.stack[parameter_name] = self.parameter_stack

    def add(self, l2i):
        """ Add a l2i data object to the stack

        Args:
            l2i (obj): l2i object (currently: pysiral.l2data.L2iNCFileImport)

        Returns:
            None
        """

        # Save the metadata from the orbit data
        self.start_time.append(l2i.timestamp[0])
        self.stop_time.append(l2i.timestamp[-1])
        self.mission.append(l2i.mission)
        self._l2i_count += 1
        self._n_records += l2i.n_records

        self._l2i_info = l2i.info

        # Get projection coordinates for l2i locations
        xi, yj = self.griddef.grid_indices(l2i.longitude, l2i.latitude)

        # Stack the l2 parameter in the corresponding grid cells
        for i in np.arange(l2i.n_records):
            # Add the surface type per default
            # (will not be gridded, therefore not in list of l2 parameter)
            x, y = int(xi[i]), int(yj[i])
            self.stack["surface_type"][y][x].append(l2i.surface_type[i])
            for parameter_name in self.l2_parameter:
                try:
                    data = getattr(l2i, parameter_name)
                    self.stack[parameter_name][y][x].append(data[i])
                except:
                    pass

    @property
    def n_total_records(self):
        return self._n_records

    @property
    def l2i_count(self):
        return self._l2i_count

    @property
    def parameter_stack(self):
        dimx, dimy = self.griddef.extent.numx, self.griddef.extent.numy
        return [[[] for _ in range(dimx)] for _ in range(dimy)]

    @property
    def l2i_info(self):
        return self._l2i_info


class L3DataGrid(DefaultLoggingClass):
    """
    Container for computing gridded data sets based on a l2i data stack
    (averaged l2i parameter, grid cell statistics)
    """

    def __init__(self, job, stack):

        super(L3DataGrid, self).__init__(self.__class__.__name__)

        self.error = ErrorStatus(caller_id=self.__class__.__name__)

        # Grid size definition
        self._griddef = job.grid
        self._l3def = job.l3def
        self._period = job.period

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

        self._metadata = None

        self.init_parameter_fields(job.l2_parameter, "l2")
        self.init_parameter_fields(job.l3_parameter, "l3")
        self.init_mandatory_parameter_fields()
        self.calculate_longitude_latitude_fields()

        # Average level-2 parameter for each grid cell
        self.set_l2i_stack(stack)
        self.log.info("Compute masks and mandatory grid statistics")
        self.compute_l3_mandatory_parameter()

        self.log.info("Grid l2i parameter")
        self.grid_l2_parameter()

        # Get the level-3 parameter
        self.log.info("Compute level-3 ouput parameter")
        self.compute_l3_output_parameter()

        # Set parameters nan if freeboard is nan
        # (list in output definition file)
        self.log.info("Apply data masks")
        for mask_def in job.l3_masks:
            self.mask_l3(mask_def)

        # Get the metadata information from the L2 stack
        self.log.info("Compile metadata")
        l3_metadata = L3MetaData()
        l3_metadata.get_missions_from_stack(stack)
        # Actual data coverage
        l3_metadata.get_data_period_from_stack(stack)
        # Requested time coverage (might not be the actual coverage)
        l3_metadata.get_time_coverage_from_period(self._period)
        l3_metadata.get_auxdata_infos(stack.l2i_info)
        l3_metadata.get_projection_parameter(job.grid)
        self.set_metadata(l3_metadata)

    def set_metadata(self, metadata):
        self._metadata = metadata

    def get_attribute(self, attribute_name, *args):
        """ Return a string for a given attribute name. This method is
        required for the output data handler """

        try:
            attr_getter = getattr(self, "_get_attr_"+attribute_name)
            attribute = attr_getter(*args)
            return attribute
        except AttributeError:
            return "attr_unavailable"
        except Exception, msg:
            print msg
            sys.exit(1)

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

        self.log.info("Adding parameter: longitude")
        self.log.info("Adding parameter: latitude")

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
        """ Geographic coordinates from GridDefinition """
        lon, lat = self.griddef.get_grid_coordinates()
        self._l3["longitude"] = lon
        self._l3["latitude"] = lat

    def set_l2i_stack(self, l2i_stack):
        """ Set the l2i data stack (list of all individual measurements
        per grid cell grouped by parameter) """
        # Input validation
        if not isinstance(l2i_stack, L2iDataStack):
            msg = "Input must be of type pysiral.l3proc.L2DataStack, was %s"
            msg = msg % type(l2i_stack)
            raise ValueError(msg)
        self._l2 = l2i_stack

    def grid_l2_parameter(self):
        """ Compute averages of all l2i parameter for each grid cell.
        The list of l2i parameter is from the output format definition
        No averages are computed for grid cells that are tagged with
        a land flag. """

        settings = self.l3def.grid_settings

        # Loop over all grid cells
        # XXX: Is there a better way?
        for xi in self.grid_xi_range:
            for yj in self.grid_yj_range:

                # Exclude land (or near land grid cells)
                if self._l3["is_land"][xi, yj] and settings.no_land_cells:
                    continue

                for name in self._l2_parameter:
                    data = np.array(self._l2.stack[name][yj][xi])

                    # nanmean needs at least 2 valid items
                    valid = np.where(np.isfinite(data))[0]
                    if len(valid) < settings.minimum_valid_grid_points:
                        continue

                    grid_method = self.l3def.l2_parameter[name].grid_method
                    if grid_method == "average":
                        self._l3[name][yj, xi] = np.nanmean(data)
                    elif grid_method == "average_uncertainty":
                        value = np.abs(np.sqrt(1./np.sum(data[valid])))
                        self._l3[name][yj, xi] = value
                    else:
                        msg = "Invalid grid method (%s) for %s"
                        msg = msg % (str(grid_method), name)
                        self.error.add_error("invalid-l3def", msg)
                        self.error.raise_on_error()

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

    def mask_l3(self, mask_def):
        """ Apply a parametrized mask to level 3 data """

        # Get the source parameter
        source = self._l3[mask_def.source]

        # Compute the masking condition
        conditions = mask_def.condition.split(";")
        n_conditions = len(conditions)

        if n_conditions == 0:
            msg = "Missing condition in %s" % str(mask_def)
            self.error.add_error("invalid-l3mask-def", msg)
            self.error.raise_on_error()

        # Start with the first (and maybe only condition)
        filter_mask = self._get_l3_mask(source, conditions[0], mask_def)

        # Add conditions
        if n_conditions >= 2:
            for i in range(1, n_conditions):
                new_filter = self._get_l3_mask(source, conditions[i], mask_def)
                if mask_def.connect_conditions == "or":
                    filter_mask = np.logical_or(filter_mask, new_filter)
                elif mask_def.connect_conditions == "and":
                    filter_mask = np.logical_and(filter_mask, new_filter)
                else:
                    msg = "Invalid l3 mask operation: %s"
                    msg = msg % mask_def.connect_conditions
                    self.error.add_error("invalid-l3mask-def", msg)
                    self.error.raise_on_error()

        self.log.info("Apply l3 mask: %s" % mask_def.branchName())

        # Apply mask
        masked_indices = np.where(filter_mask)
        for target in mask_def.targets:
            self._l3[target][masked_indices] = np.nan

    def _get_l3_mask(self, source_param, condition, options):
        """ Returna bool array based on a parameter and a predefined
        masking operation """
        if condition.strip() == "is_nan":
            return np.isnan(source_param)
        elif condition.strip() == "is_zero":
            return np.array(source_param <= 1.0e-9)
        elif condition.strip() == "is_smaller":
            return np.array(source_param < options.is_smaller_threshold)
        else:
            msg = "Unknown condition in l3 mask: %s" % condition
            self.error.add_error("invalid-l3mask-condition", msg)
            self.error.raise_on_error()

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
        try:
            parameter = self._l3[name]
        except KeyError:
            parameter = np.full(np.shape(self._l3["longitude"]), np.nan)
            self.log.warn("Parameter not availabe: %s" % name)
        except Exception, msg:
            print msg
            sys.exit(1)
        return parameter

    def _get_attr_source_mission_id(self, *args):
        mission_ids = self.metadata.mission_ids
        if args[0] == "uppercase":
            mission_ids = mission_ids.upper()
        return mission_ids

    def _get_attr_source_mission_name(self, *args):
        ids = self.metadata.mission_ids
        names = ",".join([MISSION_NAME_DICT[m] for m in ids.split(",")])
        return names

    def _get_attr_grid_id(self, *args):
        grid_id = self.griddef.grid_id
        if args[0] == "uppercase":
            grid_id = grid_id.upper()
        return grid_id

    def _get_attr_source_mission_sensor(self, *args):
        mission_sensor = self.metadata.mission_sensor
        if args[0] == "uppercase":
            mission_sensor = mission_sensor.upper()
        return mission_sensor

    def _get_attr_source_hemisphere(self, *args):
        return self.hemisphere

    def _get_attr_startdt(self, dtfmt):
        return self.metadata.start_period.strftime(dtfmt)

    def _get_attr_stopdt(self, dtfmt):
        return self.info.stop_time.strftime(dtfmt)

    def _get_attr_geospatial_lat_min(self, *args):
        latitude = self._l3["latitude"]
        return self._get_attr_geospatial_str(np.nanmin(latitude))

    def _get_attr_geospatial_lat_max(self, *args):
        latitude = self._l3["latitude"]
        return self._get_attr_geospatial_str(np.nanmax(latitude))

    def _get_attr_geospatial_lon_min(self, *args):
        longitude = self._l3["longitude"]
        return self._get_attr_geospatial_str(np.nanmin(longitude))

    def _get_attr_geospatial_lon_max(self, *args):
        longitude = self._l3["longitude"]
        return self._get_attr_geospatial_str(np.nanmax(longitude))

    def _get_attr_geospatial_str(self, value):
        return "%.4f" % value

    def _get_attr_source_auxdata_sic(self, *args):
        return self.metadata.source_auxdata_sic

    def _get_attr_source_auxdata_snow(self, *args):
        return self.metadata.source_auxdata_snow

    def _get_attr_source_auxdata_sitype(self, *args):
        return self.metadata.source_auxdata_sitype

    def _get_attr_utc_now(self, *args):
        return datetime.now().isoformat()

    def _get_attr_time_coverage_start(self, *args):
        datetime = self.metadata.time_coverage_start
        if re.match("%", args[0]):
            time_string = datetime.strftime(args[0])
        else:
            time_string = datetime.isoformat()
        return time_string

    def _get_attr_time_coverage_end(self, *args):
        datetime = self.metadata.time_coverage_end
        if re.match("%", args[0]):
            time_string = datetime.strftime(args[0])
        else:
            time_string = datetime.isoformat()
        return time_string

    def _get_attr_time_coverage_duration(self, *args):
        return self.metadata.time_coverage_duration

    def flipud(self):
        for parameter in self.parameters:
            setattr(self, parameter, np.flipud(getattr(self, parameter)))

    @property
    def metadata(self):
        return self._metadata

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

    @property
    def griddef(self):
        return self._griddef

    @property
    def l3def(self):
        return self._l3def

    @property
    def hemisphere(self):
        return self.metadata.hemisphere


class L3MetaData(object):

    """
    Container for L3S Metadata information
    (see property attribute_list for a list of attributes)
    """

    _attribute_list = [
        "mission_ids", "start_time", "stop_time", "grid_name", "period_label",
        "time_coverage_start", "time_coverage_end", "time_coverage_duration",
        "pysiral_version", "projection_str", "grid_tag", "resolution_tag",
        "hemisphere", "mission_sensor", "source_auxdata_sic",
        "source_auxdata_sitype", "source_auxdata_snow"]

    def __init__(self):
        # Init all fields
        for field in self.attribute_list:
            setattr(self, field, None)

    def get_missions_from_stack(self, stack):
        """
        Get a list of missions that went into the stack
        (must be a list, since multi-mission grids are supported)
        """
        missions = np.unique(stack.mission)
        mission_sensor = [SENSOR_NAME_DICT[mission] for mission in missions]
        self.set_attribute("mission_ids", ",".join(missions))
        self.set_attribute("mission_sensor", ",".join(mission_sensor))

    def get_data_period_from_stack(self, stack):
        """ Get the first and last timestamp """
        self.set_attribute("start_time", np.amin(stack.start_time))
        self.set_attribute("stop_time", np.amax(stack.stop_time))
        # XXX: Only monthly periods are currently supported
        self.set_attribute("period_label", self.start_time.strftime("%B %Y"))

    def get_time_coverage_from_period(self, period):
        """ Get the start and end of requested data period """
        self.set_attribute("time_coverage_start", period.start_dt)
        self.set_attribute("time_coverage_end", period.stop_dt)
        self.set_attribute("time_coverage_duration",
                           period.base_duration_isoformat)

    def get_auxdata_infos(self, l2i_info):
        """ Get information on auxiliary data sources from l2i global
        attributes """
        self.set_attribute("source_auxdata_sic", l2i_info.source_sic)
        self.set_attribute("source_auxdata_sitype", l2i_info.source_sitype)
        self.set_attribute("source_auxdata_snow", l2i_info.source_snow)

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

    def get_filename_from_data(self, l3):
        """ Return the filename for a defined level-2 data object
        based on tag filenaming in output definition file """
        filename_template = self.output_def.filenaming
        return self.fill_template_string(filename_template, l3)

    def get_directory_from_data(self, l3, create=True):
        """ Return the output directory based on information provided
        in an l2 data object """
        directory = self._get_directory_from_dt(l3.metadata.start_time)
        if create:
            self._create_directory(directory)
        return directory

    def get_fullpath_from_data(self, l3):
        """ Return export path and filename based on information
        provided in the l2 data object """
        export_directory = self.get_directory_from_data(l3)
        export_filename = self.get_filename_from_data(l3)
        return os.path.join(export_directory, export_filename)

    def get_global_attribute_dict(self, l3):
        attr_dict = {}
        for attr_name in self.output_def.global_attributes.iterkeys():
            attr_template = self.output_def.global_attributes[attr_name]
            attribute = self.fill_template_string(attr_template, l3)
            attr_dict[attr_name] = attribute
        return attr_dict

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

    def __init__(self, l3_settings_file, grid, output, period):
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
        self._period = period
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

    @property
    def period(self):
        return self._period

    @property
    def l3_masks(self):
        """ Return a sorted list of the masks applied to level 3 data """
        mask_names = sorted(self.l3def.l3_masks.keys(branch_mode="only"))
        return [self.l3def.l3_masks[name] for name in mask_names]

    @property
    def l2_parameter(self):
        """ Extract a list of paramter names to be extracted from
        l2i product files """
        l2_parameter = sorted(self.l3def.l2_parameter.keys(branch_mode="only"))
        if type(l2_parameter) is not list:
            msg = "Missing or invalid parameter in l3 settings: %s"
            msg = msg % "root.l2_parameter"
            self.error.add_error("invalid-l3-settings", msg)
            self.error.raise_on_error()
        return l2_parameter

    @property
    def l3_parameter(self):
        """ Extract a list of paramter names to be extracted from
        l2i product files """
        l3_parameter = self.l3def.l3_parameter
        if type(l3_parameter) is not list:
            msg = "Missing or invalid parameter in l3 settings: %s"
            msg = msg % "root.l3_parameter"
            self.error.add_error("invalid-l3-settings", msg)
            self.error.raise_on_error()
        return self.l3def.l3_parameter
