# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 14:04:27 2015

@author: Stefan
"""
import itertools
import re
import sys
import uuid
from collections import OrderedDict
from datetime import date, datetime
from pathlib import Path

import numpy as np
from loguru import logger
from scipy import stats
from scipy.ndimage import maximum_filter
from xarray import open_dataset

from pysiral import __version__, get_cls, psrlcfg
from pysiral.core import DefaultLoggingClass
from pysiral.core.config import get_yaml_config
from pysiral.core.errorhandler import ErrorStatus
from pysiral.core.flags import SURFACE_TYPE_DICT, ORCondition
from pysiral.core.output import Level3Output, OutputHandlerBase
from pysiral.grid import GridDefinition
from pysiral.l2data import L2iNCFileImport
from pysiral.mask import L3Mask
from pysiral.sit import frb2sit_errprop

# %% Level 3 Processor

class Level3Processor(DefaultLoggingClass):

    def __init__(self, product_def):
        super(Level3Processor, self).__init__(self.__class__.__name__)
        self.error = ErrorStatus(caller_id=self.__class__.__name__)
        self._job = product_def
        self._l3_progress_percent = 0.0
        self._l2i_files = None
        self._period = None

    def process_l2i_files(self, l2i_files, period):
        """
        The main call for the Level-3 processor
        TODO: Needs organization
        :param l2i_files:
        :param period:
        :return:
        """

        # Store l2i_files
        self._l2i_files = l2i_files

        # Store
        self._period = period

        # Initialize the stack for the l2i orbit files
        logger.info("Initialize l2i data stack")
        stack = L2iDataStack(self._job.grid, self._job.l2_parameter)

        logger.info("Parsing products (prefilter active: %s)" % (str(self._job.l3def.l2i_prefilter.active)))

        # Parse all orbit files and add to the stack
        for i, l2i_file in enumerate(l2i_files):

            self._log_progress(i)

            # Parse l2i source file
            try:
                l2i = L2iNCFileImport(l2i_file)
            except AttributeError:
                logger.warning("Attribute Error encountered in %s" % l2i_file)
                continue

            # Apply the orbit filter (for masking descending or ascending orbit segments)
            # NOTE: This tag may not be present in all level-3 settings files, as it has
            #       been added as a test case
            # TODO: Create a configurable processor item
            orbit_filter = self._job.l3def.get("orbit_filter")
            if orbit_filter is not None:
                self.apply_orbit_filter(l2i, orbit_filter)

            # Apply the orbit filter (for masking descending or ascending orbit segments)
            # NOTE: This tag may not be present in all level-3 settings files, as it has
            #       been added as a test case
            # TODO: Create a configurable processor item
            miz_filter = self._job.l3def.get("miz_filter")
            if miz_filter is not None:
                self.apply_miz_filter(l2i, miz_filter)

            # Prefilter l2i product
            # Note: In the l2i product only the minimum set of nan are used
            #       for different parameters (e.g. the radar freeboard mask
            #       does not equal the thickness mask). This leads to
            #       inconsistent results during gridding and therefore it is
            #       highly recommended to harmonize the mask for thickness
            #       and the different freeboard levels
            # TODO: Create a configurable processor item
            prefilter = self._job.l3def.l2i_prefilter
            if prefilter.active:
                l2i.transfer_nan_mask(prefilter.nan_source, prefilter.nan_targets)
            # Add to stack
            stack.add(l2i)

        # Initialize the data grid
        logger.info("Initialize l3 data grid")
        l3 = L3DataGrid(self._job, stack, period)

        # Apply the processing items
        self._apply_processing_items(l3)

        # Write output(s)
        for output_handler in self._job.outputs:
            output = Level3Output(l3, output_handler)
            logger.info("Write %s product: %s" % (output_handler.id, output.export_filename))

    def _log_progress(self, i):
        """ Concise logging on the progress of l2i stack creation """
        n = len(self._l2i_files)
        progress_percent = float(i + 1) / float(n) * 100.
        current_reminder = np.mod(progress_percent, 10)
        last_reminder = np.mod(self._l3_progress_percent, 10)
        if last_reminder > current_reminder:
            logger.info(
                'Creating l2i orbit stack: %3g%% (%g of %g)'
                % (progress_percent - current_reminder, i + 1, n)
            )

        self._l3_progress_percent = progress_percent

    @staticmethod
    def apply_orbit_filter(l2i, orbit_filter):
        """
        Apply a
        :param l2i:
        :param orbit_filter:
        :return:
        """

        # Display warning if filter is active
        logger.warning("Orbit filter is active [%s]" % str(orbit_filter.mask_orbits))

        # Get indices to filter
        if orbit_filter.mask_orbits == "ascending":
            indices = np.where(np.ediff1d(l2i.latitude) > 0.)[0]
        elif orbit_filter.mask_orbits == "descending":
            indices = np.where(np.ediff1d(l2i.latitude) < 0.)[0]
        else:
            logger.error(
                "Invalid orbit filter target, needs to be [ascending, descending], Skipping filter ...")
            indices = []

        # Filter geophysical parameters only
        targets = l2i.parameter_list
        for non_target in ["longitude", "latitude", "timestamp", "time", "surface_type"]:
            try:
                targets.remove(non_target)
            except ValueError:
                pass
        l2i.mask_variables(indices, targets)

    @staticmethod
    def apply_miz_filter(l2i, miz_filter):
        """
        Flag values based on the miz filter value
        :param l2i:
        :param miz_filter:
        :return:
        """

        flag_miz = getattr(l2i, "flag_miz", None)
        if flag_miz is None:
            return

        idx = np.where(flag_miz >= miz_filter["mask_min_value"])[0]
        l2i.mask_variables(idx, miz_filter["mask_targets"])

    def _apply_processing_items(self, l3grid):
        """
        Sequentially apply the processing items defined in the Level-3 processor definition files
        (listed under root.processing_items)
        :param l3grid:
        :return:
        """

        # Get the post processing options
        processing_items = self._job.l3def.get("processing_items", None)
        if processing_items is None:
            logger.info("No processing items defined")
            return

        # Get the list of post-processing items
        for pitem in processing_items:
            msg = "Apply Level-3 processing item: `%s`" % (pitem["label"])
            logger.info(msg)
            pp_class = get_cls(pitem["module_name"], pitem["class_name"], relaxed=False)
            processing_items = pp_class(l3grid, **pitem["options"])
            processing_items.apply()


# %% Data Containers

class L2iDataStack(DefaultLoggingClass):

    def __init__(self, griddef, l2_parameter):
        """ A container for stacking l2i variables (geophysical parameter at sensor resolution) in L3 grid cells.
        For each parameters a (numx, numy) array is created, with an list containing all l2i data points that
        fall into the grid cell area. This list can be averaged or used for histogram computation in later stages
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
        self.timeliness = []

        # Save global attributes from l2i (will be overwritten for each
        # l2i file, but it is assumed that general information, e.g.
        # on auxdata remains the same)
        self._l2i_info = None

        # Create parameter stacks
        self._initialize_stacks()

    def _initialize_stacks(self):
        """ Create all data stacks, content will be added sequentially
        with `add` method """

        # Stack dictionary that will hold the data
        self.stack = {
            parameter_name: self.parameter_stack
            for parameter_name in self.l2_parameter.keys()
        }

    def add(self, l2i):
        """ Add a l2i data object to the stack

        Args:
            l2i (obj): l2i object (currently: pysiral.l2data.L2iNCFileImport)

        Returns:
            None
        """

        # Save the metadata from the orbit data
        time = l2i.time if hasattr(l2i, "time") else l2i.timestamp
        self.start_time.append(time[0])
        self.stop_time.append(time[-1])
        self.mission.append(l2i.mission)
        self.timeliness.append(l2i.timeliness)
        self._l2i_count += 1
        self._n_records += l2i.n_records

        self._l2i_info = l2i.info

        # Get projection coordinates for l2i locations
        xi, yj = self.griddef.grid_indices(l2i.longitude, l2i.latitude)

        # Stack the l2 parameter in the corresponding grid cells
        outside_grid_flag = False
        for i in np.arange(l2i.n_records):

            # Add the surface type per default
            # (will not be gridded, therefore not in list of l2 parameter)
            x, y = int(xi[i]), int(yj[i])
            for parameter_name in self.l2_parameter.keys():
                try:
                    data = getattr(l2i, parameter_name)
                    self.stack[parameter_name][y][x].append(data[i])
                except AttributeError:
                    pass
                except IndexError:
                    outside_grid_flag = True

        if outside_grid_flag:
            logger.warning("L2 input data outside grid definition")

    @property
    def n_total_records(self):
        return int(self._n_records)

    @property
    def l2i_count(self):
        return int(self._l2i_count)

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

    def __init__(self, job, stack, period, doi=""):

        super(L3DataGrid, self).__init__(self.__class__.__name__)

        self.error = ErrorStatus(caller_id=self.__class__.__name__)

        # Grid size definition
        self._doi = doi
        self._data_record_type = "none"
        self._griddef = job.grid
        self._l3def = job.l3def
        self._period = period
        self._external_masks = {}

        # Check if any output definition requires an unlimited dimension
        # XXX: This is an unsatisfying solution resulting from the current solution that the L3DataGrid
        #      provides fix dimensions for potentially more than one output handler. This need
        #      consideration in the future
        self._time_dim_is_unlimited = [handler.time_dim_is_unlimited for handler in job.outputs]

        # Define time of dataset creation as the time of object initialization
        # to avoid slightly different timestamps for repeated calls of
        # datetime.now()
        self._creation_time = datetime.now()

        # # Shortcut to the surface type flag dictionary
        # self._surface_type_dict = SurfaceType.SURFACE_TYPE_DICT

        # List of level-2 parameter
        # (gridded parameter that are already in l2i)
        self._l2_parameter = None

        # list of stacked l2 parameters for each grid cell
        if not isinstance(stack, L2iDataStack):
            msg = "Input must be of type pysiral.l3proc.L2DataStack, was %s"
            msg %= type(stack)
            raise ValueError(msg)
        self.l2 = stack

        # Get a list of non-empty
        self._non_empty_grid_indices = None
        self._init_grid_indices_mask()

        # container for gridded parameters
        self.vars = {}

        # Product Metadata
        self._metadata = None
        self._init_metadata_from_l2()

        # Compute the longitude & latitude dimensions
        self.calculate_longitude_latitude_fields()

        # Create the parameter fields
        self._init_parameter_fields(job.l2_parameter)

        # Grid the Level-2 parameter
        logger.info("Grid Level-2 parameter")
        self.grid_l2_parameter()

    def set_doi(self, doi):
        # TODO: Move to __init__
        self._doi = doi

    def set_data_record_type(self, data_record_type):
        self._data_record_type = data_record_type

    def get_attribute(self, attribute_name, *args):
        """ Return a string for a given attribute name. This method is
        required for the output data handler """
        try:
            attr_getter = getattr(self, "_get_attr_" + attribute_name)
            return attr_getter(*args)
        except AttributeError:
            return "attr_unavailable"
        except Exception as ex:
            msg = "L3DataGrid.get_attribute Exception: " + str(ex) + " for attribute: %s" % attribute_name
            logger.error(msg)
            return "unknown"

    def add_grid_variable(self, parameter_name, fill_value, dtype):
        """
        Add a grid variable and fill with empty values
        :param parameter_name: The name of the parameter
        :param fill_value: the "empty" value assigned to all cell
        :param dtype: numpy compatible dtype
        :return:
        """

        # Check if variable already exists
        if parameter_name in self.vars:
            msg = "Variable overwrite alert: %s" % parameter_name
            self.error.add_error("l3-variable-overwrite", msg)
            self.error.raise_on_error()

        # All clear, create variable
        self.vars[parameter_name] = np.full(self.grid_shape, fill_value, dtype=dtype)

        # Log
        logger.info("Added grid parameter: %s" % parameter_name)

    def calculate_longitude_latitude_fields(self):
        """ Geographic coordinates from GridDefinition """
        lon, lat = self.griddef.get_grid_coordinates()
        self.vars["longitude"] = lon
        self.vars["latitude"] = lat

    def grid_l2_parameter(self):
        """ Compute averages of all l2i parameter for each grid cell.
        The list of l2i parameter is from the output format definition
        No averages are computed for grid cells that are tagged with
        a land flag. """

        settings = self.l3def.grid_settings

        # Loop over all parameter / grid cells
        for name in self._l2_parameter:

            # Certain parameters in the l2 stack are excluded from gridding
            # (-> indicated by grid_method: none)
            grid_method = self.l3def.l2_parameter[name]["grid_method"]
            if grid_method == "none":
                continue

            logger.info("Gridding parameter: %s [%s]" % (name, grid_method))

            for xi, yj in self.grid_indices:

                data = np.array(self.l2.stack[name][yj][xi])

                # nanmean needs at least 2 valid items
                valid = np.where(np.isfinite(data))[0]
                if len(valid) < settings.minimum_valid_grid_points:
                    continue

                # TODO: Think of a dicts with lambdas to make this more concise
                if grid_method == "average":
                    self.vars[name][yj, xi] = np.nanmean(data)
                elif grid_method == "average_uncertainty":
                    value = np.abs(np.sqrt(1. / np.sum(data[valid])))
                    self.vars[name][yj, xi] = value
                elif grid_method == "unique":
                    self.vars[name][yj, xi] = np.unique(data)
                elif grid_method == "median":
                    self.vars[name][yj, xi] = np.nanmedian(data)
                else:
                    msg = "Invalid grid method (%s) for %s"
                    msg %= (str(grid_method), name)
                    self.error.add_error("invalid-l3def", msg)
                    self.error.raise_on_error()

    def get_parameter_by_name(self, name, raise_on_error=True):
        try:
            parameter = self.vars[name]
        except KeyError:
            parameter = np.full(np.shape(self.vars["longitude"]), np.nan)
            logger.warning("Parameter not available: %s" % name)
        except Exception as ex:
            msg = "L3DataGrid.get_parameter_by_name Exception: " + str(ex)
            logger.error(msg)
            # TODO: use error handler
            if raise_on_error:
                sys.exit(1)
            else:
                parameter = np.full(np.shape(self.vars["longitude"]), np.nan)
                
        return parameter

    def set_parameter_by_name(self, name, var):
        try:
            self.vars[name] = var
        except KeyError:
            logger.warning("Parameter not available: %s" % name)
        except Exception as ex:
            print("L3DataGrid.get_parameter_by_name Exception: " + str(ex))
            sys.exit(1)

    def _init_grid_indices_mask(self) -> None:
        """
        Compute a mask of non-empty grid indices
        :return:
        """

        # Get number of items per stack
        n_records = np.ndarray(shape=self.grid_shape)
        for xi, yj in self.all_grid_indices:
            n_records[yj][xi] = len(self.l2.stack["time"][yj][xi])

        # Get indices
        self._non_empty_grid_indices = np.flip(np.array(np.where(n_records > 0))).T

    def _init_metadata_from_l2(self):
        """
        Gets metadata from Level-2 instance
        :return:
        """
        # Get the metadata information from the L2 stack
        logger.info("Compile metadata")
        self._metadata = L3MetaData()
        self._metadata.get_missions_from_stack(self.l2)
        # Actual data coverage
        self._metadata.get_data_period_from_stack(self.l2)
        # Requested time coverage (might not be the actual coverage)
        self._metadata.get_time_coverage_from_period(self._period)
        try:
            self._metadata.get_auxdata_infos(self.l2.l2i_info)
        except AttributeError:
            pass
        self._metadata.get_projection_parameter(self._griddef)

    def _init_parameter_fields(self, pardefs):
        """ Initialize output parameter fields """
        # Store the name of the parameters
        self._l2_parameter = sorted(pardefs.keys())
        for parameter_name, pardef in list(pardefs.items()):
            fillvalue = pardef["fillvalue"]
            if pardef["grid_method"] != "none":
                self.add_grid_variable(parameter_name, fillvalue, pardef["dtype"])

    def _get_attr_source_mission_id(self, *args):
        mission_ids = self.metadata.mission_ids
        if args[0] == "uppercase":
            mission_ids = mission_ids.upper()
        return mission_ids

    def _get_attr_source_mission_name(self, *args):
        ids = self.metadata.mission_ids
        return ",".join([psrlcfg.platforms.get_name(m) for m in ids.split(",")])

    def _get_attr_source_timeliness(self, *args):
        timeliness = self.metadata.source_timeliness
        if args[0] == "lowercase":
            timeliness = timeliness.lower()
        elif args[0] == "select":
            choices = {"nrt": args[1], "rep": args[2], "ntc": args[2]}
            return choices.get(timeliness.lower(), "n/a")
        return timeliness

    def _get_attr_grid_id(self, *args):
        grid_id = self.griddef.grid_id
        if args[0] == "uppercase":
            grid_id = grid_id.upper()
        return grid_id

    def _get_attr_grid_spacing_tag(self, *args):
        value = self.griddef.resolution_tag
        if args[0] == "uppercase":
            value = value.upper()
        return value

    def _get_attr_source_mission_sensor(self, *args):
        mission_sensor = self.metadata.mission_sensor
        if args[0] == "uppercase":
            mission_sensor = mission_sensor.upper()
        elif args[0] == "lower":
            mission_sensor = mission_sensor.lower()
        return mission_sensor

    def _get_attr_source_mission_sensor_fn(self, *args):
        """ Same as source mission sensor, only a sanitized version for filenames """
        mission_sensor = self.metadata.mission_sensor
        for character in ["-"]:
            mission_sensor = mission_sensor.replace(character, "")
        if args[0] == "uppercase":
            mission_sensor = mission_sensor.upper()
        return mission_sensor

    def _get_attr_source_hemisphere(self, *args):
        if args[0] != "select":
            return self.hemisphere
        choices = {"north": args[1], "south": args[2]}
        return choices.get(self.hemisphere, "n/a")

    @staticmethod
    def _get_attr_uuid(*args):
        return str(uuid.uuid4())

    def _get_attr_startdt(self, dtfmt):
        return self.metadata.start_period.strftime(dtfmt)

    def _get_attr_stopdt(self, dtfmt):
        return self.info.stop_time.strftime(dtfmt)

    def _get_attr_geospatial_lat_min(self, *args):
        latitude = self.vars["latitude"]
        return self._get_attr_geospatial_str(np.nanmin(latitude))

    def _get_attr_geospatial_lat_max(self, *args):
        latitude = self.vars["latitude"]
        return self._get_attr_geospatial_str(np.nanmax(latitude))

    def _get_attr_geospatial_lon_min(self, *args):
        longitude = self.vars["longitude"]
        return self._get_attr_geospatial_str(np.nanmin(longitude))

    def _get_attr_geospatial_lon_max(self, *args):
        longitude = self.vars["longitude"]
        return self._get_attr_geospatial_str(np.nanmax(longitude))

    @staticmethod
    def _get_attr_geospatial_str(value):
        return "%.4f" % value

    def _get_attr_source_auxdata_sic(self, *args):
        return self.metadata.source_auxdata_sic

    def _get_attr_source_auxdata_snow(self, *args):
        return self.metadata.source_auxdata_snow

    def _get_attr_source_auxdata_sitype(self, *args):
        return self.metadata.source_auxdata_sitype

    def _get_attr_utcnow(self, *args):
        dt = self._creation_time
        return dt.strftime(args[0]) if re.match("%", args[0]) else dt.isoformat()

    def _get_attr_time_coverage_start(self, *args):
        dt = self.metadata.time_coverage_start
        return dt.strftime(args[0]) if re.match("%", args[0]) else dt.isoformat()

    def _get_attr_time_coverage_end(self, *args):
        dt = self.metadata.time_coverage_end
        return dt.strftime(args[0]) if re.match("%", args[0]) else dt.isoformat()

    def _get_attr_time_coverage_duration(self, *args):
        return self.metadata.time_coverage_duration

    def _get_attr_doi(self, *args):
        return self._doi

    def _get_attr_data_record_type(self, *args):
        if args[0] != "select":
            return self._data_record_type
        choices = {"cdr": args[1], "icdr": args[2]}
        return choices.get(self._data_record_type, "n/a")

    @staticmethod
    def _get_attr_pysiral_version(*args):
        return __version__

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
    def all_grid_indices(self):
        return itertools.product(self.grid_xi_range, self.grid_yj_range)

    @property
    def grid_indices(self):
        if self._non_empty_grid_indices is None:
            return itertools.product(self.grid_xi_range, self.grid_yj_range)
        else:
            return self._non_empty_grid_indices

    # @property
    # def parameter_list(self):
    #     # TODO: Only L2 parameter for now
    #     parameter_list = list(self._l2_parameter)
    #     parameter_list.extend(self._l3_parameter)
    #     parameter_list.append("lon")
    #     parameter_list.append("lat")
    #     return parameter_list

    @property
    def dimdict(self):
        time_dim = 0 if True in self._time_dim_is_unlimited else 1
        return OrderedDict([("time", time_dim),
                            ("yc", self.griddef.extent.numx),
                            ("xc", self.griddef.extent.numy)])

    @property
    def grid_shape(self):
        return self.griddef.extent.numx, self.griddef.extent.numy

    @property
    def griddef(self):
        return self._griddef

    @property
    def l3def(self):
        return self._l3def

    @property
    def hemisphere(self):
        return self.metadata.hemisphere

    @property
    def time_bounds(self):
        return [self.metadata.time_coverage_start,
                self.metadata.time_coverage_end]


class L3MetaData(object):
    """
    Container for L3S Metadata information
    (see property attribute_list for a list of attributes)
    """

    _attribute_list = [
        "mission_ids", "start_time", "stop_time", "grid_name", "period_label",
        "time_coverage_start", "time_coverage_end", "time_coverage_duration",
        "time_coverage_resolution", "pysiral_version", "projection_str", "grid_tag", "resolution_tag",
        "hemisphere", "mission_sensor", "source_auxdata_sic",
        "source_auxdata_sitype", "source_auxdata_snow", "source_timeliness"]

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
        mission_sensor = [psrlcfg.platforms.get_sensor(mission.lower()) for mission in missions]

        try:
            self.set_attribute("mission_ids", ",".join(missions))
        except TypeError:
            self.set_attribute("mission_ids", "unkown")
        try:
            self.set_attribute("mission_sensor", ",".join(mission_sensor))
        except TypeError:
            self.set_attribute("mission_ids", "unkown")

        source_timeliness = np.unique(stack.timeliness)[0]
        self.set_attribute("source_timeliness", source_timeliness)

    def get_data_period_from_stack(self, stack):
        """ Get the first and last timestamp """
        self.set_attribute("start_time", np.amin(stack.start_time))
        self.set_attribute("stop_time", np.amax(stack.stop_time))
        # XXX: Only monthly periods are currently supported
        self.set_attribute("period_label", self.start_time.strftime("%B %Y"))

    def get_time_coverage_from_period(self, period):
        """ Get the start and end of requested data period """
        # TODO: replace with the period.get_netcdf_attributes() method
        self.set_attribute("time_coverage_start", period.tcs.dt)
        self.set_attribute("time_coverage_end", period.tce.dt)
        self.set_attribute("time_coverage_duration", period.duration.isoformat)

    def get_auxdata_infos(self, l2i_info):
        """
        Get information on auxiliary data sources from l2i global
        attributes
        TODO: This part is deprecated
        :param l2i_info:
        :return:
        """
        try:
            self.set_attribute("source_auxdata_sic", l2i_info.source_sic)
        except AttributeError:
            self.set_attribute("source_auxdata_sic", l2i_info.source_auxdata_sic)
        try:
            self.set_attribute("source_auxdata_sitype", l2i_info.source_sitype)
        except AttributeError:
            self.set_attribute("source_auxdata_sitype", l2i_info.source_auxdata_sitype)
        try:
            self.set_attribute("source_auxdata_snow", l2i_info.source_snow)
        except AttributeError:
            self.set_attribute("source_auxdata_snow", l2i_info.source_auxdata_snow)

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
        return {field: getattr(self, field) for field in self.attribute_list}

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
                 overwrite_protection=True, period="default", doi=None,
                 data_record_type="none"):

        if output_def == "default":
            output_def = self.default_output_def_filename

        super(Level3OutputHandler, self).__init__(
            output_def, applicable_data_level=3, subfolder_tags=["year"],
            default_file_location=["settings", "outputdef", "l3_default.yaml"])

        self.error.caller_id = self.__class__.__name__
        logger.name = self.__class__.__name__

        self._period = period
        self._doi = doi
        self.overwrite_protection = overwrite_protection

        self._init_product_directory(base_directory)
        self._data_record_type = data_record_type

    def get_filename_from_data(self, l3):
        """ Return the filename for a defined level-3 data object
        based on tag filenaming in output definition file """

        # Get the filenaming definition (depending on period definition)
        filename_template = ""
        try:
            template_ids = self.output_def.filenaming.keys()
            period_id = self._period

            # Add a translation for the current dissonance between dateperiod
            # period id's and the pysiral convention
            # period_id_dict = dict(month="monthly", isoweek="weekly")
            # period_id = period_id_dict.get(period_id, period_id)
            # Fall back to default if no filenaming convention for given
            # data period
            if period_id not in template_ids:
                period_id = "default"
            filename_template = self.output_def.filenaming[period_id]
        except AttributeError:
            filename_template = self.output_def.filenaming
        except KeyError:
            msg = "Missing filenaming convention for period [%s] in [%s]"
            msg %= (str(self._period), self.output_def_filename)
            self.error.add_error("invalid-outputdef", msg)
            self.error.raise_on_error()

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
        return Path(export_directory) / export_filename

    def get_global_attribute_dict(self, l3):
        attr_dict = OrderedDict()
        for attr_entry in self.output_def.global_attributes:
            attr_name, attr_template = zip(*attr_entry.items())
            attribute = self.fill_template_string(attr_template[0], l3)
            attr_dict[attr_name[0]] = attribute
        return attr_dict

    def _init_product_directory(self, base_directory_or_id):
        """ Initializes the product directory. If `base_directory` is already
        a directory, it is used as is. Else, it is treated as subfolder of
        the default pysiral product directory. Product level id and
        overwrite protection subfolders are added for both options"""
        # argument is directory
        if Path(base_directory_or_id).is_dir():
            basedir = Path(base_directory_or_id)
        # argument is id
        else:
            basedir = Path(self.pysiral_config.local_machine.product_repository) / base_directory_or_id
        # add product level subfolder
        # period_id = dict(month="monthly", isoweek="weekly")
        basedir = basedir / self.product_level_subfolder / self._period
        # optional (subfolder with current time)
        if self.overwrite_protection:
            basedir = basedir / self.now_directory
        # set the directory
        self._set_basedir(basedir)

    @property
    def flip_yc(self):
        flip_yc = self.output_def.grid_options.flip_yc
        if not isinstance(flip_yc, bool):
            flip_yc = False
        return flip_yc

    @property
    def doi(self):
        return self._doi

    @property
    def data_record_type(self):
        return self._data_record_type

    @property
    def time_dim_is_unlimited(self):

        # This property has been added. Older L3 output definitions may not have it,
        # -> Catch attribute error and return false if attribute does not exist
        if "time_dim_is_unlimited" not in self.output_def.grid_options:
            msg = "`grid_options.time_dim_is_unlimited` is missing in l3 settings file: %s (Using default: False)"
            logger.warning(msg % self.output_def_filename)
            time_dim_is_unlimited = False
        else:
            time_dim_is_unlimited = self.output_def.grid_options.time_dim_is_unlimited

        # Verification: Value must be bool
        if not isinstance(time_dim_is_unlimited, bool):
            msg = 'Invalid value type for `grid_options.time_dim_is_unlimited` in %s. ' + \
                  'Must be bool, value was %s. (Using default: False)'
            msg = msg % (self.output_def_filename, str(time_dim_is_unlimited))
            logger.error(msg)
            time_dim_is_unlimited = False

        return time_dim_is_unlimited


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
            output (Level-3 compliant output handler from pysiral.core.output)
        """
        super(Level3ProductDefinition, self).__init__(self.__class__.__name__)
        self.error = ErrorStatus(caller_id=self.__class__.__name__)
        self._l3_settings_file = l3_settings_file
        self._output = output
        self._grid = grid
        self._period = period
        self._parse_l3_settings()

        # Report settings to log handler
        logger.info("Output grid id: %s" % str(self._grid.grid_id))
        for output in self._output:
            msg = "L3 product directory (%s): %s"
            msg %= (str(output.id), str(output.basedir))
            logger.info(msg)

    def _parse_l3_settings(self):
        logger.info("Parsing settings: %s" % str(self._l3_settings_file))
        try:
            self._l3 = get_yaml_config(self._l3_settings_file)
        except Exception as ex:
            self.error.add_error("l3settings-parser-error", str(ex))
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
        try:
            mask_names = sorted(self.l3def.l3_masks.keys(branch_mode="only"))
            return [self.l3def.l3_masks[name] for name in mask_names]
        except AttributeError:
            return []

    @property
    def l3_external_masks(self):
        try:
            extmask_names = self.l3def.external_masks
        except AttributeError:
            extmask_names = []
        return extmask_names

    @property
    def l3_post_processors(self):
        try:
            names = sorted(self.l3def.l3_post_processing.keys(
                recursive=False, branch_mode="only"))
        except AttributeError:
            return []
        options = [self.l3def.l3_post_processing[n].options for n in names]
        return zip(names, options)

    @property
    def l2_parameter(self):
        """ Extract a list of paramter names to be extracted from
        l2i product files """
        return self.l3def.l2_parameter

    @property
    def l3_parameter(self):
        """ Extract a list of paramter names to be computed by the
        Level-3 processor """
        l3_parameter = sorted(self.l3def.l3_parameter.keys(branch_mode="only"))
        return [self.l3def.l3_parameter[n] for n in l3_parameter]


class Level3ProcessorItem(DefaultLoggingClass):
    """
    A parent class for processing items to be selected in the Level-3 processor settings
    and applied in the Level3Processor
    """

    def __init__(self, l3grid, **cfg):
        """
        Initizalizes the Level-3 processor item and performs checks if all option input parameters are available.
        :param l3grid: the Level3DataGrid instance to be processed
        :param cfg: The option dictionary/treedict from the config settings file
        """

        super(Level3ProcessorItem, self).__init__(self.__class__.__name__)

        # Add error handler
        self.error = ErrorStatus(caller_id=self.__class__.__name__)

        # Store the arguments with type validation
        if not isinstance(l3grid, L3DataGrid):
            msg = "Invalid data type [%s] for l3grid parameter. Must be l3proc.L3DataGrid"
            msg %= type(l3grid)
            self.error.add_error("invalid-argument", msg)
            self.error.raise_on_error()
        self.l3grid = l3grid
        self.cfg = cfg

        # run the input validation checks
        self._check_variable_dependencies()
        self._check_options()

        # Add empty parameters to the l3grid
        self._add_l3_variables()

    def _check_variable_dependencies(self):
        """
        Tests if the Level-3 data grid has all required input variables (both in the Level-2 stack as
        well as in the Level 3 parameters). All processor item classes that are inheriting this class
        require the properties `l3_variable_dependencies` & `l2_variable_dependencies` for this method to work. Both
        parameter should return a list of variable names. Empty lists should be returned in case of no
        dependency.
        :return:
        """

        # Check Level-2 stack parameter
        for l2_var_name in self.l2_variable_dependencies:
            if l2_var_name not in self.l3grid.l2.stack:
                msg = "Level-3 processor item %s requires l2 stack parameter [%s], which does not exist"
                msg %= (self.__class__.__name__, l2_var_name)
                self.error.add_error("l3procitem-missing-l2stackitem", msg)
                self.error.raise_on_error()

        # Check Level-3 grid parameter
        for l3_var_name in self.l3_variable_dependencies:
            if l3_var_name not in self.l3grid.vars:
                msg = "Level-3 processor item %s requires l3 grid parameter [%s], which does not exist"
                msg %= (self.__class__.__name__, l3_var_name)
                self.error.add_error("l3procitem-missing-l3griditem", msg)
                self.error.raise_on_error()

    def _check_options(self):
        """
        Tests if the all options are given in the Level-3 processor definition files. All processor item
        classes require the property `required_options` (list of option names) for this method to work.
        NOTE: It is in the spirit of pysiral of having all numerical values in one place only that ideally
              is not the code itself.
        :return:
        """
        for option_name in self.required_options:
            option_value = self.cfg.get(option_name, None)
            if option_value is None:
                msg = "Missing option `%s` in Level-3 processor item %s" % (option_name, self.__class__.__name__)
                self.error.add_error("l3procitem-missing-option", msg)
                self.error.raise_on_error()
            setattr(self, option_name, option_value)

    def _add_l3_variables(self):
        """
        This method initializes the output variables for a given processing item to the l3grid. All processor item
        classes require the property `l3_output_variables` for this method to work. The property should return a
        dict with variable names as keys and the value a dict with fill_value and data type.
        :return:
        """
        for variable_name in self.l3_output_variables.keys():
            vardef = self.l3_output_variables[variable_name]
            self.l3grid.add_grid_variable(variable_name, vardef["fill_value"], vardef["dtype"])


class Level3ValidSeaIceFreeboardCount(Level3ProcessorItem):
    """
    A Level-3 processor item to count valid sea ice freeboard values.
    This class should be used for very limited l2i outputs files that
    do not have the surface_type variable necessary for
    `Level3SurfaceTypeStatistics`
    """

    # Mandatory properties
    required_options = []
    l2_variable_dependencies = []
    l3_variable_dependencies = ["sea_ice_freeboard"]
    l3_output_variables = dict(n_valid_freeboards=dict(dtype="i4", fill_value=0))

    def __init__(self, *args, **kwargs):
        """
        Compute surface type statistics
        :param args:
        :param kwargs:
        """
        super(Level3ValidSeaIceFreeboardCount, self).__init__(*args, **kwargs)

    def apply(self):
        """
        Computes the number of valid sea_ice_freeboard values
        """
        for xi, yj in self.l3grid.grid_indices:

            # Extract the list of surface types inm the grid cell
            sea_ice_freeboards = np.array(self.l3grid.l2.stack["sea_ice_freeboard"][yj][xi])
            self.l3grid.vars["n_valid_freeboards"][yj, xi] = np.where(np.isfinite(sea_ice_freeboards))[0].size


class Level3SurfaceTypeStatistics(Level3ProcessorItem):
    """ A Level-3 processor item to compute surface type stastics """

    # Mandatory properties
    required_options = []
    l2_variable_dependencies = ["surface_type", "sea_ice_thickness"]
    l3_variable_dependencies = []
    l3_output_variables = dict(n_total_waveforms=dict(dtype="f4", fill_value=0.0),
                               n_valid_waveforms=dict(dtype="f4", fill_value=0.0),
                               valid_fraction=dict(dtype="f4", fill_value=0.0),
                               lead_fraction=dict(dtype="f4", fill_value=0.0),
                               seaice_fraction=dict(dtype="f4", fill_value=0.0),
                               ocean_fraction=dict(dtype="f4", fill_value=0.0),
                               negative_thickness_fraction=dict(dtype="f4", fill_value=0.0),
                               is_land=dict(dtype="i2", fill_value=-1))

    def __init__(self, *args, **kwargs):
        """
        Compute surface type statistics
        :param args:
        :param kwargs:
        """
        super(Level3SurfaceTypeStatistics, self).__init__(*args, **kwargs)

        # Init this class
        self._surface_type_dict = SURFACE_TYPE_DICT

    def apply(self):
        """
        Computes the mandatory surface type statistics on the surface type stack flag

        The current list
          - is_land (land flag exists in l2i stack)
          - n_total_waveforms (size of l2i stack)
          - n_valid_waveforms (tagged as either lead, sea ice or ocean )
          - valid_fraction (n_valid/n_total)
          - lead_fraction (n_leads/n_valid)
          - ice_fraction (n_ice/n_valid)
          - ocean_fraction (n_ocean/n_valid)
          - negative thickness fraction (n_sit<0 / n_sit)
        """

        # Loop over all grid indices
        stflags = self._surface_type_dict
        for xi, yj in self.l3grid.grid_indices:

            # Extract the list of surface types inm the grid cell
            surface_type = np.array(self.l3grid.l2.stack["surface_type"][yj][xi])

            # Stack can be empty
            if len(surface_type) == 0:
                continue

            # Create a land flag
            is_land = len(np.where(surface_type == stflags["land"])[0] > 0)
            self.l3grid.vars["is_land"][xi, yj] = is_land

            # Compute total waveforms in grid cells
            n_total_waveforms = len(surface_type)
            self.l3grid.vars["n_total_waveforms"][yj, xi] = n_total_waveforms

            # Compute valid waveforms
            # Only positively identified waveforms (either lead or ice)
            valid_waveform = ORCondition()
            valid_waveform.add(surface_type == stflags["lead"])
            valid_waveform.add(surface_type == stflags["sea_ice"])
            valid_waveform.add(surface_type == stflags["ocean"])
            n_valid_waveforms = valid_waveform.num
            self.l3grid.vars["n_valid_waveforms"][yj, xi] = n_valid_waveforms

            # Fractions of leads on valid_waveforms
            try:
                valid_fraction = float(n_valid_waveforms) / float(n_total_waveforms)
            except ZeroDivisionError:
                valid_fraction = 0
            self.l3grid.vars["valid_fraction"][yj, xi] = valid_fraction

            # Fractions of surface types with respect to valid_waveforms
            for surface_type_name in ["ocean", "lead", "sea_ice"]:
                n_wfm = len(np.where(surface_type == stflags[surface_type_name])[0])
                try:
                    detection_fraction = float(n_wfm) / float(n_valid_waveforms)
                except ZeroDivisionError:
                    detection_fraction = 0
                surface_type_id = surface_type_name.replace("_", "")
                self.l3grid.vars[f"{surface_type_id}_fraction"][yj, xi] = detection_fraction

            # Fractions of negative thickness values
            sit = np.array(self.l3grid.l2.stack["sea_ice_thickness"][yj][xi])
            n_negative_thicknesses = len(np.where(sit < 0.0)[0])
            try:
                n_ice = len(np.where(surface_type == stflags["sea_ice"])[0])
                negative_thickness_fraction = float(n_negative_thicknesses) / float(n_ice)
            except ZeroDivisionError:
                negative_thickness_fraction = np.nan
            self.l3grid.vars["negative_thickness_fraction"][yj, xi] = negative_thickness_fraction


class Level3TemporalCoverageStatistics(Level3ProcessorItem):
    """
    A Level-3 processor item to compute temporal coverage statistics of sea-ice thickness in the grid period
    """

    # Mandatory properties
    required_options = []
    l2_variable_dependencies = ["time", "sea_ice_thickness"]
    l3_variable_dependencies = []
    l3_output_variables = dict(temporal_coverage_uniformity_factor=dict(dtype="f4", fill_value=np.nan),
                               temporal_coverage_day_fraction=dict(dtype="f4", fill_value=np.nan),
                               temporal_coverage_period_fraction=dict(dtype="f4", fill_value=np.nan),
                               temporal_coverage_weighted_center=dict(dtype="f4", fill_value=np.nan))

    def __init__(self, *args, **kwargs):
        """
        Compute surface type statistics
        :param args:
        :param kwargs:
        """
        super(Level3TemporalCoverageStatistics, self).__init__(*args, **kwargs)

    def apply(self):
        """
        Computes statistics of the temporal coverage of sea ice thickness
        :return:
        """

        # Other parameter for L3DataGrid
        # All statistics are computed with respect to the temporal coverage of the grid
        # (-> the period that has been asked for, not the actual data coverage)
        # TODO: TCS & TCE as datetime should not be in the metadata, but a property of the data object
        tcs, tce = self.l3grid.metadata.time_coverage_start, self.l3grid.metadata.time_coverage_end
        start_date = date(tcs.year, tcs.month, tcs.day)
        end_date = date(tce.year, tce.month, tce.day)
        period_n_days = (end_date - start_date).days + 1

        # Links
        stack = self.l3grid.l2.stack

        # Loop over all grid cells
        for xi, yj in self.l3grid.grid_indices:

            # Get the day of observation for each entry in the Level-2 stack
            times = np.array(stack["time"][yj][xi])
            day_of_observation = np.array([date(t.year, t.month, t.day) for t in times])

            # The statistic is computed for sea ice thickness -> remove data points without valid sea ice thickness
            sea_ice_thickness = np.array(stack["sea_ice_thickness"][yj][xi])
            day_of_observation = day_of_observation[np.isfinite(sea_ice_thickness)]

            # Validity check
            #  - must have data
            if len(day_of_observation) == 0:
                continue

            # Compute the number of days for each observation with respect to the start of the period
            day_number = [(day - start_date).days for day in day_of_observation]

            # Compute the set of days with observations available
            days_with_observations = np.unique(day_number)
            first_day, last_day = np.amin(days_with_observations), np.amax(days_with_observations)

            # Compute the uniformity factor
            # The uniformity factor is derived from a Kolmogorov-Smirnov (KS) test for goodness of fit that tests
            # the list of against a uniform distribution. The definition of the uniformity factor is that is
            # reaches 1 for uniform distribution of observations and gets smaller for non-uniform distributions
            # It is therefore defined as 1-D with D being the result of KS test
            ks_test_result = stats.kstest(day_number, stats.uniform(loc=0.0, scale=period_n_days).cdf)
            uniformity_factor = 1.0 - ks_test_result[0]
            self.l3grid.vars["temporal_coverage_uniformity_factor"][yj, xi] = uniformity_factor

            # Compute the day fraction (number of days with actual data coverage/days of period)
            day_fraction = float(len(days_with_observations)) / float(period_n_days)
            self.l3grid.vars["temporal_coverage_day_fraction"][yj, xi] = day_fraction

            # Compute the period in days that is covered between the first and last day of observation
            # normed by the length of the period
            period_fraction = float(last_day - first_day + 1) / float(period_n_days)
            self.l3grid.vars["temporal_coverage_period_fraction"][yj, xi] = period_fraction

            # Compute the temporal center of the actual data coverage in units of period length
            # -> optimum 0.5
            weighted_center = np.mean(day_number) / float(period_n_days)
            self.l3grid.vars["temporal_coverage_weighted_center"][yj, xi] = weighted_center


class Level3StatusFlag(Level3ProcessorItem):
    """
    A Level-3 processor item to compute the status flag
    """

    # Mandatory properties
    required_options = ["retrieval_status_target", "sic_thrs", "flag_values"]
    l2_variable_dependencies = []
    l3_variable_dependencies = ["sea_ice_concentration", "n_valid_waveforms", "landsea"]
    l3_output_variables = dict(status_flag=dict(dtype="i1", fill_value=1))

    def __init__(self, *args, **kwargs):
        """
        Compute surface type statistics
        :param args:
        :param kwargs:
        """
        super(Level3StatusFlag, self).__init__(*args, **kwargs)

    def apply(self):
        """
        Computes the status flag
        :return:
        """

        # Get the flag values from the l3 settings file
        flag_values = self.flag_values

        # Get status flag (fill value should be set to zero)
        sf = np.copy(self.l3grid.vars["status_flag"])

        # Init the flag with not data flag value
        sf[:] = flag_values["no_data"]

        # get input parameters
        par = np.copy(self.l3grid.vars[self.retrieval_status_target])
        sic = self.l3grid.vars["sea_ice_concentration"]
        nvw = self.l3grid.vars["n_valid_waveforms"]
        lnd = self.l3grid.vars["landsea"]

        # --- Compute conditions for flags ---

        # Get sea ice mask
        is_below_sic_thrs = np.logical_and(sic >= 0., sic < self.sic_thrs)

        # Get the pole hole information
        mission_ids = self.l3grid.metadata.mission_ids.split(",")
        orbit_inclinations = [psrlcfg.platforms.get_orbit_inclination(mission_id) for mission_id in mission_ids]

        # NOTE: due to varying grid cell size, it is no sufficient to just check which grid cell coordinate
        #       is outside the orbit coverage
        is_pole_hole = np.logical_and(
            np.abs(self.l3grid.vars["latitude"]) > (np.amin(orbit_inclinations)-1.0),
            flag_values["no_data"])

        # Check where the retrieval has failed
        is_land = lnd > 0
        has_data = nvw > 0
        has_retrieval = np.isfinite(par)
        retrieval_failed = np.logical_and(
            np.logical_and(has_data, np.logical_not(is_below_sic_thrs)),
            np.logical_not(has_retrieval))

        # Set sic threshold
        sf[np.where(is_below_sic_thrs)] = flag_values["is_below_sic_thrs"]

        # Set pole hole (Antarctica: Will be overwritten below)
        sf[np.where(is_pole_hole)] = flag_values["is_pole_hole"]

        # Set land mask
        sf[np.where(is_land)] = flag_values["is_land"]

        # Set failed retrieval
        sf[np.where(retrieval_failed)] = flag_values["retrieval_failed"]

        # Set retrieval successful
        sf[np.where(has_retrieval)] = flag_values["has_retrieval"]

        # Write Status flag
        self.l3grid.vars["status_flag"] = sf


class Level3QualityFlag(Level3ProcessorItem):
    """
    A Level-3 processor item to compute the status flag
    """

    # Mandatory properties
    required_options = ["add_rule_flags", "rules"]
    l2_variable_dependencies = []
    l3_variable_dependencies = ["sea_ice_thickness", "n_valid_waveforms", "negative_thickness_fraction",
                                "lead_fraction"]
    l3_output_variables = dict(quality_flag=dict(dtype="i1", fill_value=3))

    def __init__(self, *args, **kwargs):
        """
        Compute surface type statistics
        :param args:
        :param kwargs:
        """
        super(Level3QualityFlag, self).__init__(*args, **kwargs)

    def apply(self):
        """ Computation of quality flag indicator based on several rules defined in the l3 settings file """

        # Get the quality flag indicator array
        # This array will be continously updated by the quality check rules
        qif = np.copy(self.l3grid.vars["quality_flag"])
        sit = np.copy(self.l3grid.vars["sea_ice_thickness"])
        nvw = np.copy(self.l3grid.vars["n_valid_waveforms"])
        ntf = np.copy(self.l3grid.vars["negative_thickness_fraction"])
        lfr = np.copy(self.l3grid.vars["lead_fraction"])

        # As first step set qif to 1 where data is availabe
        qif[np.where(np.isfinite(sit))] = 0

        # Get a list of all the rules
        quality_flag_rules = self.rules.keys()

        # Simple way of handling rules (for now)

        # Use the Warren99 validity maslk
        # XXX: Not implemented yet
        if "qif_warren99_valid_flag" in quality_flag_rules:
            w99 = self.l3grid.vars["warren99_is_valid"]
            # mask = 0 means warren99 is invalid
            rule_options = self.rules["qif_warren99_valid_flag"]
            flag = np.full(qif.shape, 0, dtype=qif.dtype)
            flag[np.where(w99 == 0)] = rule_options["target_flag"]
            qif = np.maximum(qif, flag)

        # Elevate the quality flag for SARin or mixed SAR/SARin regions
        # (only sensible for CryoSat-2)
        if "qif_cs2_radar_mode_is_sin" in quality_flag_rules:
            radar_modes = self.l3grid.vars["radar_mode"]
            rule_options = self.rules["qif_cs2_radar_mode_is_sin"]
            flag = np.full(qif.shape, 0, dtype=qif.dtype)
            flag[np.where(radar_modes >= 2.)] = rule_options["target_flag"]
            qif = np.maximum(qif, flag)

        # Check the number of waveforms (less valid waveforms -> higher warning flag)
        if "qif_n_waveforms" in quality_flag_rules:
            flag = np.full(qif.shape, 0, dtype=qif.dtype)
            rule_options = self.rules["qif_n_waveforms"]
            for threshold, target_flag in zip(rule_options["thresholds"], rule_options["target_flags"]):
                flag[np.where(nvw < threshold)] = target_flag
            qif = np.maximum(qif, flag)

        # Check the availiability of leads in an area adjacent to the grid cell
        if "qif_lead_availability" in quality_flag_rules:
            flag = np.full(qif.shape, 0, dtype=qif.dtype)
            rule_options = self.rules["qif_lead_availability"]
            # get the window size
            grid_res = self.l3grid.griddef.resolution
            window_size = np.ceil(rule_options["search_radius_m"] / grid_res)
            window_size = int(2 * window_size + 1)
            # Use a maximum filter to get best lead fraction in area
            area_lfr = maximum_filter(lfr, size=window_size)
            thrs = rule_options["area_lead_fraction_minimum"]
            flag[np.where(area_lfr <= thrs)] = rule_options["target_flag"]
            qif = np.maximum(qif, flag)

        if "qif_miz_flag" in quality_flag_rules:
            flag = np.full(qif.shape, 0, dtype=qif.dtype)
            rule_options = self.rules["qif_miz_flag"]
            for source_flag, target_flag in zip(rule_options["source_flags"], rule_options["source_flags"]):
                flag[np.where(self.l3grid.vars["flag_miz"] == source_flag)] = target_flag
            qif = np.maximum(qif, flag)

        # Check the negative thickness fraction (higher value -> higher warnung flag)
        if "qif_high_negative_thickness_fraction" in quality_flag_rules:
            flag = np.full(qif.shape, 0, dtype=qif.dtype)
            rule_options = self.rules["qif_high_negative_thickness_fraction"]
            for threshold, target_flag in zip(rule_options["thresholds"], rule_options["target_flags"]):
                flag[np.where(ntf > threshold)] = target_flag
            qif = np.maximum(qif, flag)

        # Set all flags with no data to last flag value again
        qif[np.where(np.isnan(sit))] = 3

        # Set flag again
        self.l3grid.vars["quality_flag"] = qif


class Level3LoadMasks(Level3ProcessorItem):
    """
    A Level-3 processor item to load external masks
    """

    # Mandatory properties
    required_options = ["mask_names"]
    l2_variable_dependencies = []
    l3_variable_dependencies = []
    # Note: the output names depend on mask name, thus these will be
    #       created in apply (works as well)
    l3_output_variables = dict()

    def __init__(self, *args, **kwargs):
        """
        Compute surface type statistics
        :param args:
        :param kwargs:
        """
        super(Level3LoadMasks, self).__init__(*args, **kwargs)

    def apply(self):
        """
        Load masks and add them as grid variable (variable name -> mask name)
        :return:
        """

        # Get the mask names and load each
        for mask_name in self.mask_names:

            # The masks are stored in external files that can be automatically
            #  found with the grid id
            mask = L3Mask(mask_name, self.l3grid.griddef.grid_id)

            # Add the mask to the l3grid as variable
            if not mask.error.status:
                self.l3grid.add_grid_variable(mask_name, np.nan, mask.mask.dtype)
                self.l3grid.vars[mask_name] = mask.mask

            # If fails, only add an empty variable
            else:
                self.l3grid.add_grid_variable(mask_name, np.nan, "f4")
                error_msgs = mask.error.get_all_messages()
                for error_msg in error_msgs:
                    logger.error(error_msg)


class Level3LoadCCILandMask(Level3ProcessorItem):
    """
    A Level-3 processor item to load the CCI land mask
    """

    # Mandatory properties
    required_options = ["local_machine_def_mask_tag", "mask_name_dict"]
    l2_variable_dependencies = []
    l3_variable_dependencies = []
    # Note: the output names depend on mask name, thus these will be
    #       created in apply (works as well)
    l3_output_variables = dict()

    def __init__(self, *args, **kwargs):
        """
        Initiate the class
        :param args:
        :param kwargs:
        """
        super(Level3LoadCCILandMask, self).__init__(*args, **kwargs)

    def apply(self):
        """
        Load masks and add them as grid variable (variable name -> mask name)
        :return:
        """

        # Short cut
        grid_id = self.l3grid.griddef.grid_id

        # Get mask target path:
        mask_tag = self.cfg["local_machine_def_mask_tag"]
        lookup_directory = psrlcfg.local_machine.auxdata_repository.mask.get(mask_tag, None)
        if lookup_directory is None:
            msg = "Missing local machine def tag: auxdata_repository.mask.{}".format(mask_tag)
            self.error.add_error("invalid-local-machine-def", msg)
            logger.error(msg)
            return
        lookup_directory = Path(lookup_directory)

        # Get the mask target filename
        filename = self.cfg["mask_name_dict"][grid_id.replace("_", "")]
        mask_filepath = lookup_directory / filename
        if not mask_filepath.is_file():
            msg = "Missing input file: {}".format(mask_filepath)
            self.error.add_error("invalid-local-machine-def", msg)
            logger.error(msg)
            return

        # Load the data and extract the flag
        nc = open_dataset(str(mask_filepath), decode_times=False)

        # The target land sea flag should be 1 for land and 0 for sea,
        # but the CCI landsea mask provides also fractional values for
        # mixed surfaces types. Thus, we add two arrays to the
        # L3grid object
        #   1. a classical land/sea mask with land:1 and sea: 0. In
        #      this notation the mixed pixels are attributed to sea
        #      because there might be some valid retrieval there
        #   2. the ocean density value as is
        density_of_ocean = np.flipud(nc.density_of_ocean.values)
        landsea_mask = density_of_ocean < 1e-5

        # Add mask to l3 grid
        mask_variable_name = self.cfg["mask_variable_name"]
        self.l3grid.add_grid_variable(mask_variable_name, np.nan, landsea_mask.dtype)
        self.l3grid.vars[mask_variable_name] = landsea_mask.astype(int)

        density_variable_name = self.cfg["density_variable_name"]
        self.l3grid.add_grid_variable(density_variable_name, np.nan, nc.density_of_ocean.values.dtype)
        self.l3grid.vars[density_variable_name] = density_of_ocean


class Level3GridUncertainties(Level3ProcessorItem):
    """
    A Level-3 processor item to compute uncertainties of key geophysical variables on a grid.
    NOTE: As a concession to backward compability: sea ice draft uncertainty will be computed, but
          the sea ice draft is not a required input parameter
    """

    # Mandatory properties
    required_options = ["water_density", "snow_depth_correction_factor", "max_l3_uncertainty"]
    l2_variable_dependencies = ["radar_freeboard_uncertainty", "sea_ice_thickness"]
    l3_variable_dependencies = ["sea_ice_thickness", "sea_ice_freeboard", "snow_depth", "sea_ice_density",
                                "snow_density", "snow_depth_uncertainty", "sea_ice_density_uncertainty",
                                "snow_density_uncertainty"]
    l3_output_variables = dict(radar_freeboard_l3_uncertainty=dict(dtype="f4", fill_value=np.nan),
                               freeboard_l3_uncertainty=dict(dtype="f4", fill_value=np.nan),
                               sea_ice_thickness_l3_uncertainty=dict(dtype="f4", fill_value=np.nan),
                               sea_ice_draft_l3_uncertainty=dict(dtype="f4", fill_value=np.nan))

    def __init__(self, *args, **kwargs):
        """
        Compute surface type statistics
        :param args:
        :param kwargs:
        """
        super(Level3GridUncertainties, self).__init__(*args, **kwargs)

    def apply(self):
        """ Compute a level 3 uncertainty. The general idea is to compute the error propagation of average
        error components, where for components for random error the error of the l2 average
        is used and for systematic error components the average of the l2 error """

        # Options
        rho_w = self.water_density
        sd_corr_fact = self.snow_depth_correction_factor

        # Loop over grid items
        for xi, yj in self.l3grid.grid_indices:

            # Check of data exists
            if np.isnan(self.l3grid.vars["sea_ice_thickness"][yj, xi]):
                continue

            # Get parameters
            frb = self.l3grid.vars["sea_ice_freeboard"][yj, xi]
            sd = self.l3grid.vars["snow_depth"][yj, xi]
            rho_i = self.l3grid.vars["sea_ice_density"][yj, xi]
            rho_s = self.l3grid.vars["snow_density"][yj, xi]

            # Get systematic error components
            sd_unc = self.l3grid.vars["snow_depth_uncertainty"][yj, xi]
            rho_i_unc = self.l3grid.vars["sea_ice_density_uncertainty"][yj, xi]
            rho_s_unc = self.l3grid.vars["snow_density_uncertainty"][yj, xi]

            # Get random uncertainty
            # Note: this applies only to the radar freeboard uncertainty.
            #       Thus we need to recalculate the sea ice freeboard uncertainty

            # Get the stack of radar freeboard uncertainty values and remove NaN's
            # rfrb_unc = self.l3["radar_freeboard_uncertainty"][yj, xi]
            rfrb_uncs = np.array(self.l3grid.l2.stack["radar_freeboard_uncertainty"][yj][xi])
            rfrb_uncs = rfrb_uncs[~np.isnan(rfrb_uncs)]

            # Compute radar freeboard uncertainty as error or the mean from values with individual
            # error components (error of a weighted mean)
            weight = np.nansum(1. / rfrb_uncs ** 2.)
            rfrb_unc = 1. / np.sqrt(weight)
            self.l3grid.vars["radar_freeboard_l3_uncertainty"][yj, xi] = rfrb_unc

            # Calculate the level-3 freeboard uncertainty with updated radar freeboard uncertainty
            deriv_snow = sd_corr_fact
            frb_unc = np.sqrt((deriv_snow * sd_unc) ** 2. + rfrb_unc ** 2.)
            self.l3grid.vars["freeboard_l3_uncertainty"][yj, xi] = frb_unc

            # Calculate the level-3 thickness uncertainty
            errprop_args = [frb, sd, rho_w, rho_i, rho_s, frb_unc, sd_unc, rho_i_unc, rho_s_unc]
            sit_l3_unc = frb2sit_errprop(*errprop_args)

            # Cap the uncertainty
            # (very large values may appear in extreme cases)
            if sit_l3_unc > self.max_l3_uncertainty:
                sit_l3_unc = self.max_l3_uncertainty

            # Assign Level-3 uncertainty
            self.l3grid.vars["sea_ice_thickness_l3_uncertainty"][yj, xi] = sit_l3_unc

            # Compute sea ice draft uncertainty
            if "sea_ice_draft" not in self.l3grid.vars:
                continue

            sid_l3_unc = np.sqrt(sit_l3_unc ** 2. + frb_unc ** 2.)
            self.l3grid.vars["sea_ice_draft_l3_uncertainty"][yj, xi] = sid_l3_unc


class Level3ParameterMask(Level3ProcessorItem):
    """
    A Level-3 processor item to load external masks
    """

    # Mandatory properties
    required_options = ["source", "condition", "targets"]
    l2_variable_dependencies = []
    l3_variable_dependencies = []
    # Note: the output names depend on mask name, thus these will be
    #       created in apply (works as well)
    l3_output_variables = dict()

    def __init__(self, *args, **kwargs):
        """
        Compute surface type statistics
        :param args:
        :param kwargs:
        """
        super(Level3ParameterMask, self).__init__(*args, **kwargs)

    def apply(self):
        """
        Mask certain parameters based on condition of one other parameter
        :return:
        """

        # Get the source parameter
        source = self.l3grid.vars[self.source]

        # Compute the masking condition
        conditions = self.condition.split(";")
        n_conditions = len(conditions)

        if n_conditions == 0:
            msg = "Missing condition in %s" % self.__class__.__name__
            self.error.add_error("invalid-l3mask-def", msg)
            return

        # Start with the first (and maybe only condition)
        filter_mask = self._get_l3_mask(source, conditions[0], self.cfg)

        # Add conditions
        if n_conditions >= 2:
            for i in range(1, n_conditions):
                new_filter = self._get_l3_mask(source, conditions[i], self.cfg)
                if self.cfg["connect_conditions"] == "or":
                    filter_mask = np.logical_or(filter_mask, new_filter)
                elif self.cfg["connect_conditions"] == "and":
                    filter_mask = np.logical_and(filter_mask, new_filter)
                else:
                    msg = "Invalid l3 mask operation: %s"
                    msg %= self.cfg["connect_conditions"]
                    self.error.add_error("invalid-l3mask-def", msg)
                    self.error.raise_on_error()

        # Apply mask
        masked_indices = np.where(filter_mask)
        for target in self.targets:
            try:
                self.l3grid.vars[target][masked_indices] = np.nan
            except ValueError:
                if self.l3grid.vars[target].dtype.kind == "i":
                    self.l3grid.vars[target][masked_indices] = -1
                else:
                    msg = "Cannot set nan (or -1) as mask value to parameter: %s " % target
                    logger.warning(msg)

    def _get_l3_mask(self, source_param, condition, options):
        """ Return bool array based on a parameter and a predefined
        masking operation """
        if condition.strip() == "is_nan":
            return np.isnan(source_param)
        elif condition.strip() == "is_zero":
            return np.array(source_param <= 1.0e-9)
        elif condition.strip() == "is_smaller":
            return np.array(source_param < options["is_smaller_threshold"])
        else:
            msg = "Unknown condition in l3 mask: %s" % condition
            self.error.add_error("invalid-l3mask-condition", msg)
            self.error.raise_on_error()


class Level3GriddedClassifiers(Level3ProcessorItem):
    """
    A Level-3 processor item to provide gridded classifiers (for different surface types)
    """

    # Mandatory properties
    required_options = ["parameters", "surface_types", "statistics"]
    l2_variable_dependencies = ["surface_type"]
    l3_variable_dependencies = []
    # Note: the output names depend on the parameters selected, thus these will be
    #       created in apply (works as well)
    l3_output_variables = dict()

    def __init__(self, *args, **kwargs):
        """
        Compute surface type statistics
        :param args:
        :param kwargs:
        """
        super(Level3GriddedClassifiers, self).__init__(*args, **kwargs)

        # Surface type dict is required to get subsets
        self._surface_type_dict = SURFACE_TYPE_DICT

        # Statistical function dictionary
        self._stat_functions = dict(mean=lambda x: np.nanmean(x), sdev=lambda x: np.nanstd(x))

    def apply(self):
        """
        Mask certain parameters based on condition of one other parameter
        :return:
        """

        # Get surface type flag
        surface_type = self.l3grid.l2.stack["surface_type"]
        target_surface_types = list(self.surface_types)
        target_surface_types.append("all")

        # Loop over all parameters
        for parameter_name in self.parameters:
            # Get the stack
            classifier_stack = None
            try:
                classifier_stack = self.l3grid.l2.stack[parameter_name]
            except KeyError:
                msg = "Level-3 processor item %s requires l2 stack parameter [%s], which does not exist"
                msg %= (self.__class__.__name__, parameter_name)
                self.error.add_error("l3procitem-missing-l2stackitem", msg)
                self.error.raise_on_error()

            # Loop over the statistical parameters (mean, sdev, ...)
            for statistic in self.statistics:
                # Loop over target surface types
                for target_surface_type in target_surface_types:
                    self._compute_grid_variable(parameter_name,
                                                classifier_stack,
                                                surface_type,
                                                target_surface_type,
                                                statistic)

    def _compute_grid_variable(self, parameter_name, classifier_stack, surface_type, target_surface_type, statistic):
        """
        Computes gridded surface type statistics for all grid cells
        :param parameter_name: The name of the classifier (for output name generation)
        :param classifier_stack: The Level-2 stack for the given classifier
        :param surface_type: The Level-2 stack of surface type
        :param target_surface_type: The name of the target surface type
        :param statistic: The name of the statistic to be computed
        :return:
        """

        # Create the output parameter name
        grid_var_name = "stat_%s_%s_%s" % (parameter_name, target_surface_type, statistic)
        self.l3grid.add_grid_variable(grid_var_name, np.nan, "f4")

        # Loop over all grid cells
        for xi, yj in self.l3grid.grid_indices:

            classifier_grid_values = np.array(classifier_stack[yj][xi])
            surface_type_flags = np.array(surface_type[yj][xi])

            # Get the surface type target subset
            if target_surface_type == "all":
                subset = np.arange(len(classifier_grid_values))
            else:
                try:
                    surface_type_target_flag = self._surface_type_dict[target_surface_type]
                    subset = np.where(surface_type_flags == surface_type_target_flag)[0]
                except KeyError:
                    msg = "Surface type %s does not exist" % target_surface_type
                    self.error.add_error("l3procitem-incorrect-option", msg)
                    self.error.raise_on_error()

            # A minimum of two values is needed to compute statistics
            if len(subset) < 2:
                continue
            result = self._stat_functions[statistic](classifier_grid_values[subset])
            self.l3grid.vars[grid_var_name][yj][xi] = result
