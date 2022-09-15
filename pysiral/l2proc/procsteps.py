# -*- coding: utf-8 -*-
"""
@author: Stefan Hendricks
"""

import numpy as np
import pandas as pd
from loguru import logger
from itertools import product

from pysiral import get_cls
from pysiral.l2data import Level2Data
from pysiral.l1bdata import Level1bData
from pysiral._class_template import DefaultLoggingClass


class Level2ProcessorStep(DefaultLoggingClass):
    """
    Parent class for any Level-2 processor step class, which may be distributed over the
    different pysiral modules.

    This class also serves as a template for all subclasses. Mandatory methods and properties
    in this class which raise a NotImplementedException must be overwritten by the subclass
    """

    def __init__(self, cfg):
        """
        Init the
        :param cfg:
        """
        super(Level2ProcessorStep, self).__init__(self.__class__.__name__)

        # -- Properties --

        # Class configuration
        self.cfg = cfg

        # Log messages
        self.msgs = []

        # Error flag dict {code: component}
        self.error_flag_bit_dict = {
            "l1b": 0,
            "l2proc": 1,
            "auxdata": 2,
            "surface_type": 3,
            "retracker": 4,
            "range_correction": 5,
            "sla": 6,
            "frb": 7,
            "sit": 8,
            "filter": 9,
            "other": 15}

    def execute(self, l1b: Level1bData, l2: Level2Data) -> None:
        """
        The main entry point for the
        :param l1b:
        :param l2:
        :return:
        """

        # Execute the method of the subclass. The class needs to
        error_status = self.execute_procstep(l1b, l2)

        # Check if the error status is correctly implemented
        if error_status is None:
            logger.warning(f"Class {self.classname} does not provide error status")
            error_status = self.get_clean_error_status(l2.n_records)

        # Update the status flag
        self.update_error_flag(l2, error_status)

    def update_error_flag(self, l2, error_status):
        """
        Add the error_flag of the the processing step to the
        :param: l2: The Level-2 data container
        :param: error_status: An array with the shape of l2.records containing the error flag
            (False: nominal, True: error)
        :return:
        """
        # Update the algorithm error flag in the Level-2 data object
        # NOTE: The bit for which the algorithm error flag is applied depends on the
        #       on the mandatory value set in each child class of Level2ProcessorStep
        flag_value = 2**self.error_bit
        flag = error_status.astype(int)*flag_value
        l2.flag = l2.flag + flag

    def execute_procstep(self, l1b, l2):
        raise NotImplementedError("This method needs to implemented in {}".format(self.classname))

    @staticmethod
    def get_clean_error_status(shape):
        """
        Return an empty (all nominal error status array)
        :param shape:
        :return:
        """
        return np.full(shape, False)

    @property
    def classname(self):
        return self.__class__.__name__

    @property
    def error_bit(self):
        raise NotImplementedError("This method needs to implemented in {}".format(self.classname))

    @property
    def l2_input_vars(self):
        raise NotImplementedError("This method needs to implemented in {}".format(self.classname))

    @property
    def l2_output_vars(self):
        raise NotImplementedError("This method needs to implemented in {}".format(self.classname))


class Level2ProcessorStepOrder(DefaultLoggingClass):
    """
    A container providing the ordered list of processing steps
    as initialized classes for each trajectory
    """

    def __init__(self, cfg):
        """
        Initialize this class
        :param cfg: the procsteps tag from the Level-2 processor definitions file
        """
        super(Level2ProcessorStepOrder, self).__init__(self.__class__.__name__)

        # Properties
        self.cfg = cfg

        # A list of the class object (not initialized!)
        self._classes = []
        self.get_procstep_classes()

    def get_procstep_classes(self):
        """
        Retrieves the required classes from the processor definition files and stores them in a list
        without initializing them. This way a freshly initialized version can be supplied to each
        l2 data object without risk of interference of class properties
        :return:
        """

        # Loop
        for procstep_def in self.cfg:

            # Get the module & pyclass
            module = procstep_def["module"]
            full_module_name = "pysiral.{}".format(module)
            obj = get_cls(full_module_name, procstep_def["pyclass"])

            # This object should not be None
            if obj is None:
                msg = "Invalid L2 processing step class: {}.{}".format(full_module_name, procstep_def["pyclass"])
                self.error.add_error("missing-class", msg)
                self.error.raise_on_error()

            # Append the class
            logger.info("Added L2 processor step: {}.{}".format(full_module_name, procstep_def["pyclass"]))
            self._classes.append(obj)

    def get_algorithm_error_flag_bit(self, procstep_module):
        """
        Return the corresponding bit (0-15) for the procstep module
        :param procstep_module: str
        :return:
        """

    def validate(self):
        """
        Checkout the difference processing steps and validate input/output variables in
        the order of the steps
        :return:
        """

        # Get the initialized classes and check if it has the required methods
        classes = self.class_instances
        for obj in classes:
            for mandatory_attr in ["execute_procstep", "l2_input_vars", "l2_output_vars", "error_bit"]:
                try:
                    has_attr = hasattr(obj, mandatory_attr)
                except NotImplementedError:
                    has_attr = False
                if not has_attr:
                    msg = "Class {} is missing the mandatory method/property: {}"
                    msg = msg.format(obj.__class__.__name__, mandatory_attr)
                    self.error.add_error("invalid-class", msg)
        self.error.raise_on_error()

    @property
    def class_list(self):
        """
        Return a list of not initialized classes (class instances with the options from the level-2
        :return:
        """
        return list(self._classes)

    @property
    def class_instances(self):
        """
        Return a list of initialized classes (class instances with the options from the level-2 processor
        definition file passed to the class instance)
        :return:
        """
        return [pyclass(opt) for pyclass, opt in zip(self.class_list, self.cfg)]


class L1BL2TransferVariables(Level2ProcessorStep):
    """
    Level-2 Processor step class to transfer variables from the
    l1b to the l2 data object
    """

    def __init__(self, *args, **kwargs):
        """
        Initiliaze the class
        :param arg:
        :param kwargs:
        """
        super(L1BL2TransferVariables, self).__init__(*args, **kwargs)

    def execute_procstep(self, l1b, l2):
        """
        Mandatory method of Level-2 processor
        :param l1b:
        :param l2:
        :return: error_status
        """

        # Get the error mandatory
        error_status = self.get_clean_error_status(l2.n_records)

        logger.info("- Transfer L1P variables to L2")
        for data_group, varlist in list(self.cfg.options.items()):

            # Get and loop over variables per data group
            l1p_variables = varlist.items()
            for var_name, vardef in list(l1p_variables):

                # Get variable via standard getter method
                # NOTE: Will return None if not found -> create an empty array
                #       -> this will be noted in the error status flag
                var = l1b.get_parameter_by_name(data_group, var_name)
                if var is None:
                    logger.warning(f"Cannot find variable {data_group}.{var_name} in l1p")
                    var = np.full(l2.n_records, np.nan)
                    error_status = np.logical_not(error_status)

                # Add variable to l2 object as auxiliary variable
                aux_id, aux_name = vardef
                l2.set_auxiliary_parameter(aux_id, aux_name, var, None)

        return error_status

    @property
    def l2_input_vars(self):
        return []

    @property
    def l2_output_vars(self):
        output_vars = []
        for data_group, varlist in list(self.cfg.options.items()):
            l1p_variables = varlist.items()
            for var_name, vardef in list(l1p_variables):
                output_vars.append(vardef[0])
        return output_vars

    @property
    def error_bit(self):
        return self.error_flag_bit_dict["l2proc"]


class L2ApplyRangeCorrections(Level2ProcessorStep):
    """
    Level-2 processing step class for applying geophysical range corrections to the elevations
    """

    def __init__(self, *args, **kwargs):
        """
        Init the class
        :param args:
        :param kwargs:
        """
        super(L2ApplyRangeCorrections, self).__init__(*args, **kwargs)

    def execute_procstep(self, l1b, l2):
        """
        Mandatory method of Level-2 processor
        :param l1b:
        :param l2:
        :return: error_status
        """

        # Get the error mandatory
        error_status = self.get_clean_error_status(l2.n_records)

        # Keep the total range correction
        total_range_correction = np.full(l2.n_records, 0.0)

        preserve_nan = self.cfg.options.get("preserve_nan", False)

        # Apply the range corrections (content of l1b data package)
        # to the l2 elevation (output of retracker) data
        for correction_name in self.cfg.options.corrections:

            # Get the range correction field
            range_delta = l1b.correction.get_parameter_by_name(correction_name)

            # Check if data is in l1p
            # -> skip with warning if not
            if range_delta is None:
                error_status[:] = True
                logger.warning("Cannot find range correction: {} - skipping".format(correction_name))
                continue

            # Check if NaN's, set values to zero and provide warning
            nans_indices = np.where(np.isnan(range_delta))[0]
            if len(nans_indices) > 0:
                if not preserve_nan:
                    logger.warning(
                        "NaNs encountered and replaced with zero in range correction parameter: %s" % correction_name)
                    range_delta[nans_indices] = 0.0
                else:
                    logger.warning(
                        "NaNs encountered and preserved in range correction parameter: %s" % correction_name)

                error_status[nans_indices] = True

            # Apply the range correction
            #
            # NOTES:
            #
            #    1. It is assumed at this point that the uncertainty of the range correction is
            #       negligible compared to the remaining uncertainty components, respectively
            #       is already included in the general range uncertainty budget
            #
            #    2. the range correction variable `range_delta` is subtracted from the elevation
            #       as the it has to be added to range:
            #
            #           elevation = altitude - (range + range_correction)
            #       ->  elevation = altitude - range - range_correction
            #                       |              |
            #                       ----------------
            #                 elevation after retracking
            #
            for target_variable in self.target_variables:
                var = l2.get_parameter_by_name(target_variable)
                var[:] = var[:] - range_delta
                l2.set_parameter(target_variable, var[:], var.uncertainty[:])

            # store the range correction in the total range correction array
            total_range_correction += range_delta

        # Add the total range correction to the l2 auxiliary variables
        l2.set_auxiliary_parameter("rctotal", "total_range_correction", total_range_correction)

        return error_status

    @property
    def target_variables(self):
        """
        Return target variables from config file or default value
        :return:
        """
        return self.cfg.options.get("target_variables", ["elevation"])

    @property
    def l2_input_vars(self):
        return self.target_variables

    @property
    def l2_output_vars(self):
        return ["rctotal"]

    @property
    def error_bit(self):
        return self.error_flag_bit_dict["range_correction"]


class CS2InstrumentModeflag(Level2ProcessorStep):
    """
    A class creating an instrument_mode flag from radar_mode:

        instrument_mode = radar_mode + 1
    """

    def __init__(self, *args, **kwargs):
        """
        Init the class
        :param args:
        :param kwargs:
        """
        super(CS2InstrumentModeflag, self).__init__(*args, **kwargs)

    def execute_procstep(self, l1b, l2):
        """
        Mandatory method of Level-2 processor
        :param l1b:
        :param l2:
        :return: error_status
        """

        # Get the error mandatory
        error_status = self.get_clean_error_status(l2.n_records)

        # Instrument mode is a variant of radar mode where the flag starts at 1 and not 0
        instrument_mode = l2.radar_mode + 1

        # Add the total range correction to the l2 auxiliary variables
        l2.set_auxiliary_parameter("imode", "instrument_mode", instrument_mode)

        return error_status

    @property
    def l2_input_vars(self):
        return ["radar_mode"]

    @property
    def l2_output_vars(self):
        return ["instrument_mode"]

    @property
    def error_bit(self):
        return self.error_flag_bit_dict["other"]


class ParameterRollingStatistics(Level2ProcessorStep):
    """
    This class add rolling statistics (mean, sdev) using the pandas rolling
    framework to the Level-2 data object.

        -   module: l2proc.procsteps
            pyclass: ParameterRollingStatistics
            options:
                window_size_m: <maximum size of the marginal ice zone>
                statistics: [ <list of statistics: "mean" | "sdev"]
                input_parameters:
                    - <l2_parameter-01>
                    - <l2_parameter-02>
                    - ...

    For each parameter & statistics combination a new parameter will be added
    with the naming convention:

                <source_parameter_name>_rolling_<statistics_id>

        e.g.

        -   module: l2proc.procsteps
            pyclass: ParameterRollingStatistics
            options:
                window_size_m: 25000.
                statistics: ["sdev"]
                input_parameters:
                    - pulse_peakiness

        will add the parameter `pulse_peakiness_rolling_sdev` to the list of l2 (auxiliary)
        data parameters.
    """

    def __init__(self, *args, **kwargs):
        super(ParameterRollingStatistics, self).__init__(*args, **kwargs)

    def execute_procstep(self,
                         l1: "Level1bData",
                         l2: "Level2Data"
                         ) -> np.ndarray:
        """
        API method for Level2ProcessorStep subclasses. Computes and add rolling statistics
        based on the processor item options in the Level-2 processor definition file.

        :param l1: Level-1 data object
        :param l2: Level-2 data object
        
        :raises None:

        :return: error_status
        """

        error_status = self.get_clean_error_status(l2.n_records)

        # The filter windows size needs to be an odd integer
        window_size_float = self.cfg.options.get("window_size_m") / l2.footprint_spacing
        window_size = None if np.isnan(window_size_float) else int(int(window_size_float) // 2 * 2 + 1)

        rolling_kwargs = dict(window=window_size, center=True, min_periods=1)
        for parameter_name, statistics_id in self.statistics_combinations:

            aux_id = f"{l2.full_variable_catalog[parameter_name]}rl"
            aux_name = f"{parameter_name}_rolling_{statistics_id}"

            if window_size is None:
                error_status[:] = True
                l2.set_auxiliary_parameter(aux_id, aux_name, np.full(l2.n_records, np.nan), None)
                continue

            # The rolling statistics are computed using the pandas rolling framework
            ts = pd.Series(l2.get_parameter_by_name(parameter_name)[:])
            if statistics_id == "sdev":
                ts_rolled = ts.rolling(**rolling_kwargs).std()
            elif statistics_id == "mean":
                ts_rolled = ts.rolling(**rolling_kwargs).mean()
            else:
                logger.error(f"{self.__class__.__name__}: Encountered invalid statistic name: {statistics_id}")
                error_status[:] = True
                ts_rolled = pd.Series(np.full(l2.n_records, np.nan))

            # Output is added as an auxiliary parameter
            l2.set_auxiliary_parameter(aux_id, aux_name, ts_rolled.values, None)

        return error_status

    @property
    def statistics_combinations(self):
        return product(self.cfg.options.input_parameters, self.cfg.options.statistics)

    @property
    def l2_input_vars(self):
        return [self.cfg.options.input_parameters]

    @property
    def l2_output_vars(self):
        return [f"{varname}_rolling_{statname}" for (varname, statname) in self.statistics_combinations]

    @property
    def error_bit(self):
        return self.error_flag_bit_dict["other"]
