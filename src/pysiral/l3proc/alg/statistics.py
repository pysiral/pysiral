# -*- coding: utf-8 -*-

"""

"""

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"

import numpy as np

from datetime import date
from scipy import stats

from pysiral.core.flags import SURFACE_TYPE_DICT, ORCondition
from pysiral.l3proc import Level3ProcessorItem


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
    l2_variable_dependencies = ["surface_type"]
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
            try:
                sit = np.array(self.l3grid.l2.stack["sea_ice_thickness"][yj][xi])
            except KeyError:
                self.l3grid.vars["negative_thickness_fraction"][yj, xi] = np.nan
                continue

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
            if len(classifier_grid_values) == 0:
                continue

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
