# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 11:25:04 2015

@author: Stefan
"""


from copy import deepcopy
import numpy as np
from attrdict import AttrDict

from pysiral.core.flags import SurfaceType
from pysiral.l2proc.procsteps import Level2ProcessorStep
from pysiral.config import RadarModes
from pysiral.core.flags import ANDCondition


class ClassifierContainer(object):

    def __init__(self):
        self._parameters = {}

    def add_parameter(self, parameter, parameter_name):
        if self.n_parameters != 0:
            if parameter.shape != self.shape:
                raise ValueError(f"Invalid dimension: {parameter.shape} for {parameter_name} [{self.shape}]")
        self._parameters[parameter_name] = parameter

    def __getattr__(self, item):
        """
        Modify the attribute getter to provide a shortcut to the data content
        :param item: Name of the parameter
        :return:
        """
        if item in list(self._parameters.keys()):
            return self._parameters[item]
        else:
            raise AttributeError()

    @property
    def n_parameters(self):
        return len(self.parameter_list)

    @property
    def parameter_list(self):
        return list(self._parameters.keys())

    @property
    def shape(self):
        if self.n_parameters == 0:
            return ()
        shape = list(set([p.shape for p in self._parameters.values()]))
        return shape[0]


class SurfaceTypeClassifier(object):
    """
    This is a parent class that allows to quickly generate surface type classication schemes
    for the Level-2 processpr
    """

    def __init__(self, *args, **kwargs):
        """
        Initializes the class
        """

        # The result of the surface type classification
        # NOTE: This instance contains the number codes for each surface type
        self.surface_type = SurfaceType()
        
        # A data container for parameters (usually waveform parameters, but also auxiliary data)
        # used to classify the surface type 
        self.classifier = ClassifierContainer()

        # This instance can be used to related the radar mode flag with the name for the radar mode
        # (which will certainly be used in the processor definition file
        self.radar_modes = RadarModes()

        # This list will need to be set by the child class
        # It should contain a list of surface types that the algorithm will detect
        # (see SurfaceType.SURFACE_TYPE_DICT for the names)
        self._classes = []

    def set_initial_classification(self, surface_type):
        """
        This method sets an initial surface type classification
        :param surface_type:
        :return:
        """
        self.surface_type = surface_type

    def transfer_l1b_classifier(self, l1b):
        """
        A standard functionality to transfer all l1b classifier
        :param l1b:
        :return:
        """
        for classifier_name in l1b.classifier.parameter_list:
            classifier = deepcopy(getattr(l1b.classifier, classifier_name))
            self.classifier.add_parameter(classifier, classifier_name)

    def set_unknown_default(self):
        """
        This method can be used to initialize the surface type with unknown values
        :return:
        """
        flag = np.ones(shape=self.classifier.n_records, dtype=np.bool)
        self.surface_type.add_flag(flag, "unknown")

    def set_l1b_land_mask(self, l1b):
        """
        Use this method to transfer the l1b land mask to the l2 surface type classification.

        NOTE: Highly recommended to do this at the end of the surface type classication
              in order to overwrite potential mis-classifications.

        :param l1b:
        :return:
        """
        l1b_land_mask = l1b.surface_type.get_by_name("land")
        self.surface_type.add_flag(l1b_land_mask.flag, "land")

    def has_class(self, name):
        return name in self._classes


# class RickerTC2014(SurfaceTypeClassifier):
#     """
#     Surface Type classification algorithm from
#
#     Ricker, R., Hendricks, S., Helm, V., Skourup, H., and Davidson, M.:
#     Sensitivity of CryoSat-2 Arctic sea-ice freeboard and thickness on
#     radar-waveform interpretation,
#     The Cryosphere, 8, 1607-1622, doi:10.5194/tc-8-1607-2014, 2014.
#     """
#
#     def __init__(self):
#         super(RickerTC2014, self).__init__()
#         self._classes = ["unkown", "ocean", "lead", "sea_ice", "land"]
#
#     def _classify(self, options):
#         self._classify_ocean(options)
#         self._classify_leads(options)
#         self._classify_sea_ice(options)
#
#     def _classify_ocean(self, options):
#         opt = options.ocean
#         parameter = self._classifier
#         ocean = ANDCondition()
#         # Mandatory radar mode flag
#         ocean.add(self._is_radar_mode)
#         # Peakiness Thresholds
#         ocean.add(parameter.peakiness >= opt.peakiness_min)
#         ocean.add(parameter.peakiness <= opt.peakiness_max)
#         # Stack Standard Deviation
#         ssd_threshold = opt.stack_standard_deviation_min
#         ocean.add(parameter.stack_standard_deviation >= ssd_threshold)
#         # Ice Concentration
#         ocean.add(parameter.sic < opt.ice_concentration_min)
#         # OCOG Width
#         ocean.add(parameter.ocog_width >= opt.ocog_width_min)
#         # Done, add flag
#         self._surface_type.add_flag(ocean.flag, "ocean")
#
#     def _classify_leads(self, options):
#         opt = options.lead
#         parameter = self._classifier
#         lead = ANDCondition()
#         # Mandatory radar mode flag
#         lead.add(self._is_radar_mode)
#         # Stack (Beam) parameters
#         lead.add(parameter.peakiness_l >= opt.peakiness_l_min)
#         lead.add(parameter.peakiness_r >= opt.peakiness_r_min)
#         lead.add(parameter.peakiness >= opt.peakiness_min)
#         lead.add(parameter.stack_kurtosis >= opt.stack_kurtosis_min)
#         ssd_threshold = opt.stack_standard_deviation_max
#         lead.add(parameter.stack_standard_deviation < ssd_threshold)
#         # Ice Concentration
#         lead.add(parameter.sic > opt.ice_concentration_min)
#         # Done, add flag
#         self._surface_type.add_flag(lead.flag, "lead")
#
#     def _classify_sea_ice(self, options):
#         opt = options.sea_ice
#         parameter = self._classifier
#         ice = ANDCondition()
#         # Mandatory radar mode flag
#         ice.add(self._is_radar_mode)
#         # Stack (Beam) parameters
#         ice.add(parameter.peakiness_r <= opt.peakiness_r_max)
#         ice.add(parameter.peakiness_l <= opt.peakiness_l_max)
#         ice.add(parameter.peakiness <= opt.peakiness_max)
#         ice.add(parameter.stack_kurtosis < opt.stack_kurtosis_max)
#         # Ice Concentration
#         ice.add(parameter.sic > opt.ice_concentration_min)
#         # Done, add flag
#         self._surface_type.add_flag(ice.flag, "sea_ice")


class SICCI2SurfaceType(Level2ProcessorStep, SurfaceTypeClassifier):
    """
    new and unified surface type classifier for cryosat2 and envisat
    based on similar (pulse peakiness, backscatter, leading edge width)

    #TODO: This could be generalized to automatically compare parameters from config file

    """

    def __init__(self, *args, **kwargs):
        """
        Init the class. All *args and **kwargs will be directed to Level2ProcessorStep as
        SurfaceTypeClassifier does not take any input at __init__
        :param args:
        :param kwargs:
        """

        # Init the parent classes
        Level2ProcessorStep.__init__(self, *args, **kwargs)
        SurfaceTypeClassifier.__init__(self)

        # Set the classes that this classifier will detect
        self._classes = ["unknown", "ocean", "lead", "sea_ice", "land"]

    def execute_procstep(self, l1b, l2):
        """
        The mandatory class for a Level2ProcessorStep.
        :param l1b:
        :param l2:
        :return:
        """

        # Step 1: Transfer classifier / sea ice concentration radar mode
        self.transfer_l1b_classifier(l1b)
        self.classifier.add_parameter(l2.sic, "sic")
        self.classifier.add_parameter(l1b.waveform.radar_mode, "radar_mode")

        # Step 2: Run the classifier for different target surface types and radar modes
        month_num = l2.time[0].month
        for radar_mode in self.cfg.options.keys():

            # Get a boolean array indicating if the observation belongs to a radar mode
            # -> Skip this round if the target radar mode is not found
            is_radar_mode = l1b.waveform.radar_mode == self.radar_modes.get_flag(radar_mode)
            if np.count_nonzero(is_radar_mode) == 0:
                continue

            # Classify ocean
            opt = AttrDict(self.cfg.options[radar_mode]["ocean"])
            self.classify_ocean(opt, is_radar_mode)

            # Classify leads
            opt = AttrDict(self.cfg.options[radar_mode]["lead"])
            self.classify_leads(opt, month_num, is_radar_mode)

            # Classify sea ice
            opt = AttrDict(self.cfg.options[radar_mode]["sea_ice"])
            self.classify_sea_ice(opt, month_num, is_radar_mode)

        # Step 3: Add the l1b land flag
        # This step is done at the end to exclude the land flag being overwritten
        # by mis-classification
        self.set_l1b_land_mask(l1b)

        # Step 4: Update the l2 data object with the classification result
        l2.surface_type = self.surface_type

        # Step 5: Generate an error flag
        # -> all surfaces that are marked as unknown
        error_flag = l2.surface_type.get_by_name("unknown")
        return error_flag.flag

    def classify_ocean(self, opt, is_radar_mode):
        """
        Classify ocean waveforms.
        :param opt:
        :param is_radar_mode:
        :return:
        """

        parameter = self.classifier
        ocean = ANDCondition()

        # Mandatory radar mode flag
        ocean.add(is_radar_mode)

        # Peakiness Thresholds
        ocean.add(parameter.peakiness <= opt.peakiness_max)

        # Ice Concentration
        ocean.add(parameter.sic < opt.ice_concentration_min)

        # Done, add flag
        self.surface_type.add_flag(ocean.flag, "ocean")

    def classify_leads(self, opt, month_num, is_radar_mode):
        """
        Classify leads in sea ice
        :param opt:
        :param month_num:
        :param is_radar_mode:
        :return:
        """

        # Pointer
        parameter = self.classifier

        # All conditions must be fulfilled
        lead = ANDCondition()

        # Mandatory radar mode flag
        lead.add(is_radar_mode)

        # Sigma0
        sea_ice_backscatter_min = self.get_threshold_value(opt, "sea_ice_backscatter_min", month_num)
        lead.add(parameter.sigma0 >= sea_ice_backscatter_min)

        # Leading Edge Width
        leading_edge_width_max = self.get_threshold_value(opt, "leading_edge_width_max", month_num)
        # TODO: Ensure leading_edge_width is always available as a parameter
        if "leading_edge_width_first_half" in self.classifier.parameter_list:
            leading_edge_width = parameter.leading_edge_width_first_half + parameter.leading_edge_width_second_half
        else:
            leading_edge_width = parameter.leading_edge_width
        lead.add(leading_edge_width <= leading_edge_width_max)

        # Pulse Peakiness
        peakiness_min = self.get_threshold_value(opt, "leading_edge_width_max", month_num)
        lead.add(parameter.peakiness >= peakiness_min)

        # Ice Concentration
        ice_concentration_min = self.get_threshold_value(opt, "ice_concentration_min", month_num)
        lead.add(parameter.sic > ice_concentration_min)

        # Done, add flag
        self.surface_type.add_flag(lead.flag, "lead")

    def classify_sea_ice(self, opt, month_num, is_radar_mode):
        """
        Classify waveforms as sea ice
        :param opt: The option attribute dictionary
        :param month_num:
        :param is_radar_mode:
        :return:
        """

        # Pointer
        parameter = self.classifier

        # All conditions must be fulfilled
        ice = ANDCondition()

        # Mandatory radar mode flag
        ice.add(is_radar_mode)

        # Sigma0
        sea_ice_backscatter_min = self.get_threshold_value(opt, "sea_ice_backscatter_min", month_num)
        sea_ice_backscatter_max = self.get_threshold_value(opt, "sea_ice_backscatter_max", month_num)
        ice.add(parameter.sigma0 >= sea_ice_backscatter_min)
        ice.add(parameter.sigma0 <= sea_ice_backscatter_max)

        # Leading Edge Width
        leading_edge_width_min = self.get_threshold_value(opt, "leading_edge_width_min", month_num)
        # TODO: Ensure leading_edge_width is always available as a parameter
        if "leading_edge_width_first_half" in self.classifier.parameter_list:
            leading_edge_width = parameter.leading_edge_width_first_half + parameter.leading_edge_width_second_half
        else:
            leading_edge_width = parameter.leading_edge_width
        ice.add(leading_edge_width >= leading_edge_width_min)

        # Pulse Peakiness
        peakiness_max = self.get_threshold_value(opt, "peakiness_max", month_num)
        ice.add(parameter.peakiness <= peakiness_max)

        # Ice Concentration
        ice_concentration_min = self.get_threshold_value(opt, "ice_concentration_min", month_num)
        ice.add(parameter.sic > ice_concentration_min)

        # Done, add flag
        self.surface_type.add_flag(ice.flag, "sea_ice")

    @staticmethod
    def get_threshold_value(options, name, month_num):
        """
        A unified method to retrieve threshold values from lists (one per month) or a scalar
        :param options: configuration object
        :param name: (str) parameter name (must be attribute of) options
        :param month_num: (int) number of the month (1 - 12)
        :return: threshold value to use for specific month
        """

        # Get the parameters
        option_value = getattr(options, name)

        # Check if list
        if isinstance(option_value, (list, tuple)):
            return option_value[month_num - 1]
        else:
            return option_value

    @property
    def l2_input_vars(self):
        return ["sic"]

    @property
    def l2_output_vars(self):
        return ["surface_type"]

    @property
    def error_bit(self):
        return self.error_flag_bit_dict["surface_type"]

# class SICCI1Envisat(Level2SurfaceTypeClassifier):
#     """
#     Surface Type classification algorithm from
#
#     SICCI code base (Envisat surface type classification)
#     """
#
#     def __init__(self):
#         super(SICCI1Envisat, self).__init__()
#         self._classes = ["unkown", "ocean", "lead", "sea_ice", "land"]
#
#     def _classify(self, options):
#         self._classify_ocean(options)
#         self._classify_leads(options)
#         self._classify_sea_ice(options)
#
#     def _classify_ocean(self, options):
#         opt = options.ocean
#         parameter = self._classifier
#         ocean = ANDCondition()
#         # Mandatory radar mode flag
#         ocean.add(self._is_radar_mode)
#         # Peakiness Thresholds
#         ocean.add(parameter.peakiness_old < opt.pulse_peakiness_max)
#         # Ice Concentration
#         ocean.add(parameter.sic < opt.ice_concentration_min)
#         # Done, add flag
#         self._surface_type.add_flag(ocean.flag, "ocean")
#
#     def _classify_leads(self, options):
#         opt = options.lead
#         parameter = self._classifier
#         lead = ANDCondition()
#         # Mandatory radar mode flag
#         lead.add(self._is_radar_mode)
#         # Stack (Beam) parameters
#         lead.add(parameter.peakiness_old > opt.pulse_peakiness_min)
#         # Ice Concentration
#         lead.add(parameter.sic > opt.ice_concentration_min)
#         # Done, add flag
#         self._surface_type.add_flag(lead.flag, "lead")
#
#     def _classify_sea_ice(self, options):
#         opt = options.sea_ice
#         parameter = self._classifier
#         ice = ANDCondition()
#         # Mandatory radar mode flag
#         ice.add(self._is_radar_mode)
#         # Stack (Beam) parameters
#         ice.add(parameter.peakiness_old < opt.pulse_peakiness_max)
#         # Ice Concentration
#         ice.add(parameter.sic > opt.ice_concentration_min)
#         # Done, add flag
#         self._surface_type.add_flag(ice.flag, "sea_ice")


# class ICESatFarellEtAl2009(Level2SurfaceTypeClassifier):
#     """
#     Surface Type classification algorithm from
#
#     SICCI code base (Envisat surface type classification)
#     """
#
#     def __init__(self):
#         super(ICESatFarellEtAl2009, self).__init__()
#         self._classes = ["unkown", "ocean", "lead", "sea_ice", "land"]
#
#     def _classify(self, options):
#         self._classify_ocean(options)
#         self._classify_leads(options)
#         self._classify_sea_ice(options)
#
#     def _classify_ocean(self, options):
#         opt = options.ocean
#         parameter = self._classifier
#         ocean = ANDCondition()
#         # Ice Concentration
#         ocean.add(parameter.sic < opt.ice_concentration_min)
#         # Done, add flag
#         self._surface_type.add_flag(ocean.flag, "ocean")
#
#     def _classify_leads(self, options):
#         opt = options.lead
#         parameter = self._classifier
#         lead = ANDCondition()
#         # Mandatory radar mode flag
#         lead.add(self._is_radar_mode)
#         # Reflectivity
#         lead.add(parameter.reflectivity <= opt.reflectivity_max)
#         # Echo Gain
#         lead.add(parameter.echo_gain <= opt.echo_gain_max)
#         lead.add(parameter.echo_gain >= opt.echo_gain_min)
# #         lead.add(parameter.echo_gain >= 150.)
#         # Ice Concentration
#         lead.add(parameter.sic > opt.ice_concentration_min)
#         # Done, add flag
#         self._surface_type.add_flag(lead.flag, "lead")
#
#     def _classify_sea_ice(self, options):
#         opt = options.sea_ice
#         parameter = self._classifier
#         ice = ANDCondition()
#         # Mandatory radar mode flag
#         ice.add(self._is_radar_mode)
#         # Reflectivity min
# #        reflectivity = np.array(parameter.reflectivity)
# #        reflectivity[np.isnan(reflectivity)] = 999.
#         ice.add(parameter.reflectivity > opt.reflectivity_min)
#
#         ice.add(parameter.echo_gain <= opt.echo_gain_max)
#
# #        ice.add(parameter.echo_gain < 100.)
#         # Ice Concentration
#         ice.add(parameter.sic > opt.ice_concentration_min)
#         # Done, add flag
#         self._surface_type.add_flag(ice.flag, "sea_ice")


# class ICESatKhvorostovskyTPEnhanced(Level2SurfaceTypeClassifier):
#     """ Classifier based on TC paper from Kirill (lead detection part)
#     which uses coincident local dips in reflectivity and elevation
#     to identify ssh tie points (see Section 3.3.2 An improved algorithm
#     for the TP method). Ocean is identified from sea ice concentration
#     and ice identified as every valid elevation that is neither ocean
#     nor lead.
#
#     Reference:
#         Khvorostovsky, K. and Rampal, P.: On retrieving sea ice freeboard from
#         ICESat laser altimeter, The Cryosphere, 10, 2329-2346,
#         https://doi.org/10.5194/tc-10-2329-2016, 2016.
#     """
#
#     REQUIRED_CLASSIFIERS = [
#             'reflectivity',                        # ICESat refletivity
#             'sea_ice_surface_elevation_corrected'  # ICESat surface elevation
#             'sic',                                 # from l2
#             'mss']                                 # from l2
#
#     CLASSES = ["unkown", "ocean", "lead", "sea_ice"]
#
#     def __init__(self):
#         super(ICESatKhvorostovskyTPEnhanced, self).__init__()
#
#     def _classify(self, options):
#         """ Class API method """
#         self._classify_ocean(options)
#         self._classify_leads(options)
#         self._classify_sea_ice(options)
#
#     def _classify_ocean(self, options):
#         """ Ocean classification based on sea ice concentration only
#         since land will be excluded anyway """
#         opt = options.ocean
#         parameter = self._classifier
#         ocean = ANDCondition()
#         # Ice Concentration
#         ocean.add(parameter.sic < opt.ice_concentration_min)
#         # Done, add flag
#         self._surface_type.add_flag(ocean.flag, "ocean")
#
#     def _classify_leads(self, options):
#         """ Follow the procedure proposed by Kirill: Identification of
#         colocated dips in local elevation & reflectivity """
#
#         # Aliases
#         opt = options.lead
#         parameter = self._classifier
#
#         # Translate window size in km to indices
#         window_size = self.get_filter_width(
#                 opt.filter_width_m,
#                 opt.footprint_spacing_m)
#
#         # get local dips in elevation
#         elevation = parameter.sea_ice_surface_elevation_corrected
#         elevation -= parameter.mss
#         hr, hr_mean, hr_sigma, index_list = self.get_elevation_parameters(
#                 elevation, window_size)
#
#         # get local dips in reflectivity
#         delta_r = self.get_delta_r(
#                 parameter.reflectivity,
#                 window_size,
#                 opt.reflectivity_offset_sdev_factor,
#                 index_list)
#
#         # All criterias needs to be fullfilled
#         lead = ANDCondition()
#
#         # Mandatory radar mode flag
#         lead.add(self._is_radar_mode)
#
#         # Local Reflectivity Minimum
#         lead.add(delta_r >= opt.reflectivity_diff_min)
#
#         # Local Elevation Minimum
#         hr_max = hr_mean - opt.elevation_offset_sdev_factor * hr_sigma
#         lead.add(hr < hr_max)
#
#         # Obligatory Ice Concentration
#         lead.add(parameter.sic > opt.ice_concentration_min)
#
#         # Done, add flag
#         self._surface_type.add_flag(lead.flag, "lead")
# #
# #        import matplotlib.pyplot as plt
# #
# #        x = np.arange(len(elevation))
# #        lead_indices = self._surface_type.lead.indices
# #
# #        plt.figure()
# #        plt.plot(x, hr, label="hr")
# #        plt.plot(x, hr_mean, label="hr_mean")
# #        plt.plot(x, hr_max, label="hr_max")
# #        plt.legend()
#
# #        plt.figure("reflectivity")
# #        plt.plot(x, parameter.reflectivity)
# #        plt.scatter(x[lead_indices], parameter.reflectivity[lead_indices],
# #                    color="red", alpha=0.5)
# #
# #        plt.figure("delta_r")
# #        plt.plot(x, delta_r)
# #        plt.scatter(x[lead_indices], delta_r[lead_indices],
# #                    color="red", alpha=0.5)
# #        plt.hlines(0.3, x[0], x[-1])
# #
# #        plt.figure("elevation")
# #        plt.plot(x, elevation, label="elevation")
# #        plt.plot(x, background_elevation, label="background_elevation")
# #        plt.scatter(x[lead_indices], elevation[lead_indices],
# #                    alpha=0.5, color="red")
# #        plt.legend()
# #
# #        plt.show()
# #        stop
#
#     def _classify_sea_ice(self, options):
#         """ Sea ice is essentially the valid rest (not-ocean & not lead) """
#         opt = options.sea_ice
#         parameter = self._classifier
#         ice = ANDCondition()
#
#         # Mandatory radar mode flag
#         ice.add(self._is_radar_mode)
#
#         # Should not be a lead
#         ice.add(np.logical_not(self._surface_type.lead.flag))
#
#         # High gain value indicates low SNR
#         ice.add(parameter.echo_gain <= opt.echo_gain_max)
#
#         # Ice Concentration
#         ice.add(parameter.sic > opt.ice_concentration_min)
#
#         # Done, add flag
#         self._surface_type.add_flag(ice.flag, "sea_ice")
#
#     def get_filter_width(self, filter_width_m, footprint_spacing_m):
#         filter_width = filter_width_m / footprint_spacing_m
#         # Make sure filter width is odd integer
#         filter_width = np.floor(filter_width) // 2 * 2 + 1
#         filter_width = filter_width.astype(int)
#         return filter_width
#
#     def get_delta_r(self, reflectivity, window_size, sdev_factor, index_list):
#         """ Compute delta_r (\Delta R) as measure for local reflectivity
#         dips """
#
#         # Compute background reflectivity
#         background_reflectivity = np.full(reflectivity.shape, np.nan)
#
#         # Filter reflectivity
#         invalid = np.where(reflectivity > 1)[0]
#         reflectivity[invalid] = np.nan
#
#         # Support Varaibles (Debug purposes only)
#         filter_mean = np.full(reflectivity.shape, np.nan)
#         filter_sdev = np.full(reflectivity.shape, np.nan)
#         filter_threshold = np.full(reflectivity.shape, np.nan)
#
#         # First try: simple loop (is there a better way?)
#         n = len(reflectivity)
#         filter_pad = int((window_size-1)/2)
#         for i in index_list:
#
#             # Get indices
#             i0, i1 = i-filter_pad, i+filter_pad+1
#             i0 = i0 if i0 >= 0 else 0
#             i1 = i1 if i1 <= n-1 else n
#
#             # Get statistics of filter subset
#             reflectivity_subset = reflectivity[i0:i1]
#             filter_mean[i] = np.nanmean(reflectivity_subset)
#             filter_sdev[i] = np.nanstd(reflectivity_subset)
#
#             # Background reflectivity is mean of filter values above
#             # certain threshold to exclude other leads
#             threshold = filter_mean[i] - sdev_factor * filter_sdev[i]
#             filter_threshold[i] = threshold
#             background_values = np.where(reflectivity_subset > threshold)[0]
#             background_reflectivity[i] = np.nanmean(
#                     reflectivity_subset[background_values])
#
#         # Compute local reflectivity offset from background reflectivity
#         delta_r = background_reflectivity - reflectivity
#
#
#
# #        import matplotlib.pyplot as plt
# #        plt.figure("delta_r debug")
# ##        plt.plot(filter_mean, label="filter_mean")
# ##        plt.plot(filter_sdev, label="filter_sdev")
# #        plt.plot(reflectivity, label="reflectivity", alpha=0.5)
# #        plt.plot(filter_threshold, label="filter_threshold")
# #        plt.plot(background_reflectivity, label="background_reflectivity")
# #        plt.plot(background_reflectivity-0.3, linestyle="dashed",
# #                 label="background_reflectivity")
# #        plt.legend()
# #        plt.show()
# #        stop
#
#         return delta_r
#
#     def get_elevation_parameters(self, elevation, window_size):
#         """ The background elevation is defined as the running mean
#         over a filter window (usually 25km) minus a certain fraction
#         of the elevation standard deviation (`sdev_factor`) of this window """
#
#         # Output variables
#         hr = np.full(elevation.shape, np.nan)
#         hr_mean = np.full(elevation.shape, np.nan)
#         hr_sigma = np.full(elevation.shape, np.nan)
#
#         # Only compute for valid elevations
#         index_list = np.where(np.isfinite(elevation))[0]
#
#         n = len(elevation)
#         filter_pad = int((window_size-1)/2)
#
#         # First pass: compute hr
#         for i in index_list:
#             # Get indices
#             i0, i1 = i-filter_pad, i+filter_pad+1
#             i0 = i0 if i0 >= 0 else 0
#             i1 = i1 if i1 <= n-1 else n
#             # Get statistics of filter subset
#             elevation_subset = elevation[i0:i1]
#             hr[i] = elevation[i] - np.nanmean(elevation_subset)
#
#         # second pass: compute hr statistics
#         for i in index_list:
#             # Get indices
#             i0, i1 = i-filter_pad, i+filter_pad+1
#             i0 = i0 if i0 >= 0 else 0
#             i1 = i1 if i1 <= n-1 else n
#             # Get statistics of filter subset
#             hr_subset = hr[i0:i1]
#             hr_mean[i] = np.nanmean(hr_subset)
#             hr_sigma[i] = np.nanstd(hr_subset)
#
#         return hr, hr_mean, hr_sigma, index_list
