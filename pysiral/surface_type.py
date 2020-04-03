# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 11:25:04 2015

@author: Stefan
"""

#TODO: options to __init__ (and not classify())

from pysiral.config import RadarModes
from pysiral.flag import FlagContainer, ANDCondition

import numpy as np
from treedict import TreeDict
from collections import OrderedDict


ESA_SURFACE_TYPE_DICT = {
    "ocean": 0,
    "closed_sea": 1,
    "land_ice": 2,
    "land": 3}


class SurfaceType(object):
    """
    Container for surface type information.

    Possible classifications (Adapted from CryoSat-2 conventions)
        - unknown
        - ocean
        - closed sea/lakes
        - lead
        - large lead/polynya
        - sea ice (general sea ice class, not to be confused with ice type)
        - continental ice
        - land
    """
    SURFACE_TYPE_DICT = {
        "unknown": 0,
        "ocean": 1,
        "lead": 2,
        "polynya": 3,
        "sea_ice": 4,
        "closed_sea": 5,
        "land_ice": 6,
        "land": 7,
        "invalid": 8}

    def __init__(self):
        self._surface_type_flags = []
        self._surface_type = None
        self._n_records = None

    @property
    def flag(self):
        return self._surface_type

    @property
    def n_records(self):
        if self._surface_type is None:
            n_records = 0
        else:
            n_records = len(self._surface_type)
        self._n_records = n_records
        return n_records

    @property
    def dimdict(self):
        """ Returns dictionary with dimensions"""
        dimdict = OrderedDict([("n_records", self.n_records)])
        return dimdict

    @property
    def parameter_list(self):
        return ["flag"]

    @property
    def lead(self):
        return self.get_by_name("lead")

    @property
    def sea_ice(self):
        return self.get_by_name("sea_ice")

    @property
    def land(self):
        return self.get_by_name("land")

    def name(self, index):
        i = self.SURFACE_TYPE_DICT.values().index(index)
        return self.SURFACE_TYPE_DICT.keys()[i]

    def set_flag(self, flag):
        self._surface_type = flag

    def add_flag(self, flag, type_str):
        """ Add a surface type flag """
        if type_str not in self.SURFACE_TYPE_DICT.keys():
            # TODO: Error Handling
            raise("surface type str %s unknown" % type_str)
        if self._invalid_n_records(len(flag)):
            raise("invalid number of records: %g (must be %g)" % (
                len(flag), self._n_records))
        # Create Flag keyword if necessary
        if self._surface_type is None:
            self._n_records = len(flag)
            self._surface_type = np.zeros(
                shape=(self._n_records), dtype=np.int8)
        # Update surface type list
        self._surface_type[np.where(flag)[0]] = self._get_type_id(type_str)
        self._surface_type_flags.append(type_str)

    def has_flag(self, type_str):
        return type_str in self._surface_type_flags

    def get_by_name(self, name):
        if name in self.SURFACE_TYPE_DICT.keys():
            type_id = self._get_type_id(name)
            return FlagContainer(self._surface_type == type_id)
        else:
            return FlagContainer(
                np.zeros(shape=(self._n_records), dtype=np.bool))

    def append(self, annex):
        self._surface_type = np.append(self._surface_type, annex.flag)

    def set_subset(self, subset_list):
        self._surface_type = self._surface_type[subset_list]
        self._n_records = len(subset_list)

    def fill_gaps(self, corrected_n_records, gap_indices, indices_map):
        """ API gap filler method. Note: Gaps will be filled with
        the nodata=unkown (8) surface_type"""
        self._n_records = corrected_n_records
        dtype = self.flag.dtype
        surface_type_corrected = np.full((corrected_n_records), 8, dtype=dtype)
        surface_type_corrected[indices_map] = self._surface_type
        self._surface_type = surface_type_corrected

    def _invalid_n_records(self, n):
        """ Check if flag array has the correct length """
        if self._n_records is None:  # New flag, ok
            return False
        elif self._n_records == n:   # New flag has correct length
            return False
        else:                        # New flag has wrong length
            return True

    def _get_type_id(self, name):
        return self.SURFACE_TYPE_DICT[name]

#    def __getattr__(self, name):
#        """
#        Return empty lists for surface type flags that have not been
#        set yet
#        """
#        print name
#        if name in self.SURFACE_TYPE_DICT.keys():
#            if self.has_flag(name):
#                type_id = self._get_type_id(name)
#                return TypeContainer(self._surface_type == type_id)
#            else:
#                return TypeContainer(
#                    np.zeros(shape=(self._n_records), dtype=np.bool))
#        # Fix of deepcopy bug
#        if name in ["__getnewargs_ex__", "__deepcopy__"]:
#            raise AttributeError("%r has no attribute %r" % (type(self), name))

class IceType(object):
    """
    Container for ice type information

    Possible classifications
        - young thin ice
        - first year ice
        - multi year ice
        - wet ice
    """
    _ICE_TYPE_DICT = {
        "thin_ice": 0,
        "first_year_ice": 1,
        "multi_year_ice": 2,
        "wet_ice": 3}

    def __init__(self):
        self._ice_type_flag = None


class ClassifierContainer(object):

    def __init__(self):
        self.parameter_list = []
        self._n_records = None

    def add_parameter(self, parameter, parameter_name):
        setattr(self, parameter_name, parameter)
        self.parameter_list.append(parameter_name)
        if self._n_records is None:
            self._n_records = len(parameter)

    @property
    def n_records(self):
        return self._n_records


class SurfaceTypeClassifier(object):
    """ Parent Class for surface type classifiers """

    def __init__(self):
        self._surface_type = SurfaceType()
        self._l1b_surface_type = None
        self._classifier = ClassifierContainer()
        self._radar_modes = RadarModes()

    @property
    def result(self):
        return self._surface_type

    def set_options(self, **opt_dict):
        self._options = TreeDict.fromdict(opt_dict, expand_nested=True)

    def add_classifiers(self, classifier, name):
        self._classifier.add_parameter(classifier, name)

    def set_l1b_surface_type(self, l1b_surface_type):
        self._l1b_surface_type = l1b_surface_type

    def set_initial_classification(self, surface_type):
        """ Overwrite classification"""
        self._surface_type = surface_type

    def classify(self, l1b, l2):

        # Add all classifiers from l1bdata
        for classifier_name in l1b.classifier.parameter_list:
            classifier = getattr(l1b.classifier, classifier_name)
            self.add_classifiers(classifier, classifier_name)

        # add year/month from L1b
        self._year = l1b.info.start_time.year
        self._month = l1b.info.start_time.month

        # add sea ice concentration
        self.add_classifiers(l2.sic, "sic")
        self.add_classifiers(l2.sic, "mss")

        # add radar mode
        self.add_classifiers(l1b.waveform.radar_mode, "radar_mode")

        # Initialize with unkown
        self.set_unknown_default()

        # loop over different radar modes
        # Note: This is necessary for CryoSat-2 with mixed SAR/SIN segments
        for radar_mode in l1b.waveform.radar_modes:

            # Obtain radar mode specific options
            # (with failsaife for older settings files)
            if radar_mode in self._options:
                options = self._options[radar_mode]
            else:
                options = self._options

            # get the radar mode flag
            radar_mode_flag = self._radar_modes.get_flag(radar_mode)

            # Create mandatory condition
            self._is_radar_mode = l1b.waveform.radar_mode == radar_mode_flag

            # Classify
            self._classify(options)

        # Keep land information
        # (also overwrite any potential impossible classifications)
        self.set_l1b_land_mask(l1b)

    def has_class(self, name):
        return name in self._classes

    def set_unknown_default(self):
        flag = np.ones(shape=(self._classifier.n_records), dtype=np.bool)
        self._surface_type.add_flag(flag, "unknown")

    def set_l1b_land_mask(self, l1b):
        l1b_land_mask = l1b.surface_type.get_by_name("land")
        self._surface_type.add_flag(l1b_land_mask.flag, "land")


class RickerTC2014(SurfaceTypeClassifier):
    """
    Surface Type classification algorithm from

    Ricker, R., Hendricks, S., Helm, V., Skourup, H., and Davidson, M.:
    Sensitivity of CryoSat-2 Arctic sea-ice freeboard and thickness on
    radar-waveform interpretation,
    The Cryosphere, 8, 1607-1622, doi:10.5194/tc-8-1607-2014, 2014.
    """

    def __init__(self):
        super(RickerTC2014, self).__init__()
        self._classes = ["unkown", "ocean", "lead", "sea_ice", "land"]

    def _classify(self, options):
        self._classify_ocean(options)
        self._classify_leads(options)
        self._classify_sea_ice(options)

    def _classify_ocean(self, options):
        opt = options.ocean
        parameter = self._classifier
        ocean = ANDCondition()
        # Mandatory radar mode flag
        ocean.add(self._is_radar_mode)
        # Peakiness Thresholds
        ocean.add(parameter.peakiness >= opt.peakiness_min)
        ocean.add(parameter.peakiness <= opt.peakiness_max)
        # Stack Standard Deviation
        ssd_threshold = opt.stack_standard_deviation_min
        ocean.add(parameter.stack_standard_deviation >= ssd_threshold)
        # Ice Concentration
        ocean.add(parameter.sic < opt.ice_concentration_min)
        # OCOG Width
        ocean.add(parameter.ocog_width >= opt.ocog_width_min)
        # Done, add flag
        self._surface_type.add_flag(ocean.flag, "ocean")

    def _classify_leads(self, options):
        opt = options.lead
        parameter = self._classifier
        lead = ANDCondition()
        # Mandatory radar mode flag
        lead.add(self._is_radar_mode)
        # Stack (Beam) parameters
        lead.add(parameter.peakiness_l >= opt.peakiness_l_min)
        lead.add(parameter.peakiness_r >= opt.peakiness_r_min)
        lead.add(parameter.peakiness >= opt.peakiness_min)
        lead.add(parameter.stack_kurtosis >= opt.stack_kurtosis_min)
        ssd_threshold = opt.stack_standard_deviation_max
        lead.add(parameter.stack_standard_deviation < ssd_threshold)
        # Ice Concentration
        lead.add(parameter.sic > opt.ice_concentration_min)
        # Done, add flag
        self._surface_type.add_flag(lead.flag, "lead")

    def _classify_sea_ice(self, options):
        opt = options.sea_ice
        parameter = self._classifier
        ice = ANDCondition()
        # Mandatory radar mode flag
        ice.add(self._is_radar_mode)
        # Stack (Beam) parameters
        ice.add(parameter.peakiness_r <= opt.peakiness_r_max)
        ice.add(parameter.peakiness_l <= opt.peakiness_l_max)
        ice.add(parameter.peakiness <= opt.peakiness_max)
        ice.add(parameter.stack_kurtosis < opt.stack_kurtosis_max)
        # Ice Concentration
        ice.add(parameter.sic > opt.ice_concentration_min)
        # Done, add flag
        self._surface_type.add_flag(ice.flag, "sea_ice")


class SICCI2(SurfaceTypeClassifier):
    """
    new and unified surface type classifier for cryosat2 and envisat
    based on similar (pulse peakiness, backscatter, leading edge width)
    """

    def __init__(self):
        super(SICCI2, self).__init__()
        self._classes = ["unkown", "ocean", "lead", "sea_ice", "land"]

    def _classify(self, options):
        self._classify_ocean(options)
        self._classify_leads(options)
        self._classify_sea_ice(options)

    def _classify_ocean(self, options):
        opt = options.ocean
        parameter = self._classifier
        ocean = ANDCondition()
        # Mandatory radar mode flag
        ocean.add(self._is_radar_mode)
        # Peakiness Thresholds
        ocean.add(parameter.peakiness <= opt.peakiness_max)
        # Ice Concentration
        ocean.add(parameter.sic < opt.ice_concentration_min)
        # Done, add flag
        self._surface_type.add_flag(ocean.flag, "ocean")

    def _classify_leads(self, options):
        """
        Classify leads in sea ice
        :param options:
        :return:
        """

        # Pointer
        opt = options.lead
        month_num = self._month
        parameter = self._classifier

        # All conditions must be fulfilled
        lead = ANDCondition()

        # Mandatory radar mode flag
        lead.add(self._is_radar_mode)

        # Sigma0
        sea_ice_backscatter_min = self.get_threshold_value(opt, "sea_ice_backscatter_min", month_num)
        lead.add(parameter.sigma0 >= sea_ice_backscatter_min)

        # Leading Edge Width
        leading_edge_width_max = self.get_threshold_value(opt, "leading_edge_width_max", month_num)
        leading_edge_width = parameter.leading_edge_width_first_half + parameter.leading_edge_width_second_half
        lead.add(leading_edge_width <= leading_edge_width_max)

        # Pulse Peakiness
        peakiness_min = self.get_threshold_value(opt, "leading_edge_width_max", month_num)
        lead.add(parameter.peakiness >= peakiness_min)

        # Ice Concentration
        ice_concentration_min = self.get_threshold_value(opt, "ice_concentration_min", month_num)
        lead.add(parameter.sic > ice_concentration_min)

        # Done, add flag
        self._surface_type.add_flag(lead.flag, "lead")

    def _classify_sea_ice(self, options):
        """
        Classify waveforms as sea ice
        :param options: The option attribute dictionary
        :return:
        """

        # Pointer
        opt = options.sea_ice
        month_num = self._month
        parameter = self._classifier

        # All conditions must be fulfilled
        ice = ANDCondition()

        # Mandatory radar mode flag
        ice.add(self._is_radar_mode)

        # Sigma0
        sea_ice_backscatter_min = self.get_threshold_value(opt, "sea_ice_backscatter_min", month_num)
        sea_ice_backscatter_max = self.get_threshold_value(opt, "sea_ice_backscatter_max", month_num)
        ice.add(parameter.sigma0 >= sea_ice_backscatter_min)
        ice.add(parameter.sigma0 <= sea_ice_backscatter_max)

        # Leading Edge Width
        leading_edge_width_min = self.get_threshold_value(opt, "leading_edge_width_min", month_num)
        leading_edge_width = parameter.leading_edge_width_first_half + parameter.leading_edge_width_second_half
        ice.add(leading_edge_width >= leading_edge_width_min)

        # Pulse Peakiness
        peakiness_max = self.get_threshold_value(opt, "peakiness_max", month_num)
        ice.add(parameter.peakiness <= peakiness_max)

        # Ice Concentration
        ice_concentration_min = self.get_threshold_value(opt, "ice_concentration_min", month_num)
        ice.add(parameter.sic > ice_concentration_min)

        # Done, add flag
        self._surface_type.add_flag(ice.flag, "sea_ice")

    def get_threshold_value(self, options, name, month_num):
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
        if type(option_value) is list:
            return option_value[month_num-1]
        else:
            return option_value


class SICCI1Envisat(SurfaceTypeClassifier):
    """
    Surface Type classification algorithm from

    SICCI code base (Envisat surface type classification)
    """

    def __init__(self):
        super(SICCI1Envisat, self).__init__()
        self._classes = ["unkown", "ocean", "lead", "sea_ice", "land"]

    def _classify(self, options):
        self._classify_ocean(options)
        self._classify_leads(options)
        self._classify_sea_ice(options)

    def _classify_ocean(self, options):
        opt = options.ocean
        parameter = self._classifier
        ocean = ANDCondition()
        # Mandatory radar mode flag
        ocean.add(self._is_radar_mode)
        # Peakiness Thresholds
        ocean.add(parameter.peakiness_old < opt.pulse_peakiness_max)
        # Ice Concentration
        ocean.add(parameter.sic < opt.ice_concentration_min)
        # Done, add flag
        self._surface_type.add_flag(ocean.flag, "ocean")

    def _classify_leads(self, options):
        opt = options.lead
        parameter = self._classifier
        lead = ANDCondition()
        # Mandatory radar mode flag
        lead.add(self._is_radar_mode)
        # Stack (Beam) parameters
        lead.add(parameter.peakiness_old > opt.pulse_peakiness_min)
        # Ice Concentration
        lead.add(parameter.sic > opt.ice_concentration_min)
        # Done, add flag
        self._surface_type.add_flag(lead.flag, "lead")

    def _classify_sea_ice(self, options):
        opt = options.sea_ice
        parameter = self._classifier
        ice = ANDCondition()
        # Mandatory radar mode flag
        ice.add(self._is_radar_mode)
        # Stack (Beam) parameters
        ice.add(parameter.peakiness_old < opt.pulse_peakiness_max)
        # Ice Concentration
        ice.add(parameter.sic > opt.ice_concentration_min)
        # Done, add flag
        self._surface_type.add_flag(ice.flag, "sea_ice")


class ICESatFarellEtAl2009(SurfaceTypeClassifier):
    """
    Surface Type classification algorithm from

    SICCI code base (Envisat surface type classification)
    """

    def __init__(self):
        super(ICESatFarellEtAl2009, self).__init__()
        self._classes = ["unkown", "ocean", "lead", "sea_ice", "land"]

    def _classify(self, options):
        self._classify_ocean(options)
        self._classify_leads(options)
        self._classify_sea_ice(options)

    def _classify_ocean(self, options):
        opt = options.ocean
        parameter = self._classifier
        ocean = ANDCondition()
        # Ice Concentration
        ocean.add(parameter.sic < opt.ice_concentration_min)
        # Done, add flag
        self._surface_type.add_flag(ocean.flag, "ocean")

    def _classify_leads(self, options):
        opt = options.lead
        parameter = self._classifier
        lead = ANDCondition()
        # Mandatory radar mode flag
        lead.add(self._is_radar_mode)
        # Reflectivity
        lead.add(parameter.reflectivity <= opt.reflectivity_max)
        # Echo Gain
        lead.add(parameter.echo_gain <= opt.echo_gain_max)
        lead.add(parameter.echo_gain >= opt.echo_gain_min)
#         lead.add(parameter.echo_gain >= 150.)
        # Ice Concentration
        lead.add(parameter.sic > opt.ice_concentration_min)
        # Done, add flag
        self._surface_type.add_flag(lead.flag, "lead")

    def _classify_sea_ice(self, options):
        opt = options.sea_ice
        parameter = self._classifier
        ice = ANDCondition()
        # Mandatory radar mode flag
        ice.add(self._is_radar_mode)
        # Reflectivity min
#        reflectivity = np.array(parameter.reflectivity)
#        reflectivity[np.isnan(reflectivity)] = 999.
        ice.add(parameter.reflectivity > opt.reflectivity_min)

        ice.add(parameter.echo_gain <= opt.echo_gain_max)

#        ice.add(parameter.echo_gain < 100.)
        # Ice Concentration
        ice.add(parameter.sic > opt.ice_concentration_min)
        # Done, add flag
        self._surface_type.add_flag(ice.flag, "sea_ice")


class ICESatKhvorostovskyTPEnhanced(SurfaceTypeClassifier):
    """ Classifier based on TC paper from Kirill (lead detection part)
    which uses coincident local dips in reflectivity and elevation
    to identify ssh tie points (see Section 3.3.2 An improved algorithm
    for the TP method). Ocean is identified from sea ice concentration
    and ice identified as every valid elevation that is neither ocean
    nor lead.

    Reference:
        Khvorostovsky, K. and Rampal, P.: On retrieving sea ice freeboard from
        ICESat laser altimeter, The Cryosphere, 10, 2329-2346,
        https://doi.org/10.5194/tc-10-2329-2016, 2016.
    """

    REQUIRED_CLASSIFIERS = [
            'reflectivity',                        # ICESat refletivity
            'sea_ice_surface_elevation_corrected'  # ICESat surface elevation
            'sic',                                 # from l2
            'mss']                                 # from l2

    CLASSES = ["unkown", "ocean", "lead", "sea_ice"]

    def __init__(self):
        super(ICESatKhvorostovskyTPEnhanced, self).__init__()

    def _classify(self, options):
        """ Class API method """
        self._classify_ocean(options)
        self._classify_leads(options)
        self._classify_sea_ice(options)

    def _classify_ocean(self, options):
        """ Ocean classification based on sea ice concentration only
        since land will be excluded anyway """
        opt = options.ocean
        parameter = self._classifier
        ocean = ANDCondition()
        # Ice Concentration
        ocean.add(parameter.sic < opt.ice_concentration_min)
        # Done, add flag
        self._surface_type.add_flag(ocean.flag, "ocean")

    def _classify_leads(self, options):
        """ Follow the procedure proposed by Kirill: Identification of
        colocated dips in local elevation & reflectivity """

        # Aliases
        opt = options.lead
        parameter = self._classifier

        # Translate window size in km to indices
        window_size = self.get_filter_width(
                opt.filter_width_m,
                opt.footprint_spacing_m)

        # get local dips in elevation
        elevation = parameter.sea_ice_surface_elevation_corrected
        elevation -= parameter.mss
        hr, hr_mean, hr_sigma, index_list = self.get_elevation_parameters(
                elevation, window_size)

        # get local dips in reflectivity
        delta_r = self.get_delta_r(
                parameter.reflectivity,
                window_size,
                opt.reflectivity_offset_sdev_factor,
                index_list)

        # All criterias needs to be fullfilled
        lead = ANDCondition()

        # Mandatory radar mode flag
        lead.add(self._is_radar_mode)

        # Local Reflectivity Minimum
        lead.add(delta_r >= opt.reflectivity_diff_min)

        # Local Elevation Minimum
        hr_max = hr_mean - opt.elevation_offset_sdev_factor * hr_sigma
        lead.add(hr < hr_max)

        # Obligatory Ice Concentration
        lead.add(parameter.sic > opt.ice_concentration_min)

        # Done, add flag
        self._surface_type.add_flag(lead.flag, "lead")
#
#        import matplotlib.pyplot as plt
#
#        x = np.arange(len(elevation))
#        lead_indices = self._surface_type.lead.indices
#
#        plt.figure()
#        plt.plot(x, hr, label="hr")
#        plt.plot(x, hr_mean, label="hr_mean")
#        plt.plot(x, hr_max, label="hr_max")
#        plt.legend()

#        plt.figure("reflectivity")
#        plt.plot(x, parameter.reflectivity)
#        plt.scatter(x[lead_indices], parameter.reflectivity[lead_indices],
#                    color="red", alpha=0.5)
#
#        plt.figure("delta_r")
#        plt.plot(x, delta_r)
#        plt.scatter(x[lead_indices], delta_r[lead_indices],
#                    color="red", alpha=0.5)
#        plt.hlines(0.3, x[0], x[-1])
#
#        plt.figure("elevation")
#        plt.plot(x, elevation, label="elevation")
#        plt.plot(x, background_elevation, label="background_elevation")
#        plt.scatter(x[lead_indices], elevation[lead_indices],
#                    alpha=0.5, color="red")
#        plt.legend()
#
#        plt.show()
#        stop

    def _classify_sea_ice(self, options):
        """ Sea ice is essentially the valid rest (not-ocean & not lead) """
        opt = options.sea_ice
        parameter = self._classifier
        ice = ANDCondition()

        # Mandatory radar mode flag
        ice.add(self._is_radar_mode)

        # Should not be a lead
        ice.add(np.logical_not(self._surface_type.lead.flag))

        # High gain value indicates low SNR
        ice.add(parameter.echo_gain <= opt.echo_gain_max)

        # Ice Concentration
        ice.add(parameter.sic > opt.ice_concentration_min)

        # Done, add flag
        self._surface_type.add_flag(ice.flag, "sea_ice")

    def get_filter_width(self, filter_width_m, footprint_spacing_m):
        filter_width = filter_width_m / footprint_spacing_m
        # Make sure filter width is odd integer
        filter_width = np.floor(filter_width) // 2 * 2 + 1
        filter_width = filter_width.astype(int)
        return filter_width

    def get_delta_r(self, reflectivity, window_size, sdev_factor, index_list):
        """ Compute delta_r (\Delta R) as measure for local reflectivity
        dips """

        # Compute background reflectivity
        background_reflectivity = np.full(reflectivity.shape, np.nan)

        # Filter reflectivity
        invalid = np.where(reflectivity > 1)[0]
        reflectivity[invalid] = np.nan

        # Support Varaibles (Debug purposes only)
        filter_mean = np.full(reflectivity.shape, np.nan)
        filter_sdev = np.full(reflectivity.shape, np.nan)
        filter_threshold = np.full(reflectivity.shape, np.nan)

        # First try: simple loop (is there a better way?)
        n = len(reflectivity)
        filter_pad = int((window_size-1)/2)
        for i in index_list:

            # Get indices
            i0, i1 = i-filter_pad, i+filter_pad+1
            i0 = i0 if i0 >= 0 else 0
            i1 = i1 if i1 <= n-1 else n

            # Get statistics of filter subset
            reflectivity_subset = reflectivity[i0:i1]
            filter_mean[i] = np.nanmean(reflectivity_subset)
            filter_sdev[i] = np.nanstd(reflectivity_subset)

            # Background reflectivity is mean of filter values above
            # certain threshold to exclude other leads
            threshold = filter_mean[i] - sdev_factor * filter_sdev[i]
            filter_threshold[i] = threshold
            background_values = np.where(reflectivity_subset > threshold)[0]
            background_reflectivity[i] = np.nanmean(
                    reflectivity_subset[background_values])

        # Compute local reflectivity offset from background reflectivity
        delta_r = background_reflectivity - reflectivity



#        import matplotlib.pyplot as plt
#        plt.figure("delta_r debug")
##        plt.plot(filter_mean, label="filter_mean")
##        plt.plot(filter_sdev, label="filter_sdev")
#        plt.plot(reflectivity, label="reflectivity", alpha=0.5)
#        plt.plot(filter_threshold, label="filter_threshold")
#        plt.plot(background_reflectivity, label="background_reflectivity")
#        plt.plot(background_reflectivity-0.3, linestyle="dashed",
#                 label="background_reflectivity")
#        plt.legend()
#        plt.show()
#        stop

        return delta_r

    def get_elevation_parameters(self, elevation, window_size):
        """ The background elevation is defined as the running mean
        over a filter window (usually 25km) minus a certain fraction
        of the elevation standard deviation (`sdev_factor`) of this window """

        # Output variables
        hr = np.full(elevation.shape, np.nan)
        hr_mean = np.full(elevation.shape, np.nan)
        hr_sigma = np.full(elevation.shape, np.nan)

        # Only compute for valid elevations
        index_list = np.where(np.isfinite(elevation))[0]

        n = len(elevation)
        filter_pad = int((window_size-1)/2)

        # First pass: compute hr
        for i in index_list:
            # Get indices
            i0, i1 = i-filter_pad, i+filter_pad+1
            i0 = i0 if i0 >= 0 else 0
            i1 = i1 if i1 <= n-1 else n
            # Get statistics of filter subset
            elevation_subset = elevation[i0:i1]
            hr[i] = elevation[i] - np.nanmean(elevation_subset)

        # second pass: compute hr statistics
        for i in index_list:
            # Get indices
            i0, i1 = i-filter_pad, i+filter_pad+1
            i0 = i0 if i0 >= 0 else 0
            i1 = i1 if i1 <= n-1 else n
            # Get statistics of filter subset
            hr_subset = hr[i0:i1]
            hr_mean[i] = np.nanmean(hr_subset)
            hr_sigma[i] = np.nanstd(hr_subset)

        return hr, hr_mean, hr_sigma, index_list


def get_surface_type_class(name):
    return globals()[name]()
