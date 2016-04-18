# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 11:25:04 2015

@author: Stefan
"""
import numpy as np
from treedict import TreeDict


class ANDCondition(object):

    def __init__(self):
        self.flag = None

    def add(self, flag):
        if self.flag is None:
            self.flag = flag
        else:
            self.flag = np.logical_and(self.flag, flag)


class TypeContainer(object):

    def __init__(self, flag):
        self._flag = flag

    @property
    def indices(self):
        return np.where(self._flag)[0]

    @property
    def flag(self):
        return self._flag

    @property
    def num(self):
        return len(self.indices)


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
    _SURFACE_TYPE_DICT = {
        "unknown": 0,
        "ocean": 1,
        "lead": 2,
        "polynya": 3,
        "sea_ice": 4,
        "closed_sea": 5,
        "land_ice": 6,
        "land": 7}

    def __init__(self):
        self._surface_type_flags = []
        self._surface_type = None
        self._n_records = None

    @property
    def flag(self):
        return self._surface_type

    @property
    def n_records(self):
        return self._n_records

    def name(self, index):
        i = self._SURFACE_TYPE_DICT.values().index(index)
        return self._SURFACE_TYPE_DICT.keys()[i]

    def add_flag(self, flag, type_str):
        """ Add a surface type flag """
        if type_str not in self._SURFACE_TYPE_DICT.keys():
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
        if name in self._SURFACE_TYPE_DICT.keys():
            type_id = self._get_type_id(name)
            return TypeContainer(self._surface_type == type_id)
        else:
            return TypeContainer(
                np.zeros(shape=(self._n_records), dtype=np.bool))

    def append(self, annex):
        self._surface_type = np.append(self._surface_type, annex.flag)

    def set_subset(self, subset_list):
        self._surface_type = self._surface_type[subset_list]
        self._n_records = len(subset_list)

    def _invalid_n_records(self, n):
        """ Check if flag array has the correct length """
        if self._n_records is None:  # New flag, ok
            return False
        elif self._n_records == n:   # New flag has correct length
            return False
        else:                        # New flag has wrong length
            return True

    def _get_type_id(self, name):
        return self._SURFACE_TYPE_DICT[name]

#    def __getattr__(self, name):
#        """
#        Return empty lists for surface type flags that have not been
#        set yet
#        """
#        print name
#        if name in self._SURFACE_TYPE_DICT.keys():
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


class SurfaceTypeClassifier(object):
    """ Parent Class for surface type classifiers """
    def __init__(self):
        self._surface_type = SurfaceType()

    @property
    def result(self):
        return self._surface_type

    def set_options(self, **opt_dict):
        self._options = TreeDict.fromdict(opt_dict, expand_nested=True)

    def set_classifiers(self, classifiers):
        self._classifier = classifiers

    def set_initial_classification(self, surface_type):
        """ Overwrite classification"""
        self._surface_type = surface_type

    def classify(self):
        self._classify()

    def has_class(self, name):
        return name in self._classes


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
        self._classes = ["unkown", "ocean", "lead", "sea_ice"]

    def _classify(self):
        self._set_unknown_default()
        self._classify_ocean()
        self._classify_leads()
        self._classify_sea_ice()

    def _set_unknown_default(self):
        flag = np.ones(shape=(self._classifier.n_records), dtype=np.bool)
        self._surface_type.add_flag(flag, "unknown")

    def _classify_ocean(self):
        opt = self._options.ocean
        l1b = self._classifier
        ocean = ANDCondition()
        # Peakiness Thresholds
        ocean.add(l1b.peakiness >= opt.peakiness_min)
        ocean.add(l1b.peakiness <= opt.peakiness_max)
        # Stack Standard Deviation
        ssd_threshold = opt.stack_standard_deviation_min
        ocean.add(l1b.stack_standard_deviation >= ssd_threshold)
        # Ice Concentration
        # ocean.add(l1b.ice_concentration < opt.ice_concentration_min)
        # OCOG Width
        ocean.add(l1b.ocog_width >= opt.ocog_width_min)
        # Done, add flag
        self._surface_type.add_flag(ocean.flag, "ocean")

    def _classify_leads(self):
        opt = self._options.lead
        l1b = self._classifier
        lead = ANDCondition()
        # Stack (Beam) parameters
        lead.add(l1b.peakiness_l >= opt.peakiness_l_min)
        lead.add(l1b.peakiness_r >= opt.peakiness_r_min)
        lead.add(l1b.peakiness >= opt.peakiness_min)
        lead.add(l1b.stack_kurtosis >= opt.stack_kurtosis_min)
        ssd_threshold = opt.stack_standard_deviation_max
        lead.add(l1b.stack_standard_deviation < ssd_threshold)
        # Ice Concentration
        # lead.add(l1b.ice_concentration > opt.ice_concentration_max)
        # Done, add flag
        self._surface_type.add_flag(lead.flag, "lead")

    def _classify_sea_ice(self):
        opt = self._options.sea_ice
        l1b = self._classifier
        ice = ANDCondition()
        # Stack (Beam) parameters
        ice.add(l1b.peakiness_r <= opt.peakiness_r_max)
        ice.add(l1b.peakiness_l <= opt.peakiness_l_max)
        ice.add(l1b.peakiness <= opt.peakiness_max)
        ice.add(l1b.stack_kurtosis < opt.stack_kurtosis_max)
        # Ice Concentration
        # ice.add(l1b.ice_concentration > opt.ice_concentration_min)
        # Done, add flag
        self._surface_type.add_flag(ice.flag, "sea_ice")
