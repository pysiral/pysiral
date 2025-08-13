# -*- coding: utf-8 -*-

"""
"""

__author__ = "Stefan Hendricks"

from collections import OrderedDict
from typing import List, Union

import numpy as np

from pysiral.core.legacy_classes import DefaultLoggingClass

from enum import Enum, StrEnum


# The standard ESA surface type flag in L1B data
ESA_SURFACE_TYPE_DICT = {
    "ocean": 0,
    "closed_sea": 1,
    "land_ice": 2,
    "land": 3
}


SURFACE_TYPE_DICT = {
    "unknown": 0,
    "ocean": 1,
    "lead": 2,
    "polynya": 3,
    "sea_ice": 4,
    "closed_sea": 5,
    "land_ice": 6,
    "land": 7,
    "invalid": 8
}

# bit values for a 16-bit integer (multiple choice)
WAVEFORM_CLASSIFICATION_BIT_DICT = {
    "valid": 0,               # Power value usable for retracking ...
    "absolute_maximum": 1,    # range bin with the global maximum
    "first_maximum": 2,       # range bin with first maximum (may also be absolute maximum)
    "leading_edge": 3,        # range bin(s) in the leading edge (excluding maximum)
    "trailing_edge": 4,       # range bin(s) in the trailing edge (excluding maximum)
    "noise_floor": 5,         # valid noise floor in the beginning of the waveform
    "side_lobe_artefact": 6,  # side lobe artefacts (if notable and detectable)
    "fft_artefact": 7,        # fft artefacts in the beginning of the waveform (older altimeters)
    "off_nadir_artefact": 8,  # off-nadir lobe artefacts (if notable and detectable)
    "unclassified": 15        # initial flag, respectively unidentified range bin(s)
}


class Hemispheres(StrEnum):
    """
    Enum for hemispheres.
    """
    NORTH = "nh"
    SOUTH = "sh"
    GLOBAL = "global"

    @classmethod
    def get_choices(cls) -> List[str]:
        return [entry.value for entry in cls]


class BasicProcessingLevels(StrEnum):
    """
    Enum for processing levels with a corresponding pysiral processor
    """
    LEVEL1 = "l1"
    LEVEL2 = "l2"
    LEVEL3 = "l3"


class ProcessingLevels(str, Enum):
    """
    Enum for processing levels.
    """
    LEVEL1 = "l1"
    LEVEL1_PREPROCESSED = "l1p"
    LEVEL2 = "l2"
    LEVEL2_INTERMEDIATE = "l2i"
    LEVEL2_PREPROCESSED = "l2p"
    LEVEL3_COLLATED = "l3c"
    LEVEL3_SUPERCOLLATED = "l3s"


class RadarModes(object):

    flag_dict = {"lrm": 0, "sar": 1, "sin": 2}

    def __init__(self):
        pass

    @classmethod
    def get_flag(cls, mode_name: str) -> Union[int, None]:
        try:
            return cls.flag_dict[mode_name]
        except KeyError:
            return None

    @classmethod
    def get_name(cls, flag: int) -> Union[str, None]:
        return next(
            (
                mode_name
                for mode_name, mode_flag in cls.flag_dict.items()
                if flag == mode_flag
            ),
            None,
        )

    def name(self, index: int) -> str:
        i = list(self.flag_dict.values()).index(index)
        return list(self.flag_dict.keys())[i]

    @property
    def num(self) -> int:
        return len(self.flag_dict.keys())


# TODO: Add marine and land/land ice flag groups
class SurfaceType(DefaultLoggingClass):
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

    def __init__(self):
        """

        """
        super(SurfaceType, self).__init__(self.__class__.__name__)

        self.surface_type_dict = dict(**SURFACE_TYPE_DICT)
        self._surface_type_flags = []
        self._surface_type = None

    def name(self, flag_value: int) -> str:
        """
        Return the flag name for a give flag value

        :param flag_value:

        :return:
        """
        i = list(self.surface_type_dict.values()).index(flag_value)
        return list(self.surface_type_dict.keys())[i]

    def set_flag(self, flag: np.ndarray) -> None:
        self._surface_type = flag

    def add_flag(self, flag: np.ndarray, type_str: str) -> None:
        """
        Add a surface type flag (boolean array)

        :param flag:
        :param type_str:

        :return:
        """

        # Add a surface type flag
        errors = []
        if type_str not in self.surface_type_dict.keys():
            errors.append(ValueError(f"surface type str {type_str} unknown"))

        if self.invalid_n_records(len(flag)):
            errors.append(ValueError(f"invalid number of records: {len(flag)} (must be {self.n_records})"))

        if errors:
            raise Exception(errors)

        # Create Flag keyword if necessary
        if self._surface_type is None:
            self._surface_type = np.zeros(shape=len(flag), dtype=np.int8)

        # Update surface type list
        indices = np.where(flag)[0]
        self._surface_type[indices] = self._get_type_id(type_str)
        self._surface_type_flags.append(type_str)

    def has_flag(self, type_str: str) -> bool:
        return type_str in self._surface_type_flags

    def get_by_name(self, name: str) -> 'FlagContainer':
        if name not in self.surface_type_dict.keys():
            return FlagContainer(np.zeros(shape=self.n_records, dtype=bool))
        type_id = self._get_type_id(name)
        flag = np.array(self._surface_type == type_id)
        return FlagContainer(flag)

    def append(self, annex: 'SurfaceType') -> None:
        self._surface_type = np.append(self._surface_type, annex.flag)

    def set_subset(self, subset_list):
        self._surface_type = self._surface_type[subset_list]

    def fill_gaps(self, corrected_n_records, indices_map):
        """ API gap filler method. Note: Gaps will be filled with
        the nodata=unkown (8) surface_type"""
        dtype = self.flag.dtype
        surface_type_corrected = np.full(corrected_n_records, 8, dtype=dtype)
        surface_type_corrected[indices_map] = self._surface_type
        self._surface_type = surface_type_corrected

    def invalid_n_records(self, n: int) -> bool:
        """ Check if flag array has the correct length """
        return self.n_records is not None and self.n_records != n

    def _get_type_id(self, name: str) -> int:
        return self.surface_type_dict[name]

    @property
    def flag(self) -> np.ndarray:
        return self._surface_type

    @property
    def n_records(self) -> int:
        return None if self._surface_type is None else len(self._surface_type)

    @property
    def dimdict(self) -> 'OrderedDict':
        """ Returns dictionary with dimensions"""
        return OrderedDict([("n_records", self.n_records)])

    @property
    def parameter_list(self) -> List[str]:
        return ["flag"]

    @property
    def lead(self) -> 'FlagContainer':
        return self.get_by_name("lead")

    @property
    def ocean(self) -> 'FlagContainer':
        return self.get_by_name("ocean")

    @property
    def sea_ice(self) -> 'FlagContainer':
        return self.get_by_name("sea_ice")

    @property
    def land(self) -> 'FlagContainer':
        return self.get_by_name("land")


# TODO: Has this been used yet?
class IceType(object):
    """
    Container for ice type information

    Possible classifications
        - young thin ice
        - first year ice
        - multi-year ice
        - wet ice
    """
    _ICE_TYPE_DICT = {
        "thin_ice": 0,
        "first_year_ice": 1,
        "multi_year_ice": 2,
        "wet_ice": 3}

    def __init__(self):
        self._ice_type_flag = None


class FlagContainer(object):

    def __init__(self, flag: np.ndarray = None) -> None:
        self._flag = flag

    def set_flag(self, flag: np.ndarray) -> None:
        self._flag = flag

    def add(self, flag: np.ndarray) -> None:
        raise NotImplementedError("FlagContainer should not be called directly")

    @property
    def indices(self) -> np.ndarray:
        return np.where(self._flag)[0]

    @property
    def flag(self) -> np.ndarray:
        return self._flag

    @property
    def num(self) -> int:
        return len(self.indices)


class ANDCondition(FlagContainer):

    def __init__(self, **kwargs):
        super(ANDCondition, self).__init__(**kwargs)

    def add(self, flag: Union[np.ndarray, bool]) -> None:
        if self._flag is None:
            self.set_flag(flag)
        else:
            self.set_flag(np.logical_and(self.flag, flag))


class ORCondition(FlagContainer):

    def __init__(self, **kwargs):
        super(ORCondition, self).__init__(**kwargs)

    def add(self, flag: Union[np.ndarray, bool]) -> None:
        if self._flag is None:
            self.set_flag(flag)
        else:
            self.set_flag(np.logical_or(self.flag, flag))
