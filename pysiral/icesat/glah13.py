# -*- coding: utf-8 -*-
"""
Created on Sun Aug 20 14:33:31 2017

@author: Stefan
"""

from pysiral.errorhandler import ErrorStatus
from pysiral.logging import DefaultLoggingClass

import numpy as np
import h5py


class GLAH13HDF(DefaultLoggingClass):

    def __init__(self, filename):

        # Init class
        class_name = self.__class__.__name__
        super(GLAH13HDF, self).__init__(class_name)
        self.error = ErrorStatus(caller_id=class_name)

        # Init class properties
        self._global_attributes = {}
        self._data_groups = {}

        # Parse filename
        self._filename = filename
        self._parse_glah13()

    def _parse_glah13(self):
        """ Simple parser: Maps alls attributes and the Data_40Hz group
        into a dictionary """

        # Open the hdf5 files
        dataset = h5py.File(self.filename, "r", libver='latest')

        # Store all global attributes in dictionary
        for name, value in dataset.attrs.items():
            self._global_attributes[str(name)] = value

        # NOTE: We are only interested in the 40Hz data
        data_40hz = dataset[u'Data_40HZ']

        for group_name in data_40hz.keys():

            parameter_dict = dict()
            group = data_40hz[group_name]

            # some children of Data_40Hz might be datasets and not
            # a group
            if isinstance(group, h5py.Group):
                for name, value in group.items():
                    parameter_dict[name] = self.get_dset_data(group, name)
            self._data_groups[group_name] = parameter_dict

        # Store additional variables from the 1Hz data
        # (needed for track allocation)
        target_dict_1Hz = {
              "timestamp_1Hz": u"Data_1HZ/Time/d_UTCTime_1",
              "longitude_1Hz": u"Data_1HZ/Geolocation/d_lon",
              "latitude_1Hz": u"Data_1HZ/Geolocation/d_lat",
              "track_id_1Hz": u"Data_1HZ/Geolocation/i_track",
              "is_icesheet_1Hz": u"Data_1HZ/Elevation_Flags/surf_is_flg",
              "is_land_1Hz": u"Data_1HZ/Elevation_Flags/surf_ld_flg",
              "is_ocean_1Hz": u"Data_1HZ/Elevation_Flags/surf_oc_flg",
              "is_seaice_1Hz": u"Data_1HZ/Elevation_Flags/surf_si_flg",
              "reflect_corr_1Hz": u"Data_1HZ/Reflectivity/d_reflCor_atm"}

        for target in target_dict_1Hz.keys():
            dataset_id = target_dict_1Hz[target]
            value = self.get_dset_data(dataset, dataset_id)
            setattr(self, target, value)

        # close dataset
        # Note: we cannot use the with h5py.File(..) as dataset, due to the
        # use of the get_dset_data method (at least I think)
        dataset.close()

    def get_attr(self, name):
        return self.global_attrs_dict.get(name, None)

    def get_parameter(self, group, name):
        try:
            return self._data_groups[group][name]
        except KeyError:
            return None

    def get_dset_data(self, dataset, data_id, fillvalue_to_nan=True):
        """ Returns the parameter value from a give dataset and
        h5 style data id. (optional: fillvalue replaced with nan) """
        dset = dataset[data_id]
        data = dset.value
        # nan's can only be set for floats
        if fillvalue_to_nan and isinstance(data[0], np.floating):
            fillvalue = dset.attrs.get("_FillValue")
            if fillvalue is not None:
                data[data == fillvalue[0]] = np.nan

        return data

    @property
    def data_groups_list(self):
        return sorted(self._data_groups.keys())

    @property
    def filename(self):
        return self._filename

    @property
    def global_attrs_dict(self):
        return self._global_attributes

    @property
    def product_version(self):
        version = self.get_attr("identifier_product_type")
        version += " v"+self.get_attr("identifier_product_format_version")
        return version
