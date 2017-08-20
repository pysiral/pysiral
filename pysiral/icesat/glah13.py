# -*- coding: utf-8 -*-
"""
Created on Sun Aug 20 14:33:31 2017

@author: Stefan
"""

from pysiral.errorhandler import ErrorStatus
from pysiral.logging import DefaultLoggingClass

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
        with h5py.File(self.filename, "r", libver='latest') as dataset:

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
                        parameter_dict[name] = value[:]
                self._data_groups[group_name] = parameter_dict

            # Store additional variables from the 1Hz data
            # (needed for track allocation)
            self.timestamp_1Hz = dataset[u"Data_1HZ/Time/d_UTCTime_1"][:]
            self.longitude_1Hz = dataset[u"Data_1HZ/Geolocation/d_lon"][:]
            self.latitude_1Hz = dataset[u"Data_1HZ/Geolocation/d_lat"][:]
            self.track_id_1Hz = dataset[u"Data_1HZ/Geolocation/i_track"][:]

    def get_attr(self, name):
        return self.global_attrs_dict.get(name, None)

    def get_parameter(self, group, name):
        try:
            return self._data_groups[group][name]
        except KeyError:
            return None

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
