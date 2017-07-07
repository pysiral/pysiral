# -*- coding: utf-8 -*-
"""
Created on Fri Sep 09 17:33:45 2016

@author: Stefan Hendricks

This module is dedicatet to convert between different time standards
"""

from pysiral.config import get_pysiral_local_path
from datetime import datetime, timedelta

import numpy as np

import os
import re


class UTCTAIConverter(object):

    def __init__(self):
        self._timestamp_units = "seconds since 1 January 1900, 00:00:00"
        self.epoch = datetime(1900, 1, 1)
        self.url = r"http://www.ietf.org/timezones/data/leap-seconds.list"
        self._get_leap_seconds_from_config_file()

    def tai2utc(self, tai_datetimes, monotonically=True, check_all=False):
        """ Converts TAI datetime into UTC datetimes """

        # prepare inputs
        tai_datetimes = np.asarray(tai_datetimes)
        utc_datetimes = np.ndarray(shape=tai_datetimes.shape, dtype=object)

        # Use leap seconds for first array entry
        if monotonically and not check_all:
            start_time = tai_datetimes[0]
            single_ls = self._get_leap_seconds_for_utc_time(start_time)
            leap_seconds = np.full(tai_datetimes.shape, single_ls)

        # use leap second from earliest timestamp
        elif not monotonically and not check_all:
            start_time = np.amin(tai_datetimes)
            single_ls = self._get_leap_seconds_for_utc_time(start_time)
            leap_seconds = np.full(tai_datetimes.shape, single_ls)

        # Compute leap second for each entry
        else:
            leap_seconds = np.ndarray(shape=tai_datetimes.shape, dtype=int)
            for i in np.arange(len(tai_datetimes)):
                leap_seconds[i] = self._get_leap_seconds_for_utc_time(
                    tai_datetimes[i])

        # Apply leap seconds
        for i, tai_time in enumerate(tai_datetimes):
            utc_datetimes[i] = tai_time - timedelta(seconds=leap_seconds[i])

        return utc_datetimes

    def update_definition(self):
        """ Get definition file from web """
        # XXX: Requires error catching
        import urllib2
        req = urllib2.Request(self.url)
        response = urllib2.urlopen(req, timeout=60)
        content = response.readlines()
        with open(self.local_ls_ietf_definition, "w") as fhandle:
            for line in content:
                fhandle.write(line)

    def _get_leap_seconds_from_config_file(self):
        """
        Read in leap seconds info from config/leap_seconds.list
        source: (IETF) http://www.ietf.org/timezones/data/leap-seconds.list
        """

        # Pull definition file if not available
        if not os.path.isfile(self.local_ls_ietf_definition):
            self.update_definition()

        # Read from local file
        with open(self.local_ls_ietf_definition, "r") as fhandle:
            content = fhandle.readlines()

        # Parse content
        self.epoch_seconds = np.array([], dtype=np.int64)
        self.leap_seconds = np.array([], dtype=np.int)
        self.leap_seconds_timestamp = np.array([], dtype=object)
        for line in content:
            if re.match("^#", line):
                continue
            arr = line.strip().split()
            self.epoch_seconds = np.append(self.epoch_seconds, int(arr[0]))
            self.leap_seconds = np.append(self.leap_seconds, int(arr[1]))
            # Compute datetime timestamp of leap seconds occurence
            timestamp = self.epoch + timedelta(seconds=int(arr[0]))
            self.leap_seconds_timestamp = np.append(
                self.leap_seconds_timestamp, timestamp)

    def _get_leap_seconds_for_utc_time(self, datetime):
        """ Returns applicable leap seconds for given datetime """
        # find the closet leap seconds change
        timedeltas = datetime - self.leap_seconds_timestamp
        positive_day_offsets = np.array([td.days > 0 for td in timedeltas])
        applicable_ls_indices = np.where(positive_day_offsets)[0]
        try:
            return self.leap_seconds[applicable_ls_indices[-1]]
        except:
            return 0

    @property
    def local_ls_ietf_definition(self):
        """ Return the local filename for the IETF leap seconds definition """
        return os.path.join(get_pysiral_local_path(), "config",
                            "leap-seconds.list")


if __name__ == "__main__":
    # XXX: test code
    converter = UTCTAIConverter()
    tai = np.array([datetime(1981, 1, 2)], dtype=object)
    utc = converter.tai2utc(tai)
    print tai
    print utc
