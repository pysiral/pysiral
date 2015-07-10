# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 20:33:51 2015

@author: Stefan
"""

from collections import namedtuple
import struct


def test_namedtuples():

    mds_def = namedtuple('mds_def', ['name', 'fmt', 'size'])
    time_orbit_group = [
        mds_def(name="day", fmt="l", size=4),
        mds_def(name="sec", fmt="l", size=4),
        mds_def(name="msec", fmt="l", size=4),
        mds_def(name="uso_corr", fmt="l", size=4),
        mds_def(name="mode_id", fmt="H", size=2),
        mds_def(name="source_sequence_counter", fmt="H", size=2),
        mds_def(name="instrument_configuration", fmt="H", size=2),
        mds_def(name="burst_counter", fmt="l", size=4),
        mds_def(name="measurement_latitude", fmt="l", size=4),
        mds_def(name="measurement_longtitude", fmt="l", size=4),
        mds_def(name="altitude_cog", fmt="l", size=4),
        mds_def(name="altitude_rate", fmt="l", size=4),
        mds_def(name="satellite_velocity", fmt="3l", size=12),
        mds_def(name="real_beam", fmt="3l", size=12),
        mds_def(name="interferometry_baseline", fmt="3l", size=12),
        mds_def(name="star_tracker_usage", fmt="H", size=2),
        mds_def(name="antenna_roll", fmt="l", size=4),
        mds_def(name="antenna_pitch", fmt="l", size=4),
        mds_def(name="antenna_yaw", fmt="l", size=4),
        mds_def(name="fbr_confidence_flag", fmt="l", size=4)]

    group_bytesize = sum([field.size for field in time_orbit_group])
    struct_parser = "".join([field.fmt for field in time_orbit_group])
    field_names = ([field.name for field in time_orbit_group])
    print group_bytesize
    print struct_parser
    print field_names

if __name__ == '__main__':
    test_namedtuples()
