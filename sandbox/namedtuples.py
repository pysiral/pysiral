# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 20:33:51 2015

@author: Stefan
"""

from collections import namedtuple
import struct


def test_namedtuples():

    Row = namedtuple('Row', ['first', 'second', 'third'])
    A = ['1', '2', '2' '3']


    mds_def = namedtuple('mds_def', ['name', 'fmt', 'size'])
    time_orbit_group = [
        mds_def(name="day", fmt="l", size=4),
        mds_def(name="sec", fmt="l", size=4),
        mds_def(name="msec", fmt="l", size=4),
        mds_def(name="uso_corr", fmt="l", size=4),
        mds_def(name="mode_id", fmt="H", size=2),
        mds_def(name="source_sequence_counter", fmt="H", size=2),
        mds_def(name="instrument_configuration", fmt="L", size=4),
        mds_def(name="burst_counter", fmt="L", size=4),
        mds_def(name="measurement_latitude", fmt="l", size=4),
        mds_def(name="measurement_longtitude", fmt="l", size=4),
        mds_def(name="altitude_cog", fmt="L", size=4),
        mds_def(name="altitude_rate", fmt="L", size=4),
        mds_def(name="satellite_velocity", fmt="3l", size=12),
        mds_def(name="real_beam", fmt="3L", size=12),
        mds_def(name="interferometry_baseline", fmt="3l", size=12),
        mds_def(name="star_tracker_usage", fmt="H", size=2),
        mds_def(name="antenna_roll", fmt="l", size=4),
        mds_def(name="antenna_pitch", fmt="l", size=4),
        mds_def(name="antenna_yaw", fmt="l", size=4),
        mds_def(name="fbr_confidence_flag", fmt="l", size=4),
        mds_def(name="spare", fmt="4c", size=4)]

    struct_parser = ">"
    bytesize = 0
    for field in time_orbit_group:
        struct_parser += field.fmt
        bytesize += field.size
        print field.fmt, field.size, struct_parser, struct.calcsize(struct_parser), bytesize

if __name__ == '__main__':
    test_namedtuples()
