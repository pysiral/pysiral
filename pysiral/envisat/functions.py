# -*- coding: utf-8 -*-


def mdsr_timestamp_to_datetime(mdsr_time):
    """ Interpret the Envisat format time into a datetime object """
    from datetime import datetime, timedelta

    epoch = 946684800  # 1:st of January 2000 in POSIX time
    return datetime.fromtimestamp(epoch) + timedelta(
        days=mdsr_time.day,
        seconds=mdsr_time.sec,
        microseconds=mdsr_time.msec)
