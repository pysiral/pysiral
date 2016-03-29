# -*- coding: utf-8 -*-

import numpy as np

NOMINAL_TRACKING_BIN = 45
BIN_WIDTH_METER = 0.4686


def mdsr_timestamp_to_datetime(mdsr_time):
    """ Interpret the Envisat format time into a datetime object """
    from datetime import datetime, timedelta

    epoch = 946684800  # 1:st of January 2000 in POSIX time
    return datetime.fromtimestamp(epoch) + timedelta(
        days=mdsr_time.day,
        seconds=mdsr_time.sec,
        microseconds=mdsr_time.msec)


def get_envisat_wfm_range(window_delay_meter, n_range_bins):
    n_records = len(window_delay_meter)
    wfm_range = np.ndarray(shape=(n_records, n_range_bins), dtype=np.float32)
    for i in range(n_records):
        wfm_range[i, :] = np.arange(n_range_bins)*BIN_WIDTH_METER + \
            window_delay_meter[i]
    return wfm_range


def get_envisat_window_delay(tracker_range, doppler_correction,
                             slope_doppler_correction):
    """
    Calculating the window delay of the first range bin:

    from SICCI-1 processor:

    window_delay = tracker_range + range_delta + doppler_correction +
                   slope_doppler_correction

    with

    tracker_range: field '18hz_tracker_range_no_doppler_ku' from
                   mds_ra2.range_information

    range_delta: - nominal_tracking_bin * bin_width

    doppler_correction: field '18Hz_ku_range_doppler' from
                        mds_ra2.range_correction

    slope_doppler_correction: field '18Hz_ku_range_doppler_slope' from
                              mds_ra2.range_correction

    """

    range_delta = -1.0 * NOMINAL_TRACKING_BIN * BIN_WIDTH_METER
    window_delay_meter = tracker_range + range_delta + doppler_correction + \
        slope_doppler_correction
    return window_delay_meter
