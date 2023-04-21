# -*- coding: utf-8 -*-
"""
Created on Fri Sep 09 17:33:45 2016

@author: Stefan Hendricks

This module is dedicatet to convert between different time standards
"""

import time
from datetime import datetime

from dateutil.relativedelta import relativedelta


class StopWatch(object):

    def __init__(self):
        self.t0 = None
        self.t1 = None
        self.reset()

    def reset(self):
        self.t0 = None
        self.t1 = None

    def start(self):
        self.t0 = time.time()
        return self

    def stop(self):
        self.t1 = time.time()

    def get_seconds(self):
        return self.t1 - self.t0

    def get_duration(self, fmt="%H:%M:%S"):
        # Example time
        datum = datetime(1900, 1, 1)
        elapsed_seconds = self.t1 - self.t0
        duration = datum + relativedelta(seconds=elapsed_seconds)
        return duration.strftime(fmt)
