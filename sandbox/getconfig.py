# -*- coding: utf-8 -*-
"""
Created on Mon Jul 06 14:41:10 2015

@author: Stefan
"""

from pysiral.config import ConfigInfo


def get_config():

    info = ConfigInfo()
    print info.mission.makeReport()


if __name__ == "__main__":
    get_config()
