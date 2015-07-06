# -*- coding: utf-8 -*-
"""
Created on Mon Jul 06 14:41:10 2015

@author: Stefan
"""

from pysiral.config import ConfigInfo


def get_config():

    info = ConfigInfo()
    print info.mission.makeReport(),"\n"
    print info.local_machine.makeReport(),"\n"
    print info.area.makeReport(),"\n"
    print info.product.makeReport(),"\n"
    print info.auxdata.makeReport(),"\n"
    print info.parameter.makeReport(),"\n"



if __name__ == "__main__":
    get_config()
