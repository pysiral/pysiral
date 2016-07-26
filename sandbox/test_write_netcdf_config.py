# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 17:10:04 2016

@author: shendric
"""

from pysiral.config import ConfigInfo, get_yaml_config

from netCDF4 import Dataset
import os

config = ConfigInfo()

# Read example configuration data
setting_file = os.path.join(
    config.pysiral_local_path, "settings", "l2", "cryosat2_north_awi.yaml")
setting = get_yaml_config(setting_file)

rootgrp = Dataset(r"test_write_netcdf_config.nc", "w")
for item in setting.iterkeys():
    rootgrp.setncattr(item, str(setting[item]))
rootgrp.close()
