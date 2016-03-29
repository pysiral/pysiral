# -*- coding: utf-8 -*-

import numpy as np
from netCDF4 import Dataset


def test_netcdf_write_in_loop():
    filename = r"D:\awi\product\altimetry\pysiral\tests\l1b_ncfiles\test.nc"
    nc = Dataset(filename, "w")
    group = nc.createGroup("testgroup")
    n_records = group.createDimension("n_records", 1000)
    for parameter in ["a", "b"]:
        data = np.arange(1000)
        dim = ("n_records",)
        var = group.createVariable(parameter, data.dtype.str, dim)
        var[:] = data
        # setattr(group, parameter, var)
    nc.close()

    nc = Dataset(filename, "r")
    a = nc.groups["testgroup"].variables["a"][:]
    print a

if __name__ == "__main__":
    test_netcdf_write_in_loop()