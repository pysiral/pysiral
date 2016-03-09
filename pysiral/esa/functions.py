# -*- coding: utf-8 -*-

import numpy as np


def get_structarr_attr(struct_arr, field):
    """
    Get all attributes from array of objects that support dict notation
    (e.g. struct_array[:].field <- does not work in python)

    """
    dtype = type(struct_arr[0][field])
    data = np.array([record[field] for record in struct_arr], dtype=dtype)
    return data
