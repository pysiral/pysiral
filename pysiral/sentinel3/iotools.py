# -*- coding: utf-8 -*-

import os


def get_sentinel3_l1b_filelist(folder, target="enhanced_measurement.nc"):
    """ Returns a list with measurement.nc files for given month """
    s3_l1b_file_list = []
    for root, dirs, files in os.walk(folder):
        for name in files:
            if name == target:
                s3_l1b_file_list.append(os.path.join(root, name))
    return s3_l1b_file_list
