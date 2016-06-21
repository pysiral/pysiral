# -*- coding: utf-8 -*-

import os
import glob


def get_sentinel3_l1b_filelist(folder, target="enhanced_measurement.nc"):
    """ Returns a list with measurement.nc files for given month """
    s3_l1b_file_list = []
    for root, dirs, files in os.walk(folder):
        for name in files:
            if name == target:
                s3_l1b_file_list.append(os.path.join(root, name))
    return s3_l1b_file_list


def get_sentinel3_sral_l1_from_l2(l2_filename, target="measurement.nc"):
    """ Returns the corresponding sral l1 file to a given l2 filename """

    # XXX: This is based on the data structure of the early access
    #      for expert users

    # Step 1: Replace product tag in folder structure
    l1nc_filename = l2_filename.replace("SR_2_LAN", "SR_1_SRA")

    # folder name for L1 and L2 data files are different, need to replace
    # one date tag with asterisk and search for match

    # split the directories
    folder, filename = os.path.split(l1nc_filename)
    directories = folder.split(os.sep)

    # split the last directory name and replace changing date with asterisk
    s3_orbit_dir = directories[-1]
    s3_orbit_dir_components = s3_orbit_dir.split("_")
    s3_orbit_dir_components[9] = "*"

    # Compile the folder again with search pattern
    search_s3_l1b_folder = "_".join(s3_orbit_dir_components)
    main_dir = os.sep.join(directories[:-1])
    sral_l1_search = os.path.join(main_dir, search_s3_l1b_folder)

    # Search and return first match. If no match is found return none
    sral_l1_folder = glob.glob(sral_l1_search)
    if len(sral_l1_folder) > 0:
        l1nc_filename = os.path.join(sral_l1_folder[0], target)
    else:
        return None
    return l1nc_filename
