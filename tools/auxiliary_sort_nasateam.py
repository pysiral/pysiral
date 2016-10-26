# -*- coding: utf-8 -*-

from pysiral.path import validate_directory, file_basename

import argparse
import glob
import sys
import os


def auxiliary_sort_nasateam():
    """
    Moves UHH nasateam myi fraction netCDF files into yyyy/mm directory
    structure
    """

    # Save space by excluding arctic summer month
    exclude_month = ["05", "06", "07", "08", "09"]

    # Get command line options (path of all myi fraction ice type files)
    path = get_args()

    # Make a file listing
    ncfiles = glob.glob(os.path.join(path, "*.nc"))

    # loop over files, detect year/month from filename and move to
    # corresponding target directory
    for ncfile in ncfiles:
        basename = file_basename(ncfile)
        yyyymmdd = basename.split("_")[8]
        yyyy, mm = yyyymmdd[0:4], yyyymmdd[4:6]
        destination = os.path.join(path, yyyy, mm)
        if mm in exclude_month:
            os.remove(ncfile)
            continue
        validate_directory(destination)
        os.rename(ncfile, os.path.join(destination, basename+".nc"))


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(action='store', dest='path')
    args = parser.parse_args()
    if os.path.isdir(args.path):
        return args.path
    else:
        msg = "Invalid directory: %s" % str(args.path)
        sys.exit(msg)


if __name__ == "__main__":
    auxiliary_sort_nasateam()
