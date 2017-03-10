# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 15:26:26 2016

@author: shendric
"""

import os
import sys
import glob
import shutil
import argparse
import zipfile

from pysiral.path import file_basename, validate_directory


def envisat_sgdr_cycles2month():
    """
    This script is used to reorder *.N1 Envisat RA-2 sgdr files
    from the cycle_xxx based sequence on th ESA ftp to
    the yyyy/mm folder structure required for pysiral

    Has also the ability to unzip all files in the cycle_xxx folders
    as some part of the data was delivered with daily zip files
    """

    # get argument parser
    parser = get_sgdr_cycles2month_argparser()
    args = parser.parse_args()

    # Check if directory
    if not os.path.isdir(args.folder):
        msg = "Invalid folder: %s" % args.folder
        sys.exit(msg)

    # get list of cycle folders
    cycle_folders = glob.glob(os.path.join(args.folder, "cycle_*"))
    print "Found %g cycle folders" % len(cycle_folders)

    for cycle_folder in sorted(cycle_folders):

        print "\nCycle Folder: %s" % cycle_folder.split(os.sep)[-1]

        # Get the number of zip files
        zip_files = glob.glob(os.path.join(cycle_folder, "*.zip"))
        n_zip_files = len(zip_files)
        print " Number of zip files: %g" % n_zip_files

        # Unzip files
        if n_zip_files > 0 and args.unzip:
            print " ...unzip files"
            for zip_file in sorted(zip_files):
                if not args.test_mode:
                    with zipfile.ZipFile(zip_file, "r") as zip_ref:
                        zip_ref.extractall(cycle_folder)
                    if args.delete_zip_files:
                        os.remove(zip_file)
                else:
                    print " %s" % zip_file
                    if args.delete_zip_files:
                        print " delete %s" % zip_file

        # Get the number of N1 files
        n1_files = glob.glob(os.path.join(cycle_folder, "*.N1"))
        n_n1_files = len(n1_files)
        print " Number of sgdr  files: %g" % n_n1_files

        print " ...moving files"

        for n1_file in sorted(n1_files):

            # get the year month for each file
            n1_filename = file_basename(n1_file)
            date_str = n1_filename.split("_")[2][6:]
            yyyy, mm = date_str[0:4], date_str[4:6]

            try:
                year = int(yyyy)
                month = int(mm)
            except:
                print " problem with filenaming for %s" % n1_filename
                continue

            if not year < 2013 and year > 2001:
                print " Year (%g) out of range for envisat" % year
                continue

            if not month < 0 and month > 12:
                print " Month (%g) out of range" % month
                continue

            # Verify output directory (cycles do not map to month)
            output_directory = os.path.join(args.folder, yyyy, mm)
            validate_directory(output_directory)

            target_file = os.path.join(output_directory, n1_filename+".N1")
            if not args.test_mode:
                try:
                    shutil.move(n1_file, target_file)
                except:
                    pass
            else:
                print " Move %s to %s" % (n1_filename, output_directory)


def get_sgdr_cycles2month_argparser():
    """ Return the arg parser for this script """

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--unzip',
        action='store_true',
        dest='unzip',
        default=False)

    parser.add_argument(
        '--delete-zip-files',
        action='store_true',
        dest='delete_zip_files',
        default=False)

    parser.add_argument(
        '--no-copy',
        action='store_true',
        dest='no_copy',
        default=False)

    parser.add_argument(
        '--test-mode',
        action='store_true',
        dest='test_mode',
        default=False)

    # Main product folder
    parser.add_argument(
        action='store',
        dest='folder',
        help='sgdr product folder')

    return parser

if __name__ == "__main__":
    envisat_sgdr_cycles2month()
