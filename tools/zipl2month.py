# -*- coding: utf-8 -*-
"""
Created on Wed May 18 12:22:13 2016

@author: shendric
"""

import re
import os
import sys
import time
import zipfile
import argparse
from logbook import Logger, StreamHandler

# set up the logging instance
StreamHandler(sys.stdout).push_application()
log = Logger("zipl2month")


def zipl2month():

    # Get command line arguments
    argparser = get_zipl2month_argparser()
    args = argparser.parse_args()

    # Validate product folder
    error_flag = validate_product_folder(args.folder)
    if error_flag:
        sys.exit("Invalid folder: %s" % args.folder)

    # Get a list of years in the product folder
    top_level_folders = os.walk(args.folder).next()[1]
    regex = re.compile(r"^[1-2]\d{3}$")
    year_folders = [d for d in top_level_folders if regex.match(d)]

    # loop over years
    for year in year_folders:

        # Extract month folders
        year_sub_folder = os.path.join(args.folder, year)
        second_level_folders = os.walk(year_sub_folder).next()[1]
        regex = re.compile(r"^[0-1][0-9]$")
        months = [d for d in second_level_folders if regex.match(d)]
        months = [d for d in months if int(d) < 13 and int(d) > 0]

        # Loop over month folders
        for month in months:
            t0 = time.time()
            log.info("+ Create archive for: %s-%s" % (year, month))
            folder = os.path.join(args.folder, year, month)
            zip_filename = os.path.join(args.folder, year, month+".zip")
            archive = zipfile.ZipFile(zip_filename, 'w', zipfile.ZIP_DEFLATED)
            zipdir(folder, archive)
            archive.close()
            log.info("- done in %.1f seconds" % (time.time() - t0))


def get_zipl2month_argparser():
    parser = argparse.ArgumentParser()
    # Mission id string: cryosat2, envisat, ...
    parser.add_argument(
        action='store',
        dest='folder',
        help='main product folder')
    return parser


def validate_product_folder(folder):
    is_folder = os.path.isdir(folder)
    return not is_folder


def zipdir(path, ziph):
    # ziph is zipfile handle
    for root, dirs, files in os.walk(path):
        log.info("- found %g files" % (len(files)))
        for file in files:
            filename = os.path.join(root, file)
            arcname = os.sep.join(filename.split(os.sep)[-2:])
            ziph.write(filename, arcname)

if __name__ == "__main__":
    zipl2month()
