# -*- coding: utf-8 -*-

from pysiral.config import ConfigInfo, get_yaml_config
from pysiral.helper import month_iterator

from pysiral.job import Level3Job
from pysiral.l3proc import Level3Processor

from datetime import datetime
import argparse
import glob
import sys
import os


def pysiral_l3proc():

    """ get the pysiral configuration info """
    config = ConfigInfo()

    """ parse command line arguments """
    parser = get_l3proc_argparser()
    args = parser.parse_args()

    """ Read the settings file """
    setting_file = os.path.join(
        config.pysiral_local_path, "settings", "l3", args.setting_id+".yaml")
    setting = get_yaml_config(setting_file)

    """ validate date values """
    validate_year_month_list(args.start_date, "start date")
    validate_year_month_list(args.stop_date, "stop date")

    """ Start the processing """
    # Assemble the job order
    job = Level3Job()
    job.set_input_directory(args.input_directory)
    job.set_grid_definition(setting.grid_definition)
    job.set_parameter(
        l2=setting.l2_parameter, l3=setting.l3_parameter,
        frb_nanmask=setting.freeboard_nan_mask_targets)
    job.validate()

    # Start the processor
    l3proc = Level3Processor(job)

    iterator = month_iterator(args.start_date[0], args.start_date[1],
                              args.stop_date[0], args.stop_date[1])
    for year, month in iterator:
        l2idata_files = job.get_monthly_l2idata_files(year, month)
        if len(l2idata_files) == 0:
            continue
        l3proc.set_l2_files(l2idata_files)
        l3proc.run()


def validate_year_month_list(year_month_list, label):
    try:
        datetime(year_month_list[0], year_month_list[1],  1)
    except ValueError:
        print "Error: Invalid "+label+" (%04g, %02g)" % (
            year_month_list[0], year_month_list[1])
        sys.exit(1)


def get_l1bdata_files(config, job, hemisphere, year, month):
    l1b_repo = config.local_machine.l1b_repository[job.mission.id].l1bdata
    directory = os.path.join(
        l1b_repo, hemisphere, "%04g" % year, "%02g" % month)
    print directory
    l1bdata_files = sorted(glob.glob(os.path.join(directory, "*.nc")))
    return l1bdata_files


def get_l3proc_argparser():
    """ Handle command line arguments """
    parser = argparse.ArgumentParser()
    # Mission id string: cryosat2, envisat, ...
    parser.add_argument(
        '-s', '-setting', action='store', dest='setting_id',
        help='setting id of yaml file in /settings/l3', required=True)
    # input directory
    parser.add_argument(
        '-i', action='store', dest='input_directory',
        help='link to directory with l2 netcdf files',
        required=True)
    # Start month as list: [yyyy, mm]
    parser.add_argument(
        '-start', action='store', dest='start_date',
        nargs='+', type=int, required=True,
        help='start date as year and month (-start yyyy mm)')
    # Stop month as list: [yyyy, mm]
    parser.add_argument(
        '-stop', action='store', dest='stop_date',
        nargs='+', type=int,  required=True,
        help='start date as year and month (-stop yyyy mm)')
    # Add an Option to skip a number of files (e.g. for a restart)
    parser.add_argument(
        "-skipmonth", action='store', type=int, const=[], nargs='?',
        dest='skip_month', help='month to skip (default: none)')
    # Show debug statements
    parser.add_argument(
        "-v", "--verbose", help="increase output verbosity",
        action="store_true")
    # show preprocessor version
    parser.add_argument(
        '--version', action='version', version='%(prog)s 0.1a')

    return parser


if __name__ == "__main__":
    pysiral_l3proc()
