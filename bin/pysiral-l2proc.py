# -*- coding: utf-8 -*-

from pysiral.config import ConfigInfo, get_yaml_config
from pysiral.helper import month_iterator

from pysiral.job import Level2Job
from pysiral.l2proc import Level2Processor
from pysiral.iotools import get_l1bdata_files

from datetime import datetime
import argparse
import sys
import os


def pysiral_l2proc():

    """ get the pysiral configuration info """
    config = ConfigInfo()

    """ parse command line arguments """
    parser = get_l2proc_argparser()
    args = parser.parse_args()

    """ Read the settings file """
    # is filename
    if os.path.isfile(args.setting_id):
        setting_file = args.setting_id
    # if not filename, than it need to be id of settings file in
    # pysiral\config\l2
    else:
        setting_file = os.path.join(
            config.pysiral_local_path, "settings", "l2",
            args.setting_id+".yaml")
        if not os.path.isfile(setting_file):
            error_message = "Unknown l2 settings file: %s" % setting_file
            sys.exit(error_message)
    setting = get_yaml_config(setting_file)

    """ Add run tag to settings """
    setting.level2.run_tag = args.run_tag

#    """ validate mission id """
#    if args.mission_id not in config.mission_ids:
#        print "Error: Invalid mission id (%s)" % args.mission_id
#        print "use: " + " | ".join(config.mission_ids)
#        sys.exit(1)

    """ validate date values """
    validate_year_month_list(args.start_date, "start date")
    validate_year_month_list(args.stop_date, "stop date")

    """ Check if data are out of bounds of mission lifetime """
    mission_info = config.get_mission_info(setting.mission.id)
    start_date = datetime(args.start_date[0], args.start_date[1], 1)
    stop_date = datetime(args.stop_date[0], args.stop_date[1], 1)
    # Set start time to mission start time (if necessary)
    if start_date < mission_info.data_period.start:
        print "Warning: start date before %s data period %s" % (
            args.mission_id, str(mission_info.data_period.start))
        args.start_date[0] = mission_info.data_period.start.year
        args.start_date[1] = mission_info.data_period.start.month
    # Stop time can be None if mission is ongoing
    if mission_info.data_period.stop is None:
        mission_stop_time = datetime.utcnow()
    else:
        mission_stop_time = mission_info.data_period.stop
    if stop_date > mission_stop_time:
        print "Warning: stopt date late then %s data period %s" % (
            args.mission_id, str(mission_stop_time))
        args.stop_date[0] = mission_stop_time.data_period.start.year
        args.stop_date[1] = mission_stop_time.data_period.start.month

    """ Start the processing """
    # Assemble the job order
    job = Level2Job()
    job.local_machine_settings(config.local_machine)
    job.mission_settings(setting.mission)
    job.roi_settings(setting.roi)
    job.l2proc_settings(setting.level2)
    job.set_overwrite_protection(args.overwrite_protection_flag)
    job.validate()

    # Start the processor
    l2proc = Level2Processor(job)
    l2proc.set_config(config)
    l2proc.skip_files(args.skip)
    l2proc.error_handling(raise_on_error=True)

    iterator = month_iterator(args.start_date[0], args.start_date[1],
                              args.stop_date[0], args.stop_date[1])
    for year, month in iterator:
        l1bdata_files = get_l1bdata_files(
            job.mission.id, setting.roi.hemisphere, year, month, config=config)
        l2proc.set_l1b_files(l1bdata_files)
        l2proc.run()


def validate_year_month_list(year_month_list, label):
    try:
        datetime(year_month_list[0], year_month_list[1],  1)
    except ValueError:
        print "Error: Invalid "+label+" (%04g, %02g)" % (
            year_month_list[0], year_month_list[1])
        sys.exit(1)


def get_l2proc_argparser():
    """ Handle command line arguments """
    parser = argparse.ArgumentParser()

    # Mission id string: cryosat2, envisat, ...
    parser.add_argument(
        '-s', '-setting', action='store', dest='setting_id',
        help='setting id of yaml file in /settings/l2', required=True)

    # run tag:
    parser.add_argument(
        '-r', '-runtag', action='store', dest='run_tag',
        help='tag (str) for the processor run (directory in products folder)',
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
        "-skip", action='store', type=int, const=0, nargs='?',
        dest='skip', help='number of files to skip')

    # Show debug statements
    parser.add_argument(
        "-v", "--verbose", help="increase output verbosity",
        action="store_true")

    # Overwrite flag
    parser.add_argument('--overwrite-protection',
                        dest='overwrite_protection_flag',
                        action='store_true')
    parser.add_argument('--no-overwrite-protection',
                        dest='overwrite_protection_flag',
                        action='store_false')
    parser.set_defaults(overwrite=True)

    # show preprocessor version
    parser.add_argument(
        '--version', action='version', version='%(prog)s 0.1a')

    return parser


if __name__ == "__main__":
    pysiral_l2proc()
