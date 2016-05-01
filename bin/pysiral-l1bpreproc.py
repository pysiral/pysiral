# -*- coding: utf-8 -*-

from pysiral.config import ConfigInfo
from pysiral.helper import month_iterator

from datetime import datetime
import argparse
import sys


def pysiral_l1bpreproc():

    """ get the pysiral configuration info """
    config = ConfigInfo()

    """ parse command line arguments """
    parser = get_l1bpreproc_argparser()
    args = parser.parse_args()

    """ validate mission id """
    if args.mission_id not in config.mission_ids:
        print "Error: Invalid mission id (%s)" % args.mission_id
        print "use: " + " | ".join(config.mission_ids)
        sys.exit(1)

    """ validate date values """
    validate_year_month_list(args.start_date, "start date")
    validate_year_month_list(args.stop_date, "stop date")

    """ Check if data are out of bounds of mission lifetime """
    mission_info = config.get_mission_info(args.mission_id)
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
    preprocessor = get_mission_preprocessor(args.mission_id)
    iterator = month_iterator(args.start_date[0], args.start_date[1],
                              args.stop_date[0], args.stop_date[1])

    for year, month in iterator:
        job = preprocessor()
        job.mission = args.mission_id
        job.year = year
        job.month = month
        job.skip = args.skip
        job.execute()


def validate_year_month_list(year_month_list, label):
    try:
        datetime(year_month_list[0], year_month_list[1], 1)
    except ValueError:
        print "Error: Invalid "+label+" (%04g, %02g)" % (
            year_month_list[0], year_month_list[1])
        sys.exit(1)


def get_l1bpreproc_argparser():
    """ Handle command line arguments """
    parser = argparse.ArgumentParser()
    # Mission id string: cryosat2, envisat, ...
    parser.add_argument(
        '-m', action='store', dest='mission_id',
        help='pysiral recognized mission id', required=True)
    # Start month as list: [yyyy, mm]
    parser.add_argument(
        '-start', action='store', dest='start_date',
        nargs='+', type=int, required=True,
        help='start date as year and month (-t0 yyyy mm)')
    # Stop month as list: [yyyy, mm]
    parser.add_argument(
        '-stop', action='store', dest='stop_date',
        nargs='+', type=int,  required=True,
        help='start date as year and month (-t0 yyyy mm)')
    # Show debug statements
    parser.add_argument(
        "-v", "--verbose", help="increase output verbosity",
        action="store_true")
    # Add an Option to skip a number of files (e.g. for a restart)
    parser.add_argument(
        "-s", "-skip", action='store', type=int, const=0, nargs='?',
        dest='skip', help='number of files to skip')
    # show preprocessor version
    parser.add_argument(
        '--version', action='version', version='%(prog)s 0.1a')

    return parser


def get_mission_preprocessor(mission_id):

    from pysiral.cryosat2.preproc import CryoSat2PreProc
    from pysiral.envisat.preproc import EnvisatPreProc
    from pysiral.ers.preproc import ERSPreProc

    if mission_id == "cryosat2":
        return CryoSat2PreProc
    elif mission_id == "envisat":
        return EnvisatPreProc
    elif mission_id == "ers2":
        return ERSPreProc
    else:
        print "error: mission %s currently not supported" % mission_id
        sys.exit(1)


if __name__ == "__main__":
    pysiral_l1bpreproc()
