# -*- coding: utf-8 -*-

from pysiral.config import DefaultCommandLineArguments
# from pysiral.config import ConfigInfo, get_yaml_config
# from pysiral.helper import month_iterator

from pysiral.l2proc import Level2Processor, L2ProcJob
from pysiral.iotools import get_l1bdata_files

import argparse
import sys


def pysiral_l2proc():

    # Collect job settings from pysiral configuration data and
    # command line arguments
    args = L2ProcArgParser()

    # Parse and validate the command line arguments
    args.parse_command_line_arguments()

    # Get confirmation for critical choices (if necessary)
    args.critical_prompt_confirmation()

    # Set up the job definition for the Level-2 processor
    jobdef = L2ProcJob()

    # configure the job definition
    jobdef.options.from_dict(args.arg_dict)

    # Parse the Level-2 settings file
    jobdef.parse_l2_settings()
    jobdef.error.raise_on_error()

    # There are several options of how the time range can be requested
    # a) start and stop month (no days specified)
    # b) full start and stop date
    # In addition, full month can be excluded from the pre-processing
    # (e.g. summer Arctic)
    jobdef.process_requested_time_range()

    # Validate given infos and availability of data files
    jobdef.validate()
    jobdef.error.raise_on_error()

    # L1b data is grouped by month
    # -> if requested time range exceeds one month, the processor
    #    will run in several iterations
    jobdef.generate_preprocessor_iterations()

    # Initialize the Level-2 Processor
    job = Level2Processor(jobdef)
    job.initialize()

    # Loop over iterations (one per month)
    for time_range in jobdef.iterations:

        job.log.info("Processing period: %s" % time_range.label)

#        # Get the list of input files from local machine def
#        version = jobdef.input_version
#        job.get_input_files_local_machine_def(time_range, version=version)
#
#        # Check error, e.g. problems with local_machine_def, ...
#        job.error.raise_on_error()
#
#        # List might be empty
#        if job.has_empty_file_list:
#            job.log.info(" Skipping period: %s" % time_range.label)
#            continue
#
#        # Empty output folder (if --remove_old is set)
#        if jobdef.remove_old and not joddef.overwrite_protection:
#            job.remove_old_l2data()

        # Pre-process data for one month
        # job.run()


class L2ProcArgParser(object):

    def __init__(self):
        self.args = None

    def parse_command_line_arguments(self):
        # use python module argparse to parse the command line arguments
        # (first validation of required options and data types)
        self.args = self.parser.parse_args()

    def critical_prompt_confirmation(self):

        # Any confirmation prompts can be overriden by --no-critical-prompt
        no_prompt = self.args.no_critical_prompt

        # if --remove_old is set, all previous l1bdata files will be
        # erased for all month
        if self.args.remove_old and not no_prompt:
            message = "You have selected to remove all previous " + \
                "l2 files for the requested period\n" + \
                "(Note: use --no-critical-prompt to skip confirmation)\n" + \
                "Enter \"YES\" to confirm and continue: "
            result = raw_input(message)

            if result != "YES":
                sys.exit(1)

    @property
    def parser(self):
        # XXX: Move back to caller

        # Take the command line options from default settings
        # -> see config module for data types, destination variables, etc.
        clargs = DefaultCommandLineArguments()

        # List of command line option required for pre-processor
        # (argname, argtype (see config module), destination, required flag)
        options = [
            ("-l2-settings", "l2-settings", "l2_settings", True),
            ("-run-tag", "run-tag", "run_tag", True),
            ("-start", "date", "start_date", True),
            ("-stop", "date", "stop_date", True),
            ("-exclude-month", "exclude-month", "exclude_month", False),
            ("-input-version", "input-version", "input_version", False),
            ("--remove-old", "remove-old", "remove_old", False),
            ("--no-critical-prompt", "no-critical-prompt",
             "no_critical_prompt", False),
            ("--no-overwrite-protection", "no-overwrite-protection",
             "overwrite_protection", False),
            ("--overwrite-protection", "overwrite-protection",
             "overwrite_protection", False)]

        # create the parser
        parser = argparse.ArgumentParser()
        for option in options:
            argname, argtype, destination, required = option
            argparse_dict = clargs.get_argparse_dict(
                argtype, destination, required)
            parser.add_argument(argname, **argparse_dict)

        parser.set_defaults(overwrite_protection=True)

        return parser

    @property
    def arg_dict(self):
        """ Return the arguments as dictionary """
        return self.args.__dict__



#    """ get the pysiral configuration info """
#    config = ConfigInfo()
#
#    """ parse command line arguments """
#    parser = get_l2proc_argparser()
#    args = parser.parse_args()
#
#    """ Read the settings file """
#    # is filename
#    if os.path.isfile(args.setting_id):
#        setting_file = args.setting_id
#    # if not filename, than it need to be id of settings file in
#    # pysiral\config\l2
#    else:
#        setting_file = os.path.join(
#            config.pysiral_local_path, "settings", "l2",
#            args.setting_id+".yaml")
#        if not os.path.isfile(setting_file):
#            error_message = "Unknown l2 settings file: %s" % setting_file
#            sys.exit(error_message)
#    setting = get_yaml_config(setting_file)
#
#    """ Add run tag to settings """
#    setting.level2.run_tag = args.run_tag
#
##    """ validate mission id """
##    if args.mission_id not in config.mission_ids:
##        print "Error: Invalid mission id (%s)" % args.mission_id
##        print "use: " + " | ".join(config.mission_ids)
##        sys.exit(1)
#
#    """ validate date values """
#    validate_year_month_list(args.start_date, "start date")
#    validate_year_month_list(args.stop_date, "stop date")
#
#    """ Check if data are out of bounds of mission lifetime """
#    mission_info = config.get_mission_info(setting.mission.id)
#    start_date = datetime(args.start_date[0], args.start_date[1], 1)
#    stop_date = datetime(args.stop_date[0], args.stop_date[1], 1)
#    # Set start time to mission start time (if necessary)
#    if start_date < mission_info.data_period.start:
#        print "Warning: start date before %s data period %s" % (
#            args.mission_id, str(mission_info.data_period.start))
#        args.start_date[0] = mission_info.data_period.start.year
#        args.start_date[1] = mission_info.data_period.start.month
#    # Stop time can be None if mission is ongoing
#    if mission_info.data_period.stop is None:
#        mission_stop_time = datetime.utcnow()
#    else:
#        mission_stop_time = mission_info.data_period.stop
#    if stop_date > mission_stop_time:
#        print "Warning: stopt date late then %s data period %s" % (
#            args.mission_id, str(mission_stop_time))
#        args.stop_date[0] = mission_stop_time.data_period.start.year
#        args.stop_date[1] = mission_stop_time.data_period.start.month
#
#    """ Start the processing """
#    # Assemble the job order
#    job = Level2Job()
#    job.local_machine_settings(config.local_machine)
#    job.mission_settings(setting.mission)
#    job.roi_settings(setting.roi)
#    job.l2proc_settings(setting.level2)
#    job.set_overwrite_protection(args.overwrite_protection_flag)
#    job.validate()
#
#    # Start the processor
#    l2proc = Level2Processor(job)
#    l2proc.set_config(config)
#    l2proc.skip_files(args.skip)
#    l2proc.error_handling(raise_on_error=True)
#
#    iterator = month_iterator(args.start_date[0], args.start_date[1],
#                              args.stop_date[0], args.stop_date[1])
#    for year, month in iterator:
#        l1bdata_files = get_l1bdata_files(
#            job.mission.id, setting.roi.hemisphere, year, month, config=config)
#        l2proc.set_l1b_files(l1bdata_files)
#        l2proc.run()


#def validate_year_month_list(year_month_list, label):
#    try:
#        datetime(year_month_list[0], year_month_list[1],  1)
#    except ValueError:
#        print "Error: Invalid "+label+" (%04g, %02g)" % (
#            year_month_list[0], year_month_list[1])
#        sys.exit(1)


#def get_l2proc_argparser():
#    """ Handle command line arguments """
#    parser = argparse.ArgumentParser()
#
#    # Mission id string: cryosat2, envisat, ...
#    parser.add_argument(
#        '-s', '-setting', action='store', dest='setting_id',
#        help='setting id of yaml file in /settings/l2', required=True)
#
#    # run tag:
#    parser.add_argument(
#        '-r', '-runtag', action='store', dest='run_tag',
#        help='tag (str) for the processor run (directory in products folder)',
#        required=True)
#
#    # Start month as list: [yyyy, mm]
#    parser.add_argument(
#        '-start', action='store', dest='start_date',
#        nargs='+', type=int, required=True,
#        help='start date as year and month (-start yyyy mm)')
#
#    # Stop month as list: [yyyy, mm]
#    parser.add_argument(
#        '-stop', action='store', dest='stop_date',
#        nargs='+', type=int,  required=True,
#        help='start date as year and month (-stop yyyy mm)')
#
#    # Add an Option to skip a number of files (e.g. for a restart)
#    parser.add_argument(
#        "-skip", action='store', type=int, const=0, nargs='?',
#        dest='skip', help='number of files to skip')
#
#    # Show debug statements
#    parser.add_argument(
#        "-v", "--verbose", help="increase output verbosity",
#        action="store_true")
#
#    # Overwrite flag
#    parser.add_argument('--overwrite-protection',
#                        dest='overwrite_protection_flag',
#                        action='store_true')
#    parser.add_argument('--no-overwrite-protection',
#                        dest='overwrite_protection_flag',
#                        action='store_false')
#    parser.set_defaults(overwrite=True)
#
#    # show preprocessor version
#    parser.add_argument(
#        '--version', action='version', version='%(prog)s 0.1a')
#
#    return parser




if __name__ == "__main__":
    pysiral_l2proc()
