# -*- coding: utf-8 -*-

from pysiral.config import (ConfigInfo, TimeRangeRequest,
                            DefaultCommandLineArguments)
from pysiral.logging import DefaultLoggingClass
from pysiral.helper import month_iterator

from datetime import datetime
import argparse
import sys
import os


def pysiral_l1bpreproc():
    """ Main call for l1b preprocessing """

    # get the pysiral configuration info
    config = ConfigInfo()

    # Collect job settings from pysiral configuration data and
    # command line arguments
    jobdef = L1bPreProcJobSettings(config)

    # Parse and validate the command line arguments
    jobdef.parse_command_line_arguments()

    # There are several options of how the time range can be requested
    # a) start and stop month (no days specified)
    # b) full start and stop date
    # In addition, full month can be excluded from the pre-processing
    # (e.g. summer Arctic)
    jobdef.process_requested_time_range()

    # L1b data is grouped by month
    # -> if requested time range exceeds one month, the preprocessor
    #    will run in several iterations
    jobdef.generate_preprocessor_iterations()

    # Get the mission-specific preprocessor class
    preprocessor = jobdef.get_mission_preprocessor()

    # Loop over iterations (one per month)
    for iteration in jobdef.iterations:

        # Get the list of input files
        #
        file_list = jobdef.get_file_list(iteration)

        #
        if len(file_list) == 0:
            continue

        job = preprocessor()

        # Though
        job.mission = jobdef.mission_id
        job.file_list = file_list
        job.execute()


class L1bPreProcJobSettings(DefaultLoggingClass):

    def __init__(self, config):

        super(L1bPreProcJobSettings, self).__init__(self.__class__.__name__)

        # Save pointer to pysiral configuration
        self.pysiral_config = config

        # Initialize the time range and set to monthly per default
        self.time_range = TimeRangeRequest()
        self.time_range.set_period("monthly")

        # Initialize job parameter
        self.args = None
        self.iterations = []

    def parse_command_line_arguments(self):

        # use python module argparse to parse the command line arguments
        # (first validation of required options and data types)
        self.args = self.parser.parse_args()

    def generate_preprocessor_iterations(self):
        pass

    def process_requested_time_range(self):

        config = self.pysiral_config
        mission_info = config.get_mission_info(self.args.mission_id)

        # Set and validate the time range
        start_date, stop_date = self.args.start_date, self.args.stop_date
        self.time_range.set_range(start_date, stop_date)

        # Check if any errors in definitions
        self.time_range.error.raise_on_error()

        # Clip time range to mission time range
        self.time_range.clip_to_range(mission_info.data_period.start,
                                      mission_info.data_period.stop)
#
#        # start stop mus be a valid year/month combination
#        validate_year_month_list(args.start_date, "start date")
#        validate_year_month_list(args.stop_date, "stop date")
#
#        # Dates should not exceed mission lifetime
#        start_date = datetime(args.start_date[0], args.start_date[1], 1)
#        stop_date = datetime(args.stop_date[0], args.stop_date[1], 1)
#
#        # Set start time to mission start time (if necessary)
#        if start_date < mission_info.data_period.start:
#            print "Warning: start date before %s data period %s" % (
#                args.mission_id, str(mission_info.data_period.start))
#            args.start_date[0] = mission_info.data_period.start.year
#            args.start_date[1] = mission_info.data_period.start.month
#        # Stop time can be None if mission is ongoing
#        if mission_info.data_period.stop is None:
#            mission_stop_time = datetime.utcnow()
#        else:
#            mission_stop_time = mission_info.data_period.stop
#        if stop_date > mission_stop_time:
#            print "Warning: stopt date late then %s data period %s" % (
#                args.mission_id, str(mission_stop_time))
#            args.stop_date[0] = mission_stop_time.data_period.start.year
#            args.stop_date[1] = mission_stop_time.data_period.start.month

    def get_mission_preprocessor(self):

        from pysiral.cryosat2.preproc import CryoSat2PreProc
        from pysiral.envisat.preproc import EnvisatPreProc
        from pysiral.ers.preproc import ERSPreProc
        from pysiral.sentinel3.preproc import Sentinel3PreProc

        if self.mission_id == "cryosat2":
            return CryoSat2PreProc
        elif self.mission_id == "envisat":
            return EnvisatPreProc
        elif self.mission_id == "ers2":
            return ERSPreProc
        elif self.mission_id == "sentinel3a":
            return Sentinel3PreProc
        else:
            print "error: mission %s currently not supported" % self.mission_id
            sys.exit(1)

    @property
    def n_l1b_files(self):
        return len(self.l1b_files)

    @property
    def parser(self):

        # Take the command line options from default settings
        # -> see config module for data types, destination variables, etc.
        clargs = DefaultCommandLineArguments()

        # List of command line option required for pre-processor
        # (argname, argtype (see config module), destination, required flag)
        options = [
            ("-mission", "mission", "mission_id", True),
            ("-start", "date", "start_date", True),
            ("-stop", "date", "stop_date", True),
            ("-hemisphere", "hemisphere", "hemisphere", False),
            ("-exclude-month", "exclude-month", "exclude-month", False)]

        # create the parser
        parser = argparse.ArgumentParser()
        for option in options:
            argname, argtype, destination, required = option
            argparse_dict = clargs.get_argparse_dict(
                argtype, destination, required)
            parser.add_argument(argname, **argparse_dict)

        return parser

    @property
    def mission_id(self):
        return self.args.mission_id

    @property
    def month_list(self):
        iterator = month_iterator(
            self.args.start_date[0], self.args.start_date[1],
            self.args.stop_date[0], self.args.stop_date[1])
        return iterator


def validate_year_month_list(year_month_list, label):
    try:
        datetime(year_month_list[0], year_month_list[1], 1)
    except ValueError:
        print "Error: Invalid "+label+" (%04g, %02g)" % (
            year_month_list[0], year_month_list[1])
        sys.exit(1)






if __name__ == "__main__":
    pysiral_l1bpreproc()
