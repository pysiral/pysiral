# -*- coding: utf-8 -*-

# from pysiral.config import (TimeRangeRequest, DefaultCommandLineArguments)
from pysiral.config import DefaultCommandLineArguments
# from pysiral.errorhandler import ErrorStatus
from pysiral.l1bpreproc import L1bPreProcJob
# from pysiral.logging import DefaultLoggingClass

import argparse


def pysiral_l1bpreproc():
    """ Main call for l1b preprocessing """

    # Collect job settings from pysiral configuration data and
    # command line arguments
    args = L1bPreProcArgParser()

    # Parse and validate the command line arguments
    args.parse_command_line_arguments()

    # Set up the pre-processing job
    jobdef = L1bPreProcJob()

    # configure the job definition
    jobdef.options.from_dict(args.arg_dict)

    # There are several options of how the time range can be requested
    # a) start and stop month (no days specified)
    # b) full start and stop date
    # In addition, full month can be excluded from the pre-processing
    # (e.g. summer Arctic)
    jobdef.process_requested_time_range()

    # Check if everything is ok
    jobdef.error.raise_on_error()

    # L1b data is grouped by month
    # -> if requested time range exceeds one month, the preprocessor
    #    will run in several iterations
    jobdef.generate_preprocessor_iterations()

    # Get the mission-specific preprocessor class
    preprocessor = jobdef.get_mission_preprocessor()

    # Loop over iterations (one per month)
    for time_range in jobdef.iterations:

        # Start the pre-processing
        job = preprocessor()
        job.log.info("Processing period: %s" % time_range.label)

        # Set the pre-processing parameter
        job.set_job_definition(jobdef)

        # Get the list of input files from local machine def
        version = jobdef.input_version
        job.set_input_files_local_machine_def(time_range, version=version)

        # Check error, e.g. problems with local_machine_def, ...
        job.error.raise_on_error()

        # List might be empty
        if job.has_empty_file_list:
            job.log.info(" Skipping period: %s" % time_range.label)
            continue

        # Empty output folder (if --remove_old is set)
        if jobdef.remove_old:
            job.remove_old_l1bdata()

        # Pre-process data for one month
        job.execute()


class L1bPreProcArgParser(object):

    def __init__(self):
        self.args = None

    def parse_command_line_arguments(self):
        # use python module argparse to parse the command line arguments
        # (first validation of required options and data types)
        self.args = self.parser.parse_args()

    @property
    def parser(self):
        # XXX: Move back to caller

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
            ("-exclude-month", "exclude-month", "exclude_month", False),
            ("-input-version", "input-version", "input_version", False),
            ("--remove-old", "remove-old", "remove_old", False)]

        # create the parser
        parser = argparse.ArgumentParser()
        for option in options:
            argname, argtype, destination, required = option
            argparse_dict = clargs.get_argparse_dict(
                argtype, destination, required)
            parser.add_argument(argname, **argparse_dict)

        return parser

    @property
    def arg_dict(self):
        """ Return the arguments as dictionary """
        return self.args.__dict__


#class L1bPreProcJobSettings(DefaultLoggingClass):
#
#    def __init__(self, config):
#
#        super(L1bPreProcJobSettings, self).__init__(self.__class__.__name__)
#
#        # Save pointer to pysiral configuration
#        self.pysiral_config = config
#
#        # Initialize the time range and set to monthly per default
#        self.time_range = TimeRangeRequest()
#
#        # Error Status
#        self.error = ErrorStatus()
#
#        # Initialize job parameter
#        self.args = None
#        self.iterations = []
#
#    def parse_command_line_arguments(self):
#
#        # use python module argparse to parse the command line arguments
#        # (first validation of required options and data types)
#        self.args = self.parser.parse_args()
#
#    def generate_preprocessor_iterations(self):
#
#        # The input data is organized in folder per month, therefore
#        # the processing period is set accordingly
#        self.time_range.set_period("monthly")
#        self.log.info("Pre-processor Period is monthly")
#        self.time_range.set_exclude_month(self.args.exclude_month)
#        self.log.info("Excluding month: %s" % str(self.args.exclude_month))
#        self.iterations = self.time_range.get_iterations()
#        self.log.info("Number of iterations: %g" % len(self.iterations))
#
#    def process_requested_time_range(self):
#
#        config = self.pysiral_config
#        mission_info = config.get_mission_info(self.args.mission_id)
#
#        # Set and validate the time range
#        start_date, stop_date = self.args.start_date, self.args.stop_date
#
#        self.time_range.set_range(start_date, stop_date)
#
#        # Check if any errors in definitions
#        self.time_range.error.raise_on_error()
#
#        self.log.info("Requested time range is: %s" % self.time_range.label)
#
#        # Clip time range to mission time range
#        is_clipped = self.time_range.clip_to_range(
#            mission_info.data_period.start, mission_info.data_period.stop)
#
#        if is_clipped:
#            self.log.info("Clipped to mission time range: %s till %s" % (
#                mission_info.data_period.start, mission_info.data_period.stop))
#            # Check if range is still valid
#            self.time_range.raise_if_empty()
#
#    def get_mission_preprocessor(self):
#
#        from pysiral.cryosat2.preproc import CryoSat2PreProc
#        from pysiral.envisat.preproc import EnvisatPreProc
#        from pysiral.ers.preproc import ERSPreProc
#        from pysiral.sentinel3.preproc import Sentinel3PreProc
#
#        if self.mission_id == "cryosat2":
#            return CryoSat2PreProc
#        elif self.mission_id == "envisat":
#            return EnvisatPreProc
#        elif self.mission_id == "ers2":
#            return ERSPreProc
#        elif self.mission_id == "sentinel3a":
#            return Sentinel3PreProc
#        else:
#            error_code = self.__class__.__name__+" (1)"
#            error_message = "Invalid mission_id: %s" % self.mission_id
#            self.error.add_error(error_code, error_message)
#
#    @property
#    def parser(self):
#
#        # Take the command line options from default settings
#        # -> see config module for data types, destination variables, etc.
#        clargs = DefaultCommandLineArguments()
#
#        # List of command line option required for pre-processor
#        # (argname, argtype (see config module), destination, required flag)
#        options = [
#            ("-mission", "mission", "mission_id", True),
#            ("-start", "date", "start_date", True),
#            ("-stop", "date", "stop_date", True),
#            ("-hemisphere", "hemisphere", "hemisphere", False),
#            ("-exclude-month", "exclude-month", "exclude_month", False),
#            ("-input-version", "input-version", "input_version", False),
#            ("--remove-old", "remove-old", "remove_old_l1bdata", False)]
#
#        # create the parser
#        parser = argparse.ArgumentParser()
#        for option in options:
#            argname, argtype, destination, required = option
#            argparse_dict = clargs.get_argparse_dict(
#                argtype, destination, required)
#            parser.add_argument(argname, **argparse_dict)
#
#        return parser
#
#    @property
#    def mission_id(self):
#        return self.args.mission_id
#
#    @property
#    def hemisphere(self):
#        return self.args.hemisphere
#
#    @property
#    def input_version(self):
#        return self.args.input_version
#
#    @property
#    def remove_old(self):
#        return self.args.remove_old_l1bdata

if __name__ == "__main__":
    pysiral_l1bpreproc()
