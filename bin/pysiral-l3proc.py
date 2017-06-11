# -*- coding: utf-8 -*-

from pysiral.config import (ConfigInfo, DefaultCommandLineArguments,
                            TimeRangeRequest)
from pysiral.errorhandler import ErrorStatus
from pysiral.grid import GridDefinition
from pysiral.logging import DefaultLoggingClass
from pysiral.l3proc import Level3Processor, Level3ProductDefinition
from pysiral.output import DefaultLevel3OutputHandler

from datetime import datetime, timedelta
import argparse
import time
import sys
import os


def pysiral_l3proc():

    # parse command line arguments
    args = Level3ProcArgParser()
    args.parse_command_line_arguments()

    # Get start time of processor run
    t0 = time.clock()

    # Get monthly iterations from requested period
    period = TimeRangeRequest(
            args.start, args.stop,
            period=args.period)

    # Get the output grid
    grid = GridDefinition()
    grid.set_from_griddef_file(args.l3_griddef)

    # Initialize the output handler
    # Currently the overwrite protection is disabled per default
    output = DefaultLevel3OutputHandler(
            output_def=args.l3_output_file,
            base_directory=args.l3_product_basedir,
            overwrite_protection=False)

    # Compile the product def
    product_def = Level3ProductDefinition(
            args.l3_settings_file, grid, output)
    product_def.validate()

    # Initialize the Processor
    l3proc = Level3Processor(product_def)

    # Loop over all iterations
    for time_range in period.iterations:
        print time_range


    t1 = time.clock()
    seconds = int(t1-t0)
    l3proc.log.info("Run completed in %s" % str(timedelta(seconds=seconds)))

    # validate date values
    # XXX: This needs to be changed to TimeRangeRequest
#    validate_year_month_list(args.start_date, "start date")
#    validate_year_month_list(args.stop_date, "stop date")

    # Assemble the job order
    # This is actually not that bad, however needs to be extended with
    # an output handler
#    job = Level3Job()
#    job.set_input_directory(args.input_directory)
#    job.set_grid_definition(setting.grid_definition)
#    job.set_parameter(
#        l2=setting.l2_parameter, l3=setting.l3_parameter,
#        frb_nanmask=setting.freeboard_nan_mask_targets,
#        sic_mask=setting.sea_ice_concentration_mask_targets)
#    job.validate()
#
#    # Start the processor
#    l3proc = Level3Processor(job)
#
#    iterator = month_iterator(args.start_date[0], args.start_date[1],
#                              args.stop_date[0], args.stop_date[1])
#    for year, month in iterator:
#        l2idata_files = job.get_monthly_l2idata_files(year, month)
#        if len(l2idata_files) == 0:
#            continue
#        l3proc.set_l2_files(l2idata_files)
#        l3proc.run()


class Level3ProcArgParser(DefaultLoggingClass):

    def __init__(self):
        super(Level3ProcArgParser, self).__init__(self.__class__.__name__)
        self.error = ErrorStatus()
        self.pysiral_config = ConfigInfo()
        self._args = None

    def parse_command_line_arguments(self):
        # use python module argparse to parse the command line arguments
        # (first validation of required options and data types)
        self._args = self.parser.parse_args()

        # Add addtional check to make sure either `l1b-files` or
        # `start ` and `stop` are set
#        l1b_file_preset_is_set = self._args.l1b_files_preset is not None
#        start_and_stop_is_set = self._args.start_date is not None and \
#            self._args.stop_date is not None
#
#        if l1b_file_preset_is_set and start_and_stop_is_set:
#            self.parser.error("-start & -stop and -l1b-files are exclusive")
#
#        if not l1b_file_preset_is_set and not start_and_stop_is_set:
#            self.parser.error("either -start & -stop or -l1b-files required")

    def critical_prompt_confirmation(self):

        # Any confirmation prompts can be overriden by --no-critical-prompt
        no_prompt = self._args.no_critical_prompt

        # if --remove_old is set, all previous l1bdata files will be
        # erased for all month
        if self._args.remove_old and not no_prompt:
            message = "You have selected to remove all previous " + \
                "l3 files for the requested period\n" + \
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
            ("-l2i-product-dir", "l2i-product-dir", "l2i_basedir", True),
            ("-l3-settings", "l3-settings", "l3_settings", False),
            ("-l3-griddef", "l3-griddef", "l3_griddef", True),
            ("-l3-output", "l3-output", "l3_output", True),
            ("-start", "date", "start_date", True),
            ("-stop", "date", "stop_date", True),
            ("-period", "period", "period", False),
            ("--remove-old", "remove-old", "remove_old", False),
            ("--no-critical-prompt", "no-critical-prompt",
             "no_critical_prompt", False)]

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
        return self._args.__dict__

    @property
    def start(self):
        return self._args.start_date

    @property
    def stop(self):
        return self._args.stop_date

    @property
    def period(self):
        return self._args.period

    @property
    def l3_settings_file(self):
        l3_settings = self._args.l3_settings
        filename = self.pysiral_config.get_settings_file("l3", l3_settings)
        if filename is None:
            msg = "Invalid l3 settings filename or id: %s\n" % l3_settings
            msg = msg + " \nRecognized Level-2 processor setting ids:\n"
            for l3_settings_id in self.pysiral_config.get_setting_ids("l3"):
                msg = msg + "  " + l3_settings_id+"\n"
            self.error.add_error("invalid-l3-settings", msg)
            self.error.raise_on_error()
        else:
            return filename

    @property
    def l3_griddef(self):
        l3_griddef = self._args.l3_griddef
        filename = self.pysiral_config.get_settings_file("griddef", l3_griddef)
        if filename is None:
            msg = "Invalid griddef filename or id: %s\n" % l3_griddef
            msg = msg + " \nRecognized grid definition ids:\n"
            for griddef_id in self.pysiral_config.get_setting_ids("griddef"):
                msg = msg + "  " + griddef_id+"\n"
            self.error.add_error("invalid-griddef", msg)
            self.error.raise_on_error()
        else:
            return filename

    @property
    def l3_output_file(self):
        l3_output = self._args.l3_output
        filename = self.pysiral_config.get_settings_file(
                "outputdef", l3_output)
        if filename is None:
            msg = "Invalid griddef filename or id: %s\n" % l3_output
            msg = msg + " \nRecognized grid definition ids:\n"
            for output_id in self.pysiral_config.get_setting_ids("outputdef"):
                msg = msg + "  " + output_id+"\n"
            self.error.add_error("invalid-outputdefdef", msg)
            self.error.raise_on_error()
        else:
            return filename

    @property
    def l3_product_basedir(self):
        """ Returns the base directory (one level below l2i) """
        # 1. Clean up the path
        product_basedir = os.path.abspath(self._args.l2i_basedir)
        dirs = os.path.split(product_basedir)
        if dirs[1] == "l2i":
            return dirs[0]
        else:
            return product_basedir

    @property
    def remove_old(self):
        return self._args.remove_old and not self._args.overwrite_protection

#def validate_year_month_list(year_month_list, label):
#    try:
#        datetime(year_month_list[0], year_month_list[1],  1)
#    except ValueError:
#        print "Error: Invalid "+label+" (%04g, %02g)" % (
#            year_month_list[0], year_month_list[1])
#        sys.exit(1)
#
#
#def get_l1bdata_files(config, job, hemisphere, year, month):
#    l1b_repo = config.local_machine.l1b_repository[job.mission.id].l1bdata
#    directory = os.path.join(
#        l1b_repo, hemisphere, "%04g" % year, "%02g" % month)
#    print directory
#    l1bdata_files = sorted(glob.glob(os.path.join(directory, "*.nc")))
#    return l1bdata_files
#
#
#def get_l3proc_argparser():
#    """ Handle command line arguments """
#    parser = argparse.ArgumentParser()
#    # Mission id string: cryosat2, envisat, ...
#    parser.add_argument(
#        '-s', '-setting', action='store', dest='setting_id',
#        help='setting id of yaml file in /settings/l3', required=True)
#    # input directory
#    parser.add_argument(
#        '-i', action='store', dest='input_directory',
#        help='link to directory with l2 netcdf files',
#        required=True)
#    # Start month as list: [yyyy, mm]
#    parser.add_argument(
#        '-start', action='store', dest='start_date',
#        nargs='+', type=int, required=True,
#        help='start date as year and month (-start yyyy mm)')
#    # Stop month as list: [yyyy, mm]
#    parser.add_argument(
#        '-stop', action='store', dest='stop_date',
#        nargs='+', type=int,  required=True,
#        help='start date as year and month (-stop yyyy mm)')
#    # Add an Option to skip a number of files (e.g. for a restart)
#    parser.add_argument(
#        "-skipmonth", action='store', type=int, const=[], nargs='?',
#        dest='skip_month', help='month to skip (default: none)')
#    # Show debug statements
#    parser.add_argument(
#        "-v", "--verbose", help="increase output verbosity",
#        action="store_true")
#    # show preprocessor version
#    parser.add_argument(
#        '--version', action='version', version='%(prog)s 0.1a')
#
#    return parser


if __name__ == "__main__":
    pysiral_l3proc()
