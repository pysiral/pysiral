# -*- coding: utf-8 -*-

"""
Main entry point for pysiral console scripts.
"""

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"

import sys
from typing import List

# Ensure that the logger is initialized
import pysiral._logger

# from pysiral.scripts.info import info, InfoScriptArguments
from pysiral.scripts.l1preproc import l1preproc, L1PreProcScriptArguments
from pysiral.scripts.l2proc import l2proc, L2ProcScriptArguments
from pysiral.scripts.l2procfiles import l2procfiles, L2ProcFilesScriptArguments
from pysiral.scripts.l2preproc import l2preproc, L2PreProcScriptArguments
from pysiral.scripts.l3proc import l3proc, L3ProcScriptArguments


def main() -> None:

    usage_str = """
    pysiral - python sea ice radar altimetry processing toolbox
    -----------------------------------------------------------
    
        Usage: pysiral(.exe) "{target_script} {script_cli_args}"

    Valid target scripts are:

        `info`        Print information about the pysiral installation.
                      (see: pysiral info --help)
                      
        `config`      Set pysiral configuration to a specific directory.
                      (see: pysiral set-cfg --help)

        `l1preproc`   Generate Level-1 files (l1p) with trajectory sensors data
                      from source files.  
                      (see: pysiral l1preproc --help)

        `l2proc`      Generate Level-2 files (l2/l2i) with geophysical information
                      from Level-1 files (l1p) and auxiliary data for a given 
                      Level-2 product definition and period. 
                      (see: pysiral l2proc --help)
                      
        `l2procfiles` Generate Level-2 files (l2/l2i) with geophysical information
                      from a list of Level-1 files (l1p) and auxiliary data for a given 
                      Level-2 product definition. 
                      (see: pysiral l2procfiles --help)                      

        `l2preproc`   Generate Level-2 files (l2p) wwith daily summaries of l2/l2i files
                      (see: pysiral l2preproc --help)

        `l3proc`      Generate Level-3 files (l3c/l3s) with gridded data from 
                      Level-2 files (l2/l2i).
                      (see: pysiral l3proc --help)
    """

    # Get the name of the script
    try:
        target_script = sys.argv[1]
    except IndexError:
        print(usage_str)
        sys.exit(1)

    # Get the script CLI arguments (if any) and run script
    args = sys.argv[2:] if len(sys.argv) > 2 else []
    try:
        func = globals()[f"{target_script}_cli"]
    except KeyError:
        print(usage_str)
        print(f"[Error] No such script: {target_script}")
        sys.exit(1)
    func(args)


# def info_cli(args_list: List = None) -> None:
#     """
#     Command-line interface entry point for the `pysiral info` script.
#
#     :param args_list: Command line arguments to be passed to the script.
#
#     :return: None
#     """
#     info(**vars(InfoScriptArguments().get(args_list)))


def l1preproc_cli(args_list: List = None) -> None:
    """
    Command-line interface entry point for the `pysiral l1preproc` script.

    :param args_list: Command line arguments to be passed to the script.

    :return: None
    """
    l1preproc(**vars(L1PreProcScriptArguments().get(args_list)))


def l2proc_cli(args_list: List = None) -> None:
    """
    Command-line interface entry point for the `pysiral l2proc` script.

    :param args_list: Command line arguments to be passed to the script.

    :return: None
    """
    l2proc(**vars(L2ProcScriptArguments().get(args_list)))


def l2preproc_cli(args_list: List = None) -> None:
    """
    Command-line interface entry point for the `pysiral l2preproc` script.

    :param args_list: Command line arguments to be passed to the script.

    :return: None
    """
    l2preproc(**vars(L2PreProcScriptArguments().get(args_list)))


def l3proc_cli(args_list: List = None) -> None:
    """
    Command-line interface entry point for the `pysiral l3proc` script.

    :param args_list: Command line arguments to be passed to the script.

    :return: None
    """
    l3proc(**vars(L3ProcScriptArguments().get(args_list)))


if __name__ == "__main__":
    main()
