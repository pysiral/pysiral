# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 18:54:15 2016

@author: Stefan

Testbed for CryoSat-2 specific pre-processor that
- creates l1bncfiles from the l1b source files
- concatenates adjacent files (e.g. at borders of SAR and SIN modes)
- trims to ocean data in Arctic and Antarctic

Filenaming:
    l1bdat_v$VERS_$REG_$MISSION_$ORBIT_$YYYYMMDDHHMISS_$YYYYMMDDHHMISS.nc

$VERS            l1bdata version       00 (beta)
$REG             region                [north | south]
$MISSION         mission short name    => cryosat2
$ORBIT           orbit/cylce number    e.g. 00026000
$YYYYMMDDHHMISS  start and end time

Definition of Arctic:
    latitude >= 50.0

Definition of Antarctic:
    latitude <= -50.0


Concatinating orbit files (CryoSat-2 specific):

- Loop over L1b files (month wise, sar and sin)
- parse l1b source files
- check if ocean content
- start a stack of l1bdata objectes
- check of adjacent to previous file
    yes: add to l1bdata stack
    no: - concatenate files in stack
        - write to l1bdata nc file
        - clear stack
        - reopen stack with current file

file where strange things are happening:

over southern ocean but no has_ocean => false
D:\awi\altim\data\altimetry\cryosat2\baseline-c\SIR_SAR_L1\2015\03\...
 ...CS_LTA__SIR_SAR_1B_20150301T025222_20150301T025526_C001.DBL

potential wrong lat/lon:
D:\awi\altim\data\altimetry\cryosat2\baseline-c\SIR_SAR_L1\2015\03\...
 ...CS_LTA__SIR_SAR_1B_20150301T032811_20150301T033251_C001.DBL

"""

from pysiral.config import ConfigInfo
from pysiral.logging import stdout_logger
from pysiral.cryosat2.iotools import CryoSat2FileListAllModes
from pysiral.cryosat2.preproc import CryoSat2PreProcJob


def preproc_cryosat2():

    """ Get the logging instance """
    log = stdout_logger("cryosat2-preproc")

    """ Get the configuration data for handling CryoSat-2 data """
    config = ConfigInfo()
    cryosat_l1b_repository = config.local_machine.l1b_repository.cryosat2

    """ Get the list of files (SAR and SIN in chronological order) """
    # for the case of this test only for one month of data
    cryosat2_files = CryoSat2FileListAllModes()
    cryosat2_files.log = log
    cryosat2_files.folder_sar = cryosat_l1b_repository.sar
    cryosat2_files.folder_sin = cryosat_l1b_repository.sin
    cryosat2_files.year = 2015
    cryosat2_files.month = 3
    cryosat2_files.search()

    """ Start the CryoSat-2 pre-processor """
    job = CryoSat2PreProcJob()
    job.config = config
    job.log = log
    job.files = cryosat2_files.sorted_list
    job.merge_and_export_polar_ocean_subsets()


if __name__ == "__main__":
    preproc_cryosat2()
