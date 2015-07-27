# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 14:09:32 2015

@author: Stefan
"""

from pysiral.config import ConfigInfo
from pysiral.job import Level2Job
from pysiral.l2proc import Level2Processor

import os
import glob


def l2_processing():

    # Get configuration
    config = ConfigInfo()

    # Get an L1B SAR file
    l1b_directory = config.local_machine.local_l1b_repository.cryosat2.sar
    l1b_directory = os.path.join(l1b_directory, "2015", "04")
    l1b_files = glob.glob(os.path.join(l1b_directory, "*.DBL"))

    # Simulate the jobconfig
    # This has to come from the job configuration file
    mission_settings = {
        "id": "cryosat2",
        "options": config.get_mission_defaults("cryosat2")}
    roi_settings = {
        "pyclass": "LowerLatLimit",
        "options": {
            "latitude_threshold": -50.0}}
    l2_settings = {
        "corrections": {},
        "surface_type": {},
        "retracker": {},
        "mss": {},
        "filter": {},
        "post_processing": {},
        "output": {}}

    # Assemble the job order
    job = Level2Job()
    job.mission_settings(mission_settings)
    job.roi_settings(roi_settings)
    job.l2proc_settings(l2_settings)
    job.validate()

    # Start the processor
    l2proc = Level2Processor(job)
    l2proc.set_config(config)
    l2proc.error_handling(raise_on_error=True)
    l2proc.set_l1b_files(l1b_files[0:1])
    l2proc.run()

if __name__ == '__main__':
    l2_processing()
