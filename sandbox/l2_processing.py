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
        "corrections": [
            "dry_troposphere",
            "wet_troposphere",
            "inverse_barometric",
            "dynamic_atmosphere",
            "ionospheric",
            "ocean_tide_elastic",
            "ocean_tide_long_period",
            "ocean_loading_tide",
            "solid_earth_tide",
            "geocentric_polar_tide"],
        "surface_type": {
            "pyclass": "RickerTC2014",
            "options": {
                "ocean": {
                    "peakiness_min": 0.0,
                    "peakiness_max": 10.0,
                    "stack_standard_deviation_min": 18.5,
                    "ice_concentration_min": 5.0,
                    "ocog_width_min": 38},
                "lead": {
                    "peakiness_l_min": 40.0,
                    "peakiness_r_min": 30.0,
                    "peakiness_min": 40.0,
                    "stack_kurtosis_min": 40.0,
                    "stack_standard_deviation_max": 4.0,
                    "ice_concentration_min": 70.0},
                "sea_ice": {
                    "peakiness_r_max": 15.0,
                    "peakiness_l_max": 20.0,
                    "peakiness_max": 30.0,
                    "stack_kurtosis_max": 8.0,
                    "ice_concentration_min": 70.0}}},
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
