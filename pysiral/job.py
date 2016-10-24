# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 14:04:11 2015

@author: Stefan
"""

from pysiral.logging import DefaultLoggingClass
from pysiral.config import td_branches, ConfigInfo
from treedict import TreeDict
from datetime import datetime
import glob
import sys
import os


class ProcJob(DefaultLoggingClass):

    def __init__(self):
        super(ProcJob, self).__init__("procjob")
        self.pysiral_def = ConfigInfo()

    def mission_settings(self, config):
        if type(config) is dict:
            self._add_option_dict("mission", config)
        elif type(config) is TreeDict:
            self.mission = config

    def roi_settings(self, config):
        if type(config) is dict:
            self._add_option_dict("roi", config)
        elif type(config) is TreeDict:
            self.roi = config

    def local_machine_settings(self, config):
        if type(config) is dict:
            self._add_option_dict("local_machine", config)
        elif type(config) is TreeDict:
            self.local_machine = config

    def _add_option_dict(self, name, opt_dict):
        setattr(self, name, TreeDict.fromdict(opt_dict, expand_nested=True))


class Level2Job(ProcJob):

    def __init__(self):
        super(Level2Job, self).__init__()
        self.log.name = "Level2Job"
        self.overwrite_protection = True

    def l2proc_settings(self, config):
        if type(config) is dict:
            self._add_option_dict("config", config)
        elif type(config) is TreeDict:
            self.config = config

    def set_overwrite_protection(self, flag):
        self.overwrite_protection = flag

    def validate(self):
        self.log.info("validate and expand auxdata")
        self._validate_and_expand_auxdata()
        self._validate_and_create_output_directory()
        self.log.info("create output directory")

    def _validate_and_expand_auxdata(self):
        """
        - Verifies auxdata information in config.auxdata with
          content of config/auxdata_def.yaml
        - Transfer relevant content of auxdata_def.yaml into configuration
          structure
        """
        auxtypes, auxinfos = td_branches(self.config.auxdata)
        for auxtype, auxinfo in zip(auxtypes, auxinfos):
            self.log.info("Validating setting for %s: %s" % (
                auxtype, auxinfo.name))
            # Locate auxdata setting in config/auxdata_def.yaml
            try:
                pysiral_def = self.pysiral_def.auxdata[auxtype][auxinfo.name]
            except:
                msg = "id %s for type %s not in auxdata_def.yaml, EXITING" % (
                    auxinfo.name, auxtype)
                self.log.error(msg)
                sys.exit(1)
            # Repace local machine directory placeholder with actual directory
            auxdata_def = self.pysiral_def.local_machine.auxdata_repository
            auxdata_id = pysiral_def.local_repository
            if auxdata_id is not None:
                local_repository = auxdata_def[auxtype][auxdata_id]
                pysiral_def.local_repository = local_repository
            # Expand the settings with information on location of files, ...
            self.config.auxdata[auxtype].update(pysiral_def)

    def _validate_and_create_output_directory(self):
        output_ids, output_defs = td_branches(self.config.output)
        export_path = self.pysiral_def.local_machine.product_repository
        export_path = os.path.join(export_path, self.config.run_tag)
        time = datetime.now()
        tstamp = time.strftime("%Y%m%dT%H%M%S")
        for output_id, output_def in zip(output_ids, output_defs):
            if self.overwrite_protection:
                product_export_path = os.path.join(
                   export_path, tstamp, output_id)
            else:
                product_export_path = os.path.join(export_path, output_id)
            self.config.output[output_id].path = product_export_path
            self.log.info("Exporting %s data in directory %s" % (
                output_id, product_export_path))


class Level3Job(ProcJob):

    def __init__(self):
        super(Level3Job, self).__init__()
        self.log.name = "Level3Job"
        self.year = None
        self.month = None
        self.period = "month"
        self._input_directory = None
        self._grid_def = None
        self._l2_file_filter = "l2i*nc"
        self._l2_parameter = []
        self._l2_parameter = []
        self._frb_nanmask = []

    def set_input_directory(self, directory):
        self._input_directory = directory

    def set_grid_definition(self, setting):
        self._grid_def = setting

    def set_parameter(self, l2=[], l3=[], frb_nanmask=[]):
        self._l2_parameter = l2
        self._l3_parameter = l3
        self._frb_nanmask = frb_nanmask

    def get_monthly_l2idata_files(self, year, month, l2tag="l2i"):
        """
        Returns a sorted list of level-2 files
        (set_input_directory needs to be called first)
        """
        if self._input_directory is None:
            self.log.warning(
                "get_monthly_l2idata_files called without calling " +
                "\'set_input_directory\' first, return list is empty")
            return []
        self.year = year
        self.month = month
        # Get the file search pattern
        l2search = self._input_directory
        l2search = os.path.join(l2search, l2tag, "%04g" % year, "%02g" % month)
        l2search = os.path.join(l2search, self._l2_file_filter)
        # Query the folder
        l2files = sorted(glob.glob(l2search))
        return l2files

    def validate(self):
        pass

    @property
    def grid(self):
        return self._grid_def

    @property
    def l2_parameter(self):
        return self._l2_parameter

    @property
    def l3_parameter(self):
        return self._l3_parameter

    @property
    def hemisphere(self):
        return self._grid_def.hemisphere

    @property
    def resolution_tag(self):
        return self._grid_def.resolution_tag

    @property
    def grid_tag(self):
        return self._grid_def.grid_tag

    @property
    def freeboard_nan_mask_targets(self):
        return self._frb_nanmask

    @property
    def export_folder(self):
        export_folder = self._input_directory
        export_folder = os.path.join(export_folder, "l3s", "%04g" % self.year)
        return export_folder
