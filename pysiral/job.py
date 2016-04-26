# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 14:04:11 2015

@author: Stefan
"""

from pysiral.logging import DefaultLoggingClass
from pysiral.config import td_branches, ConfigInfo
from treedict import TreeDict
from datetime import datetime
import sys
import os


class ProcJob(DefaultLoggingClass):

    def __init__(self):
        super(ProcJob, self).__init__("procjob")
        self.pysiral_def = ConfigInfo()

    def mission_settings(self, config):
        self._add_option_dict("mission", config)

    def roi_settings(self, config):
        self._add_option_dict("roi", config)

    def local_machine_settings(self, config):
        self._add_option_dict("local_machine", config)

    def _add_option_dict(self, name, opt_dict):
        setattr(self, name, TreeDict.fromdict(opt_dict, expand_nested=True))


class Level2Job(ProcJob):

    def __init__(self):
        super(Level2Job, self).__init__()
        self.log.name = "Level2Job"

    def l2proc_settings(self, config):
        self._add_option_dict("config", config)

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
            product_export_path = os.path.join(export_path, tstamp, output_id)
            self.config.output[output_id].path = product_export_path
            self.log.info("Exporting %s data in directory %s" % (
                output_id, product_export_path))
