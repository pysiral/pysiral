# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 14:04:11 2015

@author: Stefan
"""

from treedict import TreeDict


class ProcJob(object):

    def __init__(self):
        pass

    def mission_settings(self, config):
        self._add_option_dict("mission", config)

    def roi_settings(self, config):
        self._add_option_dict("roi", config)

    def _add_option_dict(self, name, opt_dict):
        setattr(self, name, TreeDict.fromdict(opt_dict, expand_nested=True))


class Level2Job(ProcJob):

    def __init__(self):
        super(Level2Job, self).__init__()

    def l2proc_settings(self, config):
        self._add_option_dict("config", config)

    def validate():
        pass
