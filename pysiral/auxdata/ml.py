# -*- coding: utf-8 -*-
"""

module for ingesting models from machine learning

Important Note:

    All mdt data handlers must be subclasses of pysiral.auxdata.AuxdataBaseClass in order to work
    for the Level-2 Processor. If the auxiliary class is based on a static dataset, this should be parsed
    in `__init__`.

    Please review the variables and properties in the parent class, as well as the correspodning config and
    support classes for grid track interpolation in the pysiral.auxdata module for additional guidance.

    The only other hard requirements is the presence of on specific method in order to be a valid subclass of
    AuxdataBaseClass:


        get_l2_track_vars(l2)

            This method will be called during the Level-2 processor. The argument is the Level-2 data object and
            the purpose of the method is to compute the auxilary variable(s) and associated uncertainty. These
            variable need to be registered using the `register_auxvar(id, name, value, uncertainty)` method of
            the base class. All MDT subclasses need to register at minimum the following variable:

            mean dynamic topography (relative to MSS):
                id: mdt
                name: mean_dynamic_topography

            e.g., this code line is mandatory for `get_l2_track_vars` (uncertainty can be None):

                # Register Variables
                self.register_auxvar("mdt", "mean_dynamic_topography", value, uncertainty)

"""
import re
from pathlib import Path
from typing import Iterable, Any, Union

import xgboost as xgb

from pysiral.auxdata import AuxdataBaseClass
from pysiral.l2data import Level2Data
from pysiral.l1bdata import L1bdataNCFile

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"


class RetrackerThresholdModel(AuxdataBaseClass):

    def __init__(self, *args: Iterable[Any], **kwargs: Iterable[Any]) -> None:
        """
        Initialiaze the class. This step includes establishing the model by parsing the
        model parameter file as specified in the Level-2 processor definition file.
        :param args:
        :param kwargs:
        """
        super(RetrackerThresholdModel, self).__init__(*args, **kwargs)

        # Query available model file
        suffixes = self.cfg.options.get("suffixes", [])
        model_files = self.get_available_model_files(self.cfg.local_repository, suffixes)

        # Retrieve requested model files
        model_id = self.cfg.options.get("model_id", None)
        if model_id is None:
            msg = f"Missing option `model_id` in auxiliary data configuration {self.cfg.options}"
            self.error.add_error("missing-option", msg)
            self.error.raise_on_error()
        model_filepath = [model_file for model_file in model_files if re.search(model_id, str(model_files))]

        # Verify input
        if len(model_filepath) != 1:
            msg = f"No or multiple model files found for model_id = {model_id}: {model_filepath}"
            self.error.add_error("ambigous-input", msg)
            self.error.raise_on_error()

        # Save input
        self.model_filepath = model_filepath[0]

        # NOTE: Danger zone - allowing multi-threading at which level?
        self.model = xgb.Booster({'nthread': 1})
        self.model.load_model(str(self.model_filepath))

    def receive_l1p_input(self, l1p: 'L1bdataNCFile') -> None:
        """
        Optional method to add l1p variables to this class before `get_l2_track_vars()` is
        called
        :param l1p:
        :return:
        """

        breakpoint()

    def get_l2_track_vars(self, l2: 'Level2Data') -> None:
        """
        TODOC:
        :param l2:
        :return:
        """
        breakpoint()

    def get_available_model_files(self, lookup_dir: Union[str, Path], suffixes: Iterable[str]) -> Iterable[Path]:
        """
        Check the avaiable model files
        :param lookup_dir:
        :param suffixes:
        :return: List of available model files
        """

        # Test if directory exists
        lookup_dir = Path(lookup_dir)
        if not lookup_dir.is_dir():
            msg = f"Directory for machine learned models does not exist: {lookup_dir}"
            self.error.add_error("directory-missing", msg)
            self.error.raise_on_error()

        # Find files
        model_files = []
        for suffix in suffixes:
            model_file_list = list(lookup_dir.rglob(f"*{suffix}"))
            model_files.extend(model_file_list)

        if len(model_files) == 0:
            msg = f"Did not find any machine learned files: {lookup_dir}/*{suffixes}"
            self.error.add_error("files-missing", msg)
            self.error.raise_on_error()

        return model_files
