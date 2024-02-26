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
from typing import Any, Iterable, Union

import bottleneck as bn
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as torch_nn_functional
import xgboost as xgb

from pysiral import get_cls
from pysiral.auxdata import AuxdataBaseClass
from pysiral.l1data import L1bdataNCFile
from pysiral.l2data import Level2Data

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"


class RetrackerThresholdModel(AuxdataBaseClass):
    """ NOTE: This class is very likely obsolete due to switch to pytorch """

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
        model_filepath = [model_file for model_file in model_files if re.search(model_id, str(model_file))]

        # At this point there should only be one file
        if len(model_filepath) != 1:
            msg = f"No or multiple model files found for model_id = {model_id}: {self.cfg.local_repository}"
            self.error.add_error("ambigous-input", msg)
            self.error.raise_on_error()

        # Save input
        self.model_filepath = model_filepath[0]

        # NOTE: Danger zone - allowing multi-threading at which level?
        self.model_input = None
        self.model = xgb.Booster({'nthread': 1})
        self.model.load_model(str(self.model_filepath))

    def receive_l1p_input(self, l1p: 'L1bdataNCFile') -> None:
        """
        Optional method to add l1p variables to this class before `get_l2_track_vars()` is
        called. This method here grabs the necessary classifier parameter for the model
        prediction. The parameter names and their data group in the l1p data structure need to
        be specified in the Level-2 processor definition file as a list of [data_group, parameter_name],
        e.g.:

            parameter:
                - ["classifier", "leading_edge_width"]
                - ["classifier", "sigma0"]

        For this example the parameter list will be reorderd to a feature vector:

            [ [leading_edge_width[0], sigma0[0]], [leading_edge_width[1], sigma0[1]], ... ]

        :param l1p:
        :return:
        """

        # Retrieve required parameters for the model prediction
        parameter_list = self.cfg.options.get("parameter", [])
        if len(parameter_list) == 0:
            msg = f"Missing or empty option value parameter: {self.cfg.options}"
            self.error.add_error("missing-option", msg)
            self.error.raise_on_error()

        # Retrieve data from l1p data object
        model_input_list = []
        for data_group, parameter_name in parameter_list:
            variable = l1p.get_parameter_by_name(data_group, parameter_name)
            variable[np.isinf(variable)] = np.nan
            model_input_list.append(variable)

        # Reform the input data for the model prediction
        self.model_input = np.array(model_input_list).T

    def get_l2_track_vars(self, l2: 'Level2Data') -> None:
        """
        [Mandatory class method] Add the model prediction for the tfmra retracker threshold to the
        Level-2 data object. The model evaluation depends solely on waveform shape parameters that
        are not part of the Level-2 data object, because auxiliary data classes are called before the
        processor steps.
        :param l2:
        :return:
        """

        # Evaluate the model for this case
        model_data = xgb.DMatrix(self.model_input, missing=np.nan)
        prediction = self.model.predict(model_data)

        # Add prediction to the Level-2 data object
        var_id, var_name = self.cfg.options.get("output_parameter", ["tfmrathrs_ml", "tfmra_threshold_ml"])
        l2.set_auxiliary_parameter(var_id, var_name, prediction)

        # Remove the model input (safety measure)
        self.model_input = None

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

        if not model_files:
            msg = f"Did not find any machine learned files: {lookup_dir}/*{suffixes}"
            self.error.add_error("files-missing", msg)
            self.error.raise_on_error()

        return model_files


class RetrackerThresholdModelTorch(AuxdataBaseClass):

    def __init__(self, *args: Iterable[Any], **kwargs: Iterable[Any]) -> None:
        """
        Initialiaze the class. This step includes establishing the model by parsing the
        model parameter file as specified in the Level-2 processor definition file.
        :param args:
        :param kwargs:
        """
        super(RetrackerThresholdModelTorch, self).__init__(*args, **kwargs)

        # Retrieve requested model files
        self.model_file = self.cfg.options.get("model_file", None)
        if self.model_file is None:
            msg = f"Missing option `model_file` in auxiliary data configuration {self.cfg.options}"
            self.error.add_error("missing-option", msg)
            self.error.raise_on_error()

        self.model_filepath = Path(self.cfg.local_repository) / self.model_file
        if self.model_file is None:
            msg = f"Model file {self.model_filepath} not found"
            self.error.add_error("missing-file", msg)
            self.error.raise_on_error()

        # The model uses waveform power as input
        self.waveform_for_prediction = None

        # Get and initialize the required torch model class
        # REQ: Needs to be part of this file
        torch_model_class = self.cfg.options.get("torch_class", None)
        if torch_model_class is None:
            msg = "PyTorch model class not specified (options.torch_class missing)"
            self.error.add_error("missing-option", msg)
            self.error.raise_on_error()
        model_class, err = get_cls("pysiral.auxdata.ml", torch_model_class)
        #retrieve number of input layer neurons, defaults to 45 (the Envisat standard)
        input_neurons = self.cfg.options.get("input_neurons", 45)
        if model_class is None:
            msg = f"PyTorch model class not found: pysiral.auxdata.ml.{torch_model_class}"
            self.error.add_error("class-not-found", msg)
            self.error.raise_on_error()
        self.model = model_class(input_neurons)
        self.model.load_state_dict(torch.load(self.model_filepath,
                                              map_location=torch.device('cpu')))
        self.model.eval()

    def receive_l1p_input(self, l1p: 'L1bdataNCFile') -> None:
        """
        Optional method to add l1p variables to this class before `get_l2_track_vars()` is
        called. This method here store the waveform power arrays, which are input to the
        model, normalizes the waveforms and converts to the necessary data type.

        :param l1p:

        :return:
        """
        
        #function to identify the first maximum index (fmi)
        def id_fmi(wfm, start=0, bins=128, spacing=2):
            xshape = bins+2*spacing
            x = np.ndarray(shape=(xshape))

            x[:spacing] = wfm[0]-1.e-6
            x[-spacing:] = wfm[-1]-1.e-6
            x[spacing:spacing+bins] = wfm
            
            h_b = x[start:start + bins]  # before
            start = spacing
            h_c = x[start:start + bins]  # central
            start = spacing + 1
            h_a = x[start:start + bins]  # after
            
            peak_candidate = np.logical_and(h_c > h_b, h_c >= h_a)
            peak_candidate = np.logical_and(peak_candidate, np.arange(bins) > 10)
            peak_candidate = np.logical_and(peak_candidate, wfm > 0.5)

            return np.where(peak_candidate)[0][0]

        #function to retrieve the sub waveform to be used
        def get_sub_wf(wfm, fmi, i0=10, i1=35):
            #zero-pad waveforms
            wfm = np.concatenate((wfm,np.zeros(50)))
            return wfm[(fmi-i0):(fmi+i1)]
        
        #get normalized waveform power
        waveform_norms = bn.nanmax(l1p.waveform.power, axis=1)
        normed_waveform_power = l1p.waveform.power[:]/waveform_norms[:, np.newaxis]
                
        #get waveform/sensor meta
        n_bins = normed_waveform_power.shape[1]
        n_wf = normed_waveform_power.shape[0]
        if n_bins>=128:
            i0 = 10
            i1 = 35
        else:
            i0 = 10
            i1 = 25
        
        #subset waveform
        sub_waveform_power = np.zeros((n_wf,i0+i1))
        c = 0
        for x in normed_waveform_power:
            try:
                wf_fmi = id_fmi(x, bins=n_bins)
            except IndexError:
                c += 1
                continue
            if wf_fmi >= 10:
                sub_waveform_power[c] = get_sub_wf(x,wf_fmi,i0,i1)
            c += 1
        self.waveform_for_prediction = torch.tensor(sub_waveform_power.astype('float32')).unsqueeze(1)
#        self.waveform_for_prediction = torch.tensor(normed_waveform_power.astype('float32'))


    def get_l2_track_vars(self, l2: 'Level2Data') -> None:
        """
        [Mandatory class method] Add the model prediction for the tfmra retracker threshold to the
        Level-2 data object. The model evaluation depends solely on waveform power.

        :param l2: `Level2Data` container

        :return: None

        :raises: None
        """

        # Predict the waveform range
        with torch.no_grad():
            opt = self.model(self.waveform_for_prediction)
        tfmra_threshold_predicted = opt.numpy().flatten()

        # Limit threshold range to pre-defined range (or at least [0-1])
        valid_min, valid_max = self.cfg.options.get("valid_range", [0.0, 1.0])
        tfmra_threshold_predicted[tfmra_threshold_predicted < valid_min] = valid_min
        tfmra_threshold_predicted[tfmra_threshold_predicted > valid_max] = valid_max

        # Add prediction to the Level-2 data object
        var_id, var_name = self.cfg.options.get("output_parameter", ["tfmrathrs_ml", "tfmra_threshold_ml"])
        l2.set_auxiliary_parameter(var_id, var_name, tfmra_threshold_predicted)


# class TorchFunctionalWaveformModel(nn.Module):
#     """
#     Create Feed-Forward Neural Network architecture
#     REQ: Required for RetrackerThresholdModelTorch
#     """
#     def __init__(self):
#         super(TorchFunctionalWaveformModel, self).__init__()
#         self.fc1 = nn.Linear(128, 256)
#         self.bn1 = nn.BatchNorm1d(256)
#         self.fc2 = nn.Linear(256, 512)
#         self.bn2 = nn.BatchNorm1d(512)
#         self.fc3 = nn.Linear(512, 256)
#         self.bn3 = nn.BatchNorm1d(256)
#         self.fc4 = nn.Linear(256, 128)
#         self.bn4 = nn.BatchNorm1d(128)
#         self.fc5 = nn.Linear(128, 1)
#
#     def forward(self, x):
#         x = self.fc1(x)
#         x = self.bn1(x)
#         x = torch_nn_functional.leaky_relu(x)
#         x = self.fc2(x)
#         x = self.bn2(x)
#         x = torch_nn_functional.leaky_relu(x)
#         x = self.fc3(x)
#         x = self.bn3(x)
#         x = torch_nn_functional.leaky_relu(x)
#         x = self.fc4(x)
#         x = self.bn4(x)
#         x = torch_nn_functional.leaky_relu(x)
#         x = self.fc5(x)
#         return x

class TorchFunctionalWaveformModelFNN(nn.Module):
    """
    Create Feed-Forward Neural Network architecture
    REQ: Required for RetrackerThresholdModelTorch
    """
    def __init__(self, fc1_input=45):
        super(TorchFunctionalWaveformModelFNN, self).__init__()
        self.fc1 = nn.Linear(fc1_input, 2048)
        #self.bn1 = nn.BatchNorm1d(256)
        self.fc2 = nn.Linear(2048, 2048)
        #self.bn2 = nn.BatchNorm1d(512)
        self.fc3 = nn.Linear(2048, 2048)
        #self.bn3 = nn.BatchNorm1d(256)
        self.fc4 = nn.Linear(2048, 2048)
        #self.bn4 = nn.BatchNorm1d(128)
        self.fc5 = nn.Linear(2048, 2048)
        self.fc6 = nn.Linear(2048, 1)

    def forward(self, x):
        x = self.fc1(x)
        #x = self.bn1(x)
        x = torch_nn_functional.relu(x)
        x = self.fc2(x)
        #x = self.bn2(x)
        x = torch_nn_functional.relu(x)
        x = self.fc3(x)
        #x = self.bn3(x)
        x = torch_nn_functional.relu(x)
        x = self.fc4(x)
        #x = self.bn4(x)
        x = torch_nn_functional.relu(x)
        x = self.fc5(x)
        x = torch_nn_functional.relu(x)
        x = self.fc6(x)
        x = torch.sigmoid(x)
        return x


class TorchFunctionalWaveformModelSNN(nn.Module):
    """
    Create Self-Normalizing Neural Network architecture
    REQ: Required for RetrackerThresholdModelTorch
    """
    def __init__(self, fc1_input=45):
        super(TorchFunctionalWaveformModelSNN, self).__init__()
        self.fc1 = nn.Linear(fc1_input, 512)
        self.fc2 = nn.Linear(512, 1024)
        self.fc3 = nn.Linear(1024, 1024)
        self.fc4 = nn.Linear(1024, 1024)
        self.fc5 = nn.Linear(1024, 512)
        self.fc6 = nn.Linear(512, 1)
    def forward(self, x): 
        x = self.fc1(x)
        x = torch_nn_functional.selu(x)
        x = self.fc2(x)
        x = torch_nn_functional.selu(x)
        x = self.fc3(x)
        x = torch_nn_functional.selu(x)
        x = self.fc4(x)
        x = torch_nn_functional.selu(x)
        x = self.fc5(x)
        x = torch_nn_functional.selu(x)
        x = self.fc6(x)
        return x


class TorchFunctionalWaveformModelCNN(nn.Module):
    '''
    Create Convolutional Neural Network architecture
    REQ: Required for RetrackerThresholdModelTorch
    '''
    def __init__(self):
        super(TorchFunctionalWaveformModelCNN, self).__init__()
        self.cv1_1 = nn.Conv1d( 1, 32, kernel_size=(5,), stride=(1,), padding=True)
        self.cv1_2 = nn.Conv1d(32, 64, kernel_size=(5,), stride=(1,), padding=True)
        self.cv1_3 = nn.Conv1d(64, 128, kernel_size=(5,), stride=(1,), padding=True)
        self.mp1 = nn.MaxPool1d(3, stride=2)
        self.fc1 = nn.Linear(1792, 4096)
        self.fc2 = nn.Linear(4096, 1)
        #self.cv1_1 = nn.Conv1d(1, 64, kernel_size=5, stride=1, padding=True)
        #self.cv1_2 = nn.Conv1d(64, 64, kernel_size=5, stride=1, padding=True)
        #self.cv1_3 = nn.Conv1d(64, 64, kernel_size=5, stride=1, padding=True)
        #self.mp1 = nn.MaxPool1d(3, stride=2)
        #self.cv2_1 = nn.Conv1d(64, 128, kernel_size=5, stride=1, padding=True)
        #self.cv2_2 = nn.Conv1d(128, 128, kernel_size=5, stride=1, padding=True)
        #self.cv2_3 = nn.Conv1d(128, 128, kernel_size=5, stride=1, padding=True)
        #self.mp2 = nn.MaxPool1d(3, stride=2)
        #self.fc1 = nn.Linear(384, 2048)
        #self.fc2 = nn.Linear(2048, 1)
        
    def forward(self, x): 
        x = self.cv1_1(x)
        x = torch_nn_functional.relu(x)
        x = self.cv1_2(x)
        x = torch_nn_functional.relu(x)
        x = self.cv1_3(x)
        x = torch_nn_functional.relu(x)
        x = self.mp1(x)
        #x = self.cv2_1(x)
        #x = torch_nn_functional.relu(x)
        #x = self.cv2_2(x)
        #x = torch_nn_functional.relu(x)
        #x = self.cv2_3(x)
        #x = torch_nn_functional.relu(x)
        #x = self.mp2(x)
        x = x.view(x.size(0), -1)
        x = self.fc1(x)
        x = torch_nn_functional.relu(x)
        x = self.fc2(x)
        x = torch.sigmoid(x)
        return x


