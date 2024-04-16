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
import torch.nn.init as torch_nn_init
import torch.nn.functional as torch_nn_functional

from pysiral import get_cls
from pysiral.auxdata import AuxdataBaseClass
from pysiral.l1data import L1bdataNCFile
from pysiral.l2data import Level2Data

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"


class RetrackerThresholdModel(AuxdataBaseClass):
    """
    This class replaces the previous same-name version but utilizizes solely pytorch
    as pysiral's ML/AI framework for retracker-threshold models
    """
    
    def __init__(self, *args: Iterable[Any], **kwargs: Iterable[Any]) -> None:
        """
        Initialiaze the class. This step includes establishing the model by parsing the
        model parameter file as specified in the Level-2 processor definition file.
        :param args:
        :param kwargs:
        """
        super(RetrackerThresholdModel, self).__init__(*args, **kwargs)

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
        torch_model_class = self.cfg.options.get("model_class", None)
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
        self.model = model_class()
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

        # function to retrieve the sub waveform to be used
        def get_subset_wfm(wfm: np.array, 
                           fmi: int, 
                           i0: int = 5, 
                           i1: int = 30) -> np.array:
            #check for sufficient fmi position
            if fmi >= 30:
                #zero-pad waveforms
                wfm = np.concatenate((wfm, np.zeros(50)))
                try:
                    sub = wfm[(fmi-i0):(fmi+i1)]
                except TypeError:
                    sub = np.repeat(np.nan, i0+i1)
            else:
                sub = np.repeat(np.nan, i0+i1)
            return sub
        
        # get normalized waveform power
        power_max = bn.nanmax(l1p.waveform.power, axis=1)
        normed_waveform_power = l1p.waveform.power[:]/power_max[:, np.newaxis]

        # get l1p parameter
        eps = l1p.classifier.epsilon_sec
        eps = eps * 0.5 * 299792458.
        fmi = l1p.classifier.first_maximum_index
        
        # get settings for leading/trailing bins of fmi
        i0 = self.cfg.options.fmi_leading_bins
        i1 = self.cfg.options.fmi_trailing_bins

        # get nornmalization parameters based on TDS
        p_eps = self.cfg.options.classifiers.get("epsilon", [0, 1])
        p_pmax = self.cfg.options.classifiers.get("waveform_max_power", [0, 1])
        
        # subset waveforms
        sub_waveform_power = np.array([get_subset_wfm(x, i, i0, i1) 
                                       for x,i in zip(normed_waveform_power, fmi)])
        
        # only use valid waveforms w/ a FMI >= 30
        self.valid_waveforms_idx = np.where(fmi>=30)[0]
        # limit waveforms and parameters to valid ones
        valid_waveforms = sub_waveform_power[self.valid_waveforms_idx]
        valid_eps = eps[self.valid_waveforms_idx]
        valid_pmax = power_max[self.valid_waveforms_idx]
        # keep total number of waveforms
        self.n_total_waveforms = sub_waveform_power.shape[0]

        if len(self.valid_waveforms_idx)>=5:
            # create stack of five waveforms
            window_size = 5
            wfm_roll = []
            wfm_roll.extend([valid_waveforms[i:(i + window_size)] 
                             for i in range(len(valid_waveforms) - window_size + 1)])
            # keep for L2 processing
            self.waveform_for_prediction = torch.from_numpy(np.array(wfm_roll)).to(torch.float32)
    
            # make sure the indices are correct by putting them through the same procedure
            idx_roll = []
            idx_roll.extend([self.valid_waveforms_idx[i:(i + window_size)] 
                             for i in range(len(self.valid_waveforms_idx) - window_size + 1)])
            self.valid_waveforms_idx = np.array(idx_roll)[:,2]

            # stack as well as for the parameters
            eps_roll = []
            pmax_roll = []
            eps_roll.extend([valid_eps[i:(i + window_size)] 
                             for i in range(len(valid_eps) - window_size + 1)])
            pmax_roll.extend([valid_pmax[i:(i + window_size)] 
                             for i in range(len(valid_pmax) - window_size + 1)])
            eps = torch.from_numpy(np.array(eps_roll)).to(torch.float32)
            pmax = torch.from_numpy(np.array(pmax_roll).astype(np.float32)/1000)
            # normalize parameter
            eps = (eps - p_eps[0]) / (p_eps[1] - p_eps[0])
            pmax = (pmax - p_pmax[0]) / (p_pmax[1] - p_pmax[0])
            # keep for L2 processing       
            self.parameters_for_prediction = torch.cat([eps, pmax], dim=1)
        else:
            self.waveform_for_prediction = None
            self.parameters_for_prediction = None


    def get_l2_track_vars(self, l2: 'Level2Data') -> None:
        """
        [Mandatory class method] Add the model prediction for the tfmra retracker threshold to the
        Level-2 data object. The model evaluation depends solely on waveform power.

        :param l2: `Level2Data` container

        :return: None

        :raises: None
        """

        # Predict the waveform range
        if self.waveform_for_prediction is not None:  
            with torch.no_grad():
                if hasattr(self, 'parameters_for_prediction'):
                    opt = self.model(self.waveform_for_prediction, self.parameters_for_prediction)
                else:
                    opt = self.model(self.waveform_for_prediction)
            tfmra_threshold_predicted = opt.numpy().flatten()
    
            # Limit threshold range to pre-defined range (or at least [0-1])
            valid_min, valid_max = self.cfg.options.get("valid_range", [0.0, 1.0])
            tfmra_threshold_predicted[tfmra_threshold_predicted < valid_min] = valid_min
            tfmra_threshold_predicted[tfmra_threshold_predicted > valid_max] = valid_max
    
            # Return predicted thresholds back on original trajectory
            tfmra_on_trajectory = np.full(self.n_total_waveforms, np.nan)
            tfmra_on_trajectory[self.valid_waveforms_idx] = tfmra_threshold_predicted
        else:
            # Create placeholder
            tfmra_on_trajectory = np.full(self.n_total_waveforms, np.nan)
        # Add prediction to the Level-2 data object
        var_id, var_name = self.cfg.options.get("output_parameter", ["tfmrathrs_ml", "tfmra_threshold_ml"])
        l2.set_auxiliary_parameter(var_id, var_name, tfmra_on_trajectory)


class ERS2_TestCandidate_001_FNN_LeakyRelu(nn.Module):
    '''
    Creates Feed-Forward Neural Network architecture using leaky_relu activation function
    using two separate branches for 1) the input of the subset waveform power of a five 
    waveform stack and 2) the parameter input of Epsilon and Max Power for the respective
    waveforms of the stack. 
    
    Both branches are then input and feed through a common two-layer network and output 
    decimal optimal retracker thresholds through a sigmoid activation finish. 
    
    Weights of all layers are initialized using LeCun uniform (kaiming in pytorch) with 
    corresponding activaion function non-linearity.

    Desgined/Trained for ERS2 subwaveform input using 5 bins leading and 30 bins trailing 
    the identified first-maximum index (fmi).

    REQ: Required for RetrackerThresholdModel
    '''
    def __init__(self, n_in: int = 5, n_out: int = 1, n_par: int = 2):
        super(ERS2_TestCandidate_001_FNN_LeakyRelu, self).__init__()
        # number of input/output channels and parameters
        self.n_in = n_in
        self.n_out = n_out
        self.n_par = n_par
        # waveform branch
        self.fc1 = nn.Linear(self.n_in*35, 512)
        self.fc2 = nn.Linear(512, 512)
        self.fc3 = nn.Linear(512, 512)
        self.fc4 = nn.Linear(512, 512)
        self.fc5 = nn.Linear(512, 512)
        self.fc6 = nn.Linear(512, 512)
        # eps branch
        self.fc1_par = nn.Linear(self.n_in*self.n_par, 32)
        self.fc2_par = nn.Linear(32, 32)
        self.fc2_par = nn.Linear(32, 32)
        # combo
        self.fc7 = nn.Linear(544, 1024)
        self.fc8 = nn.Linear(1024, 1024)
        self.fc9 = nn.Linear(1024,  self.n_out)
        # initialize weights using LeCun initialization with relu nonlinearity
        torch.manual_seed(27570)
        torch_nn_init.kaiming_uniform_(self.fc1.weight, nonlinearity='leaky_relu')
        torch_nn_init.kaiming_uniform_(self.fc2.weight, nonlinearity='leaky_relu')
        torch_nn_init.kaiming_uniform_(self.fc3.weight, nonlinearity='leaky_relu')
        torch_nn_init.kaiming_uniform_(self.fc4.weight, nonlinearity='leaky_relu')
        torch_nn_init.kaiming_uniform_(self.fc5.weight, nonlinearity='leaky_relu')
        torch_nn_init.kaiming_uniform_(self.fc6.weight, nonlinearity='leaky_relu')
        torch_nn_init.kaiming_uniform_(self.fc7.weight, nonlinearity='leaky_relu')
        torch_nn_init.kaiming_uniform_(self.fc8.weight, nonlinearity='leaky_relu')
        torch_nn_init.kaiming_uniform_(self.fc9.weight, nonlinearity='leaky_relu')
        torch_nn_init.kaiming_uniform_(self.fc1_par.weight, nonlinearity='leaky_relu')
        torch_nn_init.kaiming_uniform_(self.fc2_par.weight, nonlinearity='leaky_relu')
        
    def forward(self, x, par):
        # waveform part
        x = x.view(-1, self.n_in*35)
        x = self.fc1(x)
        x = torch_nn_functional.leaky_relu(x)
        x = self.fc2(x)
        x = torch_nn_functional.leaky_relu(x)
        x = self.fc3(x)
        x = torch_nn_functional.leaky_relu(x)
        x = self.fc4(x)
        x = torch_nn_functional.leaky_relu(x)
        x = self.fc5(x)
        x = torch_nn_functional.leaky_relu(x)
        x = self.fc6(x)
        x = torch_nn_functional.leaky_relu(x)
        # inflow of parameter
        par_x = par.view(-1, self.n_in*self.n_par)
        par_x = self.fc1_par(par_x)        
        par_x = torch_nn_functional.leaky_relu(par_x)
        par_x = self.fc2_par(par_x)
        par_x = torch_nn_functional.leaky_relu(par_x)
        # combine both
        x = torch.cat([x, par_x], dim=1)
        x = self.fc7(x)
        x = torch_nn_functional.leaky_relu(x)
        x = self.fc8(x)
        x = torch_nn_functional.leaky_relu(x)
        x = self.fc9(x)
        x = torch.sigmoid(x)
        return x



class ERS2_TestCandidate_002_FNN_TanH(nn.Module):
    '''
    Creates Feed-Forward Neural Network architecture using tanh activation function
    using two separate branches for 1) the input of the subset waveform power of a five 
    waveform stack and 2) the parameter input of Epsilon and Max Power for the respective
    waveforms of the stack. 
    
    Both branches are then input and feed through a common two-layer network and output 
    decimal optimal retracker thresholds through a sigmoid activation finish. 
    
    Weights of all layers are initialized using Xavier uniform with corresponding activaion
    function gain.

    Desgined/Trained for ERS2 subwaveform input using 5 bins leading and 30 bins trailing 
    the identified first-maximum index (fmi).

    REQ: Required for RetrackerThresholdModel
    '''
    def __init__(self, n_in: int = 5, n_out: int = 1, n_par: int = 2):
        super(ERS2_TestCandidate_001_FNN_TanH, self).__init__()
        # number of input channels
        self.n_in = n_in
        self.n_out = n_out
        self.n_par = n_par
        # waveform branch
        self.fc1 = nn.Linear(self.n_in*35, 512)
        self.fc2 = nn.Linear(512, 512)
        self.fc3 = nn.Linear(512, 512)
        self.fc4 = nn.Linear(512, 512)
        self.fc5 = nn.Linear(512, 512)
        self.fc6 = nn.Linear(512, 512)
        # eps branch
        self.fc1_par = nn.Linear(self.n_in*self.n_par, 32)
        self.fc2_par = nn.Linear(32, 32)
        # combo
        self.fc7 = nn.Linear(544, 1024)
        self.fc8 = nn.Linear(1024, 1024)
        self.fc9 = nn.Linear(1024, self.n_out)
        # initialize weights using LeCun initialization with relu nonlinearity
        torch.manual_seed(27570)
        torch_nn_init.xavier_uniform_(self.fc1.weight, gain=init.calculate_gain('tanh'))
        torch_nn_init.xavier_uniform_(self.fc2.weight, gain=init.calculate_gain('tanh'))
        torch_nn_init.xavier_uniform_(self.fc3.weight, gain=init.calculate_gain('tanh'))
        torch_nn_init.xavier_uniform_(self.fc4.weight, gain=init.calculate_gain('tanh'))
        torch_nn_init.xavier_uniform_(self.fc5.weight, gain=init.calculate_gain('tanh'))
        torch_nn_init.xavier_uniform_(self.fc6.weight, gain=init.calculate_gain('tanh'))
        torch_nn_init.xavier_uniform_(self.fc7.weight, gain=init.calculate_gain('tanh'))
        torch_nn_init.xavier_uniform_(self.fc8.weight, gain=init.calculate_gain('tanh'))
        torch_nn_init.xavier_uniform_(self.fc9.weight, gain=init.calculate_gain('tanh'))
        torch_nn_init.xavier_uniform_(self.fc1_par.weight, gain=init.calculate_gain('tanh'))
        torch_nn_init.xavier_uniform_(self.fc2_par.weight, gain=init.calculate_gain('tanh'))
        
    def forward(self, x, par):
        # waveform part
        x = x.view(-1, self.n_in*35)
        x = self.fc1(x)
        x = torch_nn_functional.tanh(x)
        x = self.fc2(x)
        x = torch_nn_functional.tanh(x)
        x = self.fc3(x)
        x = torch_nn_functional.tanh(x)
        x = self.fc4(x)
        x = torch_nn_functional.tanh(x)
        x = self.fc5(x)
        x = torch_nn_functional.tanh(x)
        x = self.fc6(x)
        x = torch_nn_functional.tanh(x)
        # inflow of parameter
        par_x = par.view(-1, self.n_in*self.n_par)
        par_x = self.fc1_par(par_x)        
        par_x = torch_nn_functional.tanh(par_x)
        par_x = self.fc2_par(par_x)
        par_x = torch_nn_functional.tanh(par_x)
        # combine both
        x = torch.cat([x, par_x], dim=1)
        x = self.fc7(x)
        x = torch_nn_functional.tanh(x)
        x = self.fc8(x)
        x = torch_nn_functional.tanh(x)
        x = self.fc9(x)
        x = torch.sigmoid(x)
        return x
        

class TorchFunctionalWaveformModelFNN(nn.Module):
    """
    Create Feed-Forward Neural Network architecture
    REQ: Required for RetrackerThresholdModel
    
    Note: Legacy version of the current ENVISAT implementation until rework
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
    REQ: Required for RetrackerThresholdModel

    Note: Deprecated!
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
    REQ: Required for RetrackerThresholdModel

    Note: Deprecated!
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


