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
        
        # function to identify the first maximum index (fmi)
        def id_fmi(wfm: np.array, 
                   start: int = 0, 
                   bins: int = 128, 
                   spacing: int = 2) -> int:
            xshape = n_bins + 2 * spacing
            x = np.ndarray(shape=(xshape))
        
            x[:spacing] = wfm[0]-1.e-6
            x[-spacing:] = wfm[-1]-1.e-6
            x[spacing : spacing + n_bins] = wfm
        
            #before
            h_b = x[start : start + n_bins]
            start = spacing
            #central
            h_c = x[start : start + n_bins]
            start = spacing + 1
            #after
            h_a = x[start : start + n_bins]
        
            #find peak candidates
            peak_candidate = np.logical_and(h_c > h_b, h_c > h_a)
            #exclude first 10 bins
            peak_candidate = np.logical_and(peak_candidate, np.arange(n_bins) > 10)
            #at least half power is needed to be considered
            peak_candidate = np.logical_and(peak_candidate, wfm > 0.5)
        
            try: 
                fmi = np.where(peak_candidate)[0][0]
            except IndexError:
                fmi = np.nan
            return fmi

        # function to retrieve the sub waveform to be used
        def get_subset_wfm(wfm: np.array, 
                           fmi: int, 
                           i0: int = 5, 
                           i1: int = 30) -> np.array:
            #zero-pad waveforms
            wfm = np.concatenate((wfm, np.zeros(50)))
            try:
                sub = wfm[(fmi-i0):(fmi+i1)]
            except TypeError:
                sub = np.repeat(np.nan, i0+i1)
            return sub
        
        # get normalized waveform power
        waveform_norms = bn.nanmax(l1p.waveform.power, axis=1)
        normed_waveform_power = l1p.waveform.power[:]/waveform_norms[:, np.newaxis]
                
        #get waveform/sensor meta
        n_bins = normed_waveform_power.shape[1]
        n_wf = normed_waveform_power.shape[0]
        if n_bins>=128:
            i0 = 10
            i1 = 35
        else:
            i0 = 5
            i1 = 35
        
        # identify first-maximum index
        fmi = np.array([id_fmi(x) for x in normed_waveform_power])
        # subset waveforms
        sub_waveform_power = np.array([get_subset_wfm(x,i) 
                                       for x,i in zip(normed_waveform_power, fmi)])

        # create stack of five waveforms
        window_size = 5
        wfm_roll = []
        wfm_roll.extend([sub_waveform_power[i:(i + window_size)] 
                         for i in range(len(sub_waveform_power) - window_size + 1)])
        self.waveform_for_prediction = torch.stack(wfm_roll.astype('float32').astype('float32'))

        # as well as for the parameters


        #sub_waveform_power = np.zeros((n_wf,i0+i1))
        #c = 0
        #for x in normed_waveform_power:
        #    try:
        #        wf_fmi = id_fmi(x, bins=n_bins)
        #    except IndexError:
        #        c += 1
        #        continue
        #    if wf_fmi >= 10:
        #        sub_waveform_power[c] = get_subset_wfm(x,wf_fmi,i0,i1)
        #    c += 1
        #self.waveform_for_prediction = torch.tensor(sub_waveform_power.astype('float32')).unsqueeze(1)
        #self.waveform_for_prediction = torch.tensor(normed_waveform_power.astype('float32'))


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
        init.kaiming_uniform_(self.fc1.weight, nonlinearity='leaky_relu')
        init.kaiming_uniform_(self.fc2.weight, nonlinearity='leaky_relu')
        init.kaiming_uniform_(self.fc3.weight, nonlinearity='leaky_relu')
        init.kaiming_uniform_(self.fc4.weight, nonlinearity='leaky_relu')
        init.kaiming_uniform_(self.fc5.weight, nonlinearity='leaky_relu')
        init.kaiming_uniform_(self.fc6.weight, nonlinearity='leaky_relu')
        init.kaiming_uniform_(self.fc7.weight, nonlinearity='leaky_relu')
        init.kaiming_uniform_(self.fc8.weight, nonlinearity='leaky_relu')
        init.kaiming_uniform_(self.fc9.weight, nonlinearity='leaky_relu')
        init.kaiming_uniform_(self.fc1_par.weight, nonlinearity='leaky_relu')
        init.kaiming_uniform_(self.fc2_par.weight, nonlinearity='leaky_relu')
        
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



class ERS2_TestCandidate_001_FNN_TanH(nn.Module):
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
        init.xavier_uniform_(self.fc1.weight, gain=init.calculate_gain('tanh'))
        init.xavier_uniform_(self.fc2.weight, gain=init.calculate_gain('tanh'))
        init.xavier_uniform_(self.fc3.weight, gain=init.calculate_gain('tanh'))
        init.xavier_uniform_(self.fc4.weight, gain=init.calculate_gain('tanh'))
        init.xavier_uniform_(self.fc5.weight, gain=init.calculate_gain('tanh'))
        init.xavier_uniform_(self.fc6.weight, gain=init.calculate_gain('tanh'))
        init.xavier_uniform_(self.fc7.weight, gain=init.calculate_gain('tanh'))
        init.xavier_uniform_(self.fc8.weight, gain=init.calculate_gain('tanh'))
        init.xavier_uniform_(self.fc9.weight, gain=init.calculate_gain('tanh'))
        init.xavier_uniform_(self.fc1_par.weight, gain=init.calculate_gain('tanh'))
        init.xavier_uniform_(self.fc2_par.weight, gain=init.calculate_gain('tanh'))
        
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


