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
from pathlib import Path
from typing import Any, Iterable

import bottleneck as bn
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as torch_nn_functional
from loguru import logger

from pysiral import get_cls
from pysiral.auxdata import AuxdataBaseClass

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
        self.valid_waveforms_idx = None
        self.n_total_waveforms = None
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
        self.parameters_for_prediction = None

        # Get and initialize the required torch model class
        # REQ: Needs to be part of this file
        torch_model_class = self.cfg.options.get("model_class", None)
        if torch_model_class is None:
            msg = "PyTorch model class not specified (options.torch_class missing)"
            self.error.add_error("missing-option", msg)
            self.error.raise_on_error()
        model_class, err = get_cls("pysiral.auxdata.ml", torch_model_class)
        # retrieve number of input layer neurons, defaults to 45 (the Envisat standard)
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
        def get_subset_wfm(wfm: np.ndarray,
                           fmi: int,
                           i0: int = 5,
                           i1: int = 30) -> np.ndarray:
            # check for sufficient fmi position
            if fmi >= 30:
                # zero-pad waveforms
                wfm = np.concatenate((wfm, np.zeros(50)))
                try:
                    sub = wfm[(fmi - i0):(fmi + i1)]
                except TypeError:
                    sub = np.repeat(np.nan, i0 + i1)
            else:
                sub = np.repeat(np.nan, i0 + i1)
            return sub

        # get normalized waveform power
        power_max = bn.nanmax(l1p.waveform.power, axis=1)
        normed_waveform_power = l1p.waveform.power[:]/power_max[:, np.newaxis]

        # get normalization parameters based on TDS from settings file
        if 'classifiers' in self.cfg.options.keys():
            normed_parameters = np.empty((normed_waveform_power.shape[0], 0))
            classifiers = self.cfg.options.classifiers.keys()
            for key in classifiers:
                logger.debug(f'- Imported Parameter: {key}')
                # get classifier parameter
                if key == 'latitude' or key == 'longitude':
                    par = getattr(l1p.time_orbit, key)
                else:
                    par = getattr(l1p.classifier, key)
                if key == 'epsilon_sec':
                    par = par * 0.5 * 299792458.
                # get normalization parameters from TDS (mean/std)
                norm_pars = self.cfg.options.classifiers.get(key, [0, 1])
                # normalize the classifier parameter
                # par_normed = (par - norm_pars[0]) / norm_pars[1]
                par_normed = (par - norm_pars[0]) / (norm_pars[1] - norm_pars[0])
                normed_parameters = np.column_stack((normed_parameters, par_normed))

        # get reference l1p parameter
        fmi = l1p.classifier.first_maximum_index.astype(int)
        ami = bn.nanargmax(l1p.waveform.power, axis=1).astype(int)

        # get settings for leading/trailing bins of fmi
        i0 = self.cfg.options.fmi_leading_bins
        i1 = self.cfg.options.fmi_trailing_bins

        # subset waveforms
        sub_waveform_power = np.array([get_subset_wfm(x, i, i0, i1)
                                       for x, i in zip(normed_waveform_power, fmi)])
        # only use valid waveforms w/ a FMI >= 30
        self.valid_waveforms_idx = np.where((fmi >= 30) &
                                            # (l1p.classifier.peakiness>3.5) &
                                            # (l1p.classifier.late_tail_to_peak_power < 0.35) &
                                            #  (l1p.classifier.trailing_edge_slope < -0.015) &
                                            # (l1p.classifier.leading_edge_width < 1.10) &
                                            (ami == fmi)
                                           )[0]
        logger.debug(f'- Number of valid waveforms: {len(self.valid_waveforms_idx)}')

        # check for invalid parameters (mainly on non-sea-ice waveforms)
        if 'classifiers' in self.cfg.options.keys():
            invalid_parameters_idx = np.unique(np.where(np.isnan(normed_parameters)))
            self.valid_waveforms_idx = np.setdiff1d(self.valid_waveforms_idx,
                                                    invalid_parameters_idx)

        # validate there a sea-ice waveforms
        if len(self.valid_waveforms_idx) > 0:
            # limit waveforms to valid ones
            valid_waveforms = sub_waveform_power[self.valid_waveforms_idx]
            # keep for L2 processing
            self.waveform_for_prediction = torch.from_numpy(valid_waveforms).to(torch.float32)
            # logger.debug(f'- Shape of Waveform Tensor: {self.waveform_for_prediction.size()}')

            if 'classifiers' in self.cfg.options.keys():
                # limit also parameters to valid ones
                valid_parameters = normed_parameters[self.valid_waveforms_idx, :]
                # keep for L2 processing
                self.parameters_for_prediction = torch.from_numpy(valid_parameters).to(torch.float32)
                # logger.debug(f'- Shape of Parameter Tensor: {self.parameters_for_prediction.size()}')

        else:
            self.waveform_for_prediction = None
            self.parameters_for_prediction = None
        # keep total number of waveforms
        self.n_total_waveforms = sub_waveform_power.shape[0]

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
            tfmra_threshold_predicted = opt.flatten().numpy()
            logger.debug(f'Predicted Thresholds Percentiles 10%: {np.around(np.quantile(tfmra_threshold_predicted,.1),2)}'+\
                         f'/ 90%: {np.around(np.quantile(tfmra_threshold_predicted,.9),2)}')

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

# --- Model Setups ---


class AutoEncoderERS2(nn.Module):
    def __init__(self, n_in: int = 35, n_bn: int = 3):
        super(AutoEncoderERS2, self).__init__()
        # number of input channels and bottleneck layer size
        self.n_in = n_in
        self.n_bn = n_bn
        # encoder
        self.encoder = nn.Sequential(
            nn.Linear(self.n_in, self.n_in*3),
            nn.LeakyReLU(),
            nn.Linear(self.n_in*3, self.n_in*2),
            nn.LeakyReLU(),
            nn.Linear(self.n_in*2, self.n_in*1),
            nn.LeakyReLU(),
            nn.Linear(self.n_in*1, self.n_in//2),
            nn.LeakyReLU(),
            nn.Linear(self.n_in//2, self.n_bn),
            nn.LeakyReLU()
        )

        # decoder
        self.decoder = nn.Sequential(
            nn.Linear(self.n_bn, self.n_in//2),
            nn.LeakyReLU(),
            nn.Linear(self.n_in//2, self.n_in*1),
            nn.LeakyReLU(),
            nn.Linear(self.n_in*1, self.n_in*2),
            nn.LeakyReLU(),
            nn.Linear(self.n_in*2, self.n_in*3),
            nn.LeakyReLU(),
            nn.Linear(self.n_in*3, self.n_in),
            nn.Sigmoid()
        )

    def forward(self, x):
        x = self.encoder(x)
        x = self.decoder(x)
        return x


class ERS2_TestCandidate_006_FNN(nn.Module):
    def __init__(self, n_in: int = 35, n_par=8):
        super(ERS2_TestCandidate_006_FNN, self).__init__()
        # number of input channels
        self.n_in = n_in
        self.n_par = n_par
        # model
        self.model = nn.Sequential(
            nn.Linear(self.n_in+self.n_par, 2048),
            nn.LeakyReLU(),
            nn.Linear(2048, 2048),
            nn.LeakyReLU(),
            nn.Linear(2048, 2048),
            nn.LeakyReLU(),
            nn.Linear(2048, 2048),
            nn.LeakyReLU(),
            nn.Linear(2048, 1),
            nn.Sigmoid(),
        )

    def forward(self, x, par):
        x = self.model(torch.cat([x, par], dim=1))
        return x


class ERS2_TestCandidate_005_FNN(nn.Module):
    def __init__(self, n_in: int = 35, n_par=6):
        super(ERS2_TestCandidate_005_FNN, self).__init__()
        # number of input channels
        self.n_in = n_in
        self.n_par = n_par
        # model
        self.model = nn.Sequential(
            nn.Linear(self.n_in+self.n_par, 2048),
            nn.LeakyReLU(),
            nn.Linear(2048, 2048),
            nn.LeakyReLU(),
            nn.Linear(2048, 2048),
            nn.LeakyReLU(),
            nn.Linear(2048, 1),
            nn.Sigmoid(),
        )

    def forward(self, x, par):
        x = self.model(torch.cat([x, par], dim=1))
        return x


class ERS2_TestCandidate_004_FNN(nn.Module):

    def __init__(self, n_in: int = 35, n_par: int = 6):
        super(ERS2_TestCandidate_004_FNN, self).__init__()
        # number of input channels
        self.n_in = n_in
        self.n_par = n_par
        # waveform model branch
        self.wfm = nn.Sequential(
            nn.Linear(self.n_in, self.n_in*5),
            nn.Tanh(),
            nn.Linear(self.n_in*5, self.n_in*5),
            nn.Tanh(),
            nn.Linear(self.n_in*5, self.n_in*5),
            nn.Tanh(),
            nn.Linear(self.n_in*5, self.n_in*5),
            nn.Tanh(),
        )
        # parameter model branch
        self.par = nn.Sequential(
            nn.Linear(self.n_par, self.n_par*5),
            nn.Tanh(),
            nn.Linear(self.n_par*5, self.n_par*5),
            nn.Tanh(),
            nn.Linear(self.n_par*5, self.n_par*5),
            nn.Tanh(),
            nn.Linear(self.n_par*5, self.n_par*5),
            nn.Tanh(),
        )
        # fusion model branch
        self.cmb = nn.Sequential(
            nn.Linear((self.n_in+self.n_par)*5, (self.n_in+self.n_par)*5),
            nn.Tanh(),
            nn.Linear((self.n_in+self.n_par)*5, (self.n_in+self.n_par)*5),
            nn.Tanh(),
            nn.Linear((self.n_in+self.n_par)*5, 1),
            nn.Sigmoid(),
        )

    def forward(self, wfm, par):
        wfm_x = self.wfm(wfm)
        par_x = self.par(par)  # .view(-1, self.n_par))
        x = self.cmb(torch.cat([wfm_x, par_x], dim=1))
        return x


class ERS2_TestCandidate_003_FNN_TanHLeakyRelu(nn.Module):
    """
    Create Feed-Forward Neural Network architecture
    """
    def __init__(self, n_in: int = 5, n_out: int = 1, n_par: int = 3):
        super(ERS2_TestCandidate_003_FNN_TanHLeakyRelu, self).__init__()
        # base setup
        self.n_in = n_in
        self.n_out = n_out
        self.n_par = n_par
        # layer setup
        in_lay_wfm = self.n_in*35
        lay_wfm = in_lay_wfm * 3
        in_lay_par = self.n_in*self.n_par
        lay_par = in_lay_par * 5
        in_lay_cmb = lay_wfm + lay_par
        lay_cmb = in_lay_cmb * 2
        # waveform branch
        self.fc1 = nn.Linear(in_lay_wfm, lay_wfm)
        self.fc2 = nn.Linear(lay_wfm, lay_wfm)
        self.fc3 = nn.Linear(lay_wfm, lay_wfm)
        self.fc4 = nn.Linear(lay_wfm, lay_wfm)
        self.fc5 = nn.Linear(lay_wfm, lay_wfm)
        # eps branch
        self.fc1_par = nn.Linear(in_lay_par, lay_par)
        self.fc2_par = nn.Linear(lay_par, lay_par)
        self.fc3_par = nn.Linear(lay_par, lay_par)
        self.fc4_par = nn.Linear(lay_par, lay_par)
        self.fc5_par = nn.Linear(lay_par, lay_par)
        # combo
        self.fc7 = nn.Linear(in_lay_cmb, lay_cmb)
        self.fc8 = nn.Linear(lay_cmb, lay_cmb)
        self.fc9 = nn.Linear(lay_cmb, self.n_out)

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
        # inflow of parameter
        par_x = par.view(-1, self.n_in*self.n_par)
        par_x = self.fc1_par(par_x)
        par_x = torch_nn_functional.tanh(par_x)
        par_x = self.fc2_par(par_x)
        par_x = torch_nn_functional.tanh(par_x)
        par_x = self.fc3_par(par_x)
        par_x = torch_nn_functional.tanh(par_x)
        par_x = self.fc4_par(par_x)
        par_x = torch_nn_functional.tanh(par_x)
        par_x = self.fc5_par(par_x)
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


class ERS2_TestCandidate_001_FNN_LeakyRelu(nn.Module):
    """
    Creates Feed-Forward Neural Network architecture using leaky_relu activation function
    using two separate branches for 1) the input of the subset waveform power of a five
    waveform stack and 2) the parameter input of Epsilon and Max Power for the respective
    waveforms of the stack.

    Both branches are then input and feed through a common two-layer network and output
    decimal optimal retracker thresholds through a sigmoid activation finish.

    Desgined/Trained for ERS2 subwaveform input using 5 bins leading and 30 bins trailing
    the identified first-maximum index (fmi).

    REQ: Required for RetrackerThresholdModel
    """
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
    """
    Creates Feed-Forward Neural Network architecture using tanh activation function
    using two separate branches for 1) the input of the subset waveform power of a five
    waveform stack and 2) the parameter input of Epsilon and Max Power for the respective
    waveforms of the stack.

    Both branches are then input and feed through a common two-layer network and output
    decimal optimal retracker thresholds through a sigmoid activation finish.

    Desgined/Trained for ERS2 subwaveform input using 5 bins leading and 30 bins trailing
    the identified first-maximum index (fmi).

    REQ: Required for RetrackerThresholdModel
    """
    def __init__(self, n_in: int = 5, n_out: int = 1, n_par: int = 2):
        super(ERS2_TestCandidate_002_FNN_TanH, self).__init__()
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


class ERS2_TestCandidate_003_LSTM_LeakyRelu(nn.Module):

    def __init__(self, n_in: int = 5, n_out: int = 1, n_par: int = 2):
        super(ERS2_TestCandidate_003_LSTM_LeakyRelu, self).__init__()
        # number of input channels
        self.n_in = n_in
        self.n_out = n_out
        self.n_par = n_par
        # waveform branch with lstm/fc [input_size, hidden_size, num_layers]
        self.lstm = nn.LSTM(self.n_in*35, 128, 2, batch_first=True)
        self.fc1 = nn.Linear(128, 512)
        self.fc2 = nn.Linear(512, 512)
        # eps branch
        self.fc1_par = nn.Linear(self.n_in*self.n_par, 32)
        self.fc2_par = nn.Linear(32, 32)
        self.fc2_par = nn.Linear(32, 32)
        # combo
        self.fc3 = nn.Linear(544, 1024)
        self.fc4 = nn.Linear(1024, self.n_out)
        # flatten LSTM parameters for memory efficiency
        self.lstm.flatten_parameters()

    def forward(self, x, par):
        # flatten the input data
        batch_size, seq_len, input_size = x.size()
        x = x.view(batch_size, -1)
        # LSTM input [batch_size, seq_len, input_size]
        # LSTM output [batch_size, seq_len, hidden_size]
        lstm_out, _ = self.lstm(x)
        # flatten the LSTM output
        lstm_out_flat = lstm_out.contiguous().view(batch_size, -1)
        # further waveform processing
        x = self.fc1(lstm_out_flat)
        x = torch_nn_functional.leaky_relu(x)
        x = self.fc2(x)
        x = torch_nn_functional.leaky_relu(x)
        # inflow of parameter
        par_x = par.view(-1, self.n_in*self.n_par)
        par_x = self.fc1_par(par_x)
        par_x = torch_nn_functional.leaky_relu(par_x)
        par_x = self.fc2_par(par_x)
        par_x = torch_nn_functional.leaky_relu(par_x)
        # combine both
        x = torch.cat([x, par_x], dim=1)
        x = self.fc3(x)
        x = torch_nn_functional.leaky_relu(x)
        x = self.fc4(x)
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
        # self.bn1 = nn.BatchNorm1d(256)
        self.fc2 = nn.Linear(2048, 2048)
        # self.bn2 = nn.BatchNorm1d(512)
        self.fc3 = nn.Linear(2048, 2048)
        # self.bn3 = nn.BatchNorm1d(256)
        self.fc4 = nn.Linear(2048, 2048)
        # self.bn4 = nn.BatchNorm1d(128)
        self.fc5 = nn.Linear(2048, 2048)
        self.fc6 = nn.Linear(2048, 1)

    def forward(self, x):
        x = self.fc1(x)
        # x = self.bn1(x)
        x = torch_nn_functional.relu(x)
        x = self.fc2(x)
        # x = self.bn2(x)
        x = torch_nn_functional.relu(x)
        x = self.fc3(x)
        # x = self.bn3(x)
        x = torch_nn_functional.relu(x)
        x = self.fc4(x)
        # x = self.bn4(x)
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
    """"
    Create Convolutional Neural Network architecture
    REQ: Required for RetrackerThresholdModel

    Note: Deprecated!
    """
    def __init__(self):
        super(TorchFunctionalWaveformModelCNN, self).__init__()
        self.cv1_1 = nn.Conv1d( 1, 32, kernel_size=(5,), stride=(1,), padding=True)
        self.cv1_2 = nn.Conv1d(32, 64, kernel_size=(5,), stride=(1,), padding=True)
        self.cv1_3 = nn.Conv1d(64, 128, kernel_size=(5,), stride=(1,), padding=True)
        self.mp1 = nn.MaxPool1d(3, stride=2)
        self.fc1 = nn.Linear(1792, 4096)
        self.fc2 = nn.Linear(4096, 1)
        # self.cv1_1 = nn.Conv1d(1, 64, kernel_size=5, stride=1, padding=True)
        # self.cv1_2 = nn.Conv1d(64, 64, kernel_size=5, stride=1, padding=True)
        # self.cv1_3 = nn.Conv1d(64, 64, kernel_size=5, stride=1, padding=True)
        # self.mp1 = nn.MaxPool1d(3, stride=2)
        # self.cv2_1 = nn.Conv1d(64, 128, kernel_size=5, stride=1, padding=True)
        # self.cv2_2 = nn.Conv1d(128, 128, kernel_size=5, stride=1, padding=True)
        # self.cv2_3 = nn.Conv1d(128, 128, kernel_size=5, stride=1, padding=True)
        # self.mp2 = nn.MaxPool1d(3, stride=2)
        # self.fc1 = nn.Linear(384, 2048)
        # self.fc2 = nn.Linear(2048, 1)

    def forward(self, x):
        x = self.cv1_1(x)
        x = torch_nn_functional.relu(x)
        x = self.cv1_2(x)
        x = torch_nn_functional.relu(x)
        x = self.cv1_3(x)
        x = torch_nn_functional.relu(x)
        x = self.mp1(x)
        # x = self.cv2_1(x)
        # x = torch_nn_functional.relu(x)
        # x = self.cv2_2(x)
        # x = torch_nn_functional.relu(x)
        # x = self.cv2_3(x)
        # x = torch_nn_functional.relu(x)
        # x = self.mp2(x)
        x = x.view(x.size(0), -1)
        x = self.fc1(x)
        x = torch_nn_functional.relu(x)
        x = self.fc2(x)
        x = torch.sigmoid(x)
        return x


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
        model_class, _ = get_cls("pysiral.auxdata.ml", torch_model_class)
        # retrieve number of input layer neurons, defaults to 45 (the Envisat standard)
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
        def id_fmi(wfm, start=0, bins=128, spacing=2):

            xshape = bins + 2 * spacing
            x = np.ndarray(shape=(xshape,))

            x[:spacing] = wfm[0] - 1.e-6
            x[-spacing:] = wfm[-1] - 1.e-6
            x[spacing:spacing + bins] = wfm

            h_b = x[start:start + bins]  # before
            start = spacing
            h_c = x[start:start + bins]  # central
            start = spacing + 1
            h_a = x[start:start + bins]  # after

            peak_candidate = np.logical_and(h_c > h_b, h_c >= h_a)
            peak_candidate = np.logical_and(peak_candidate, np.arange(bins) > 10)
            peak_candidate = np.logical_and(peak_candidate, wfm > 0.5)

            return np.where(peak_candidate)[0][0]

        # function to retrieve the sub waveform to be used
        def get_sub_wf(wfm, fmi, i0=10, i1=35):
            # zero-pad waveforms
            wfm = np.concatenate((wfm, np.zeros(50)))
            return wfm[(fmi - i0):(fmi + i1)]

        # get normalized waveform power
        waveform_norms = bn.nanmax(l1p.waveform.power, axis=1)
        normed_waveform_power = l1p.waveform.power[:] / waveform_norms[:, np.newaxis]

        # get waveform/sensor meta
        n_bins = normed_waveform_power.shape[1]
        n_wf = normed_waveform_power.shape[0]
        if n_bins >= 128:
            i0 = 10
            i1 = 35
        else:
            i0 = 10
            i1 = 25

        # subset waveform
        sub_waveform_power = np.zeros((n_wf, i0 + i1))
        c = 0
        for x in normed_waveform_power:
            try:
                wf_fmi = id_fmi(x, bins=n_bins)
            except IndexError:
                c += 1
                continue
            if wf_fmi >= 10:
                sub_waveform_power[c] = get_sub_wf(x, wf_fmi, i0, i1)
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
