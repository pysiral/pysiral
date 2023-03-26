# -*- coding: utf-8 -*-
"""

@author: Stefan Hendricks

"""


import time

import bottleneck as bn
import numpy as np
from typing import Dict
from attrdict import AttrDict
from loguru import logger

from pysiral.core.flags import FlagContainer
from pysiral.l2proc.procsteps import Level2ProcessorStep


__all__ = ["ccilead", "corrections", "mock", "ocog", "samosa", "tfmra",
           "BaseRetracker"]


class BaseRetracker(object):
    """
    Main Retracker Class (all retrackers must be of instance BaseRetracker)
    # TODO: API clean-up is sorely needed.
    """

    def __init__(self):

        self._indices = None
        self._classifier = None
        self._l1b = None
        self._l2 = None
        self._range = None
        self._power = None
        self._options = AttrDict()

        # Dictionary containing potential auxiliary output variables of
        # the retracker algorithm that will be transferred to the l2
        # data object.
        # NOTE: This is a bit clunky because retracker algorithm will
        # not be run on all the waveforms, so a bit extra work is needed
        # that the output a) has the correct dimension and b) does not
        # overwrite itself, if the same algorithm is called consecutively.
        self.auxdata_output = []

    def set_options(self, **opt_dict):
        # TODO: Create options object, respectively use __init__
        self._options = AttrDict(opt_dict)

    def set_indices(self, indices):
        # TODO: Validation
        self._indices = indices

    def set_classifier(self, classifier):
        self._classifier = classifier

    def init(self, n_records):
        # TODO: Move to __init__
        self._create_default_properties(n_records)

    def register_auxdata_output(self, var_id, var_name, value, uncertainty=None):
        """
        Add an auxiliary parameter, that will be transferred to the l2data object after retracking
        """
        self.auxdata_output.append([var_id, var_name, value, uncertainty])

    def retrack(self, l1b, l2):

        # Store the pointer to l1b and l2 data objects
        self._l1b = l1b
        self._l2 = l2

        # Initialize the retracked range with an NaN array
        # -> only waveforms of type "surface_type" will retracked
        self._create_default_properties(l1b.n_records)

        # Give the Retracker the possibility to create additional
        # data arrays (All retracker must have this method)
        self.create_retracker_properties(l1b.n_records)

        # Set the indices for the given surface type
        # Check if there is something to do
        if self._indices is None:
            self._indices = np.arange(l1b.n_records)
        if len(self._indices) == 0:
            return False

        # Loop over each waveform of given surface type
        self.l2_retrack(l1b.waveform.range, l1b.waveform.power, self._indices,
                        l1b.waveform.radar_mode, l1b.waveform.is_valid)
        return True

    def l2_retrack(self, rng, pwr, indices, radar_mode, is_valid):
        """
        Abstract method, not to be called directly but expected to be overwritten
        by the child class
        :return:
        """
        raise NotImplementedError("BaseRetracker.l2_retrack should not be called directly")

    def get_l1b_parameter(self, data_group, parameter_name):
        """ Get any valid level-2 paremeter name """
        try:
            return self._l1b.get_parameter_by_name(data_group, parameter_name)
        except (KeyError, AttributeError):
            return None

    def get_l2_parameter(self, parameter_name):
        """ Get any valid level-2 paremeter name """
        try:
            return self._l2.get_parameter_by_name(parameter_name)
        except (KeyError, AttributeError):
            return None

    def _create_default_properties(self, n_records):
        # XXX: Currently only range and status (False: ok)
        for parameter in ["_range", "_power"]:
            setattr(self, parameter, np.full(n_records, np.nan))
        self._uncertainty = np.full(n_records, 0.0, dtype=np.float32)
        self._flag = np.zeros(shape=n_records, dtype=bool)

    def create_retracker_properties(self, n_records):
        # Will have to be overwritten
        pass

    @property
    def range(self):
        return self._range

    @property
    def uncertainty(self):
        return self._uncertainty

    @property
    def power(self):
        return self._power

    @property
    def indices(self):
        return self._indices

    @property
    def error_flag(self):
        return FlagContainer(self._flag)


class Level2RetrackerContainer(Level2ProcessorStep):
    """
    The interface for the Level-2 processor for all retrackers
    """

    def __init__(self, *args, **kwargs):
        """
        Initialize the instance
        :param args:
        :param kwargs:
        """
        super(Level2RetrackerContainer, self).__init__(*args, **kwargs)

    def execute_procstep(self, l1b, l2):
        """
        Mandatory class
        :param l1b:
        :param l2:
        :return:
        """

        # Get the error status
        error_status = self.get_clean_error_status(l2.n_records)

        # Retracker are surface type dependent
        # -> Loop over all requested surface types in the
        #    l2 processor definition file
        for surface_type, retracker_def in list(self.cfg.options.items()):

            # Check first if there are any waveforms of the requested
            # surface type
            surface_type_flag = l2.surface_type.get_by_name(surface_type)
            if surface_type_flag.num == 0:
                logger.info(f"- no waveforms of type {surface_type}")
                continue

            # Benchmark retracker performance
            timestamp = time.time()

            # Retrieve the retracker associated with surface type from the l2 settings
            retracker = get_retracker_class(retracker_def["pyclass"])

            # Set options (if any)
            if retracker_def["options"] is not None:
                retracker.set_options(**retracker_def["options"])

            # set subset of waveforms
            retracker.set_indices(surface_type_flag.indices)

            # Add classifier data (some retracker need that)
            retracker.set_classifier(l1b.classifier)

            # Start the retracker procedure
            retracker.retrack(l1b, l2)

            # Update range value of the Level-2 data object
            l2.update_retracked_range(retracker)

            # XXX: Let the retracker return other parameters?
            l2.radar_mode = l1b.waveform.radar_mode

            # retrieve potential error status and update surface type flag
            if retracker.error_flag.num > 0:
                l2.surface_type.add_flag(retracker.error_flag.flag, "invalid")
            logger.info("- Retrack class %s with %s in %.3f seconds" % (
                surface_type, retracker_def["pyclass"],
                time.time()-timestamp))

        return error_status

    @property
    def l2_input_vars(self):
        return []

    @property
    def l2_output_vars(self):
        return ["radar_mode", "range", "elev"]

    @property
    def error_bit(self):
        return self.error_flag_bit_dict["retracker"]


def wfm_get_noise_level(wfm, oversample_factor):
    """ According to CS2AWI TFMRA implementation """
    return bn.nanmean(wfm[:5*oversample_factor])


def smooth(x, window):
    """ Numpy implementation of the IDL SMOOTH function """
    return np.convolve(x, np.ones(window)/window, mode='same')


def peakdet(v, delta, x=None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html

    Returns two arrays

    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %
    %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %        in MAXTAB and MINTAB are replaced with the corresponding
    %        X-values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.

    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.

    """
    maxtab = []
    mintab = []

    if x is None:
        x = np.arange(len(v))

    v = np.asarray(v)

    mn, mx = np.Inf, -np.Inf
    mnpos, mxpos = np.NaN, np.NaN

    lookformax = True

    for i in np.arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]

        if lookformax:
            if this < mx-delta:
                maxtab.append(mxpos)
                mn = this
                mnpos = x[i]
                lookformax = False
        elif this > mn+delta:
            mintab.append(mnpos)
            mx = this
            mxpos = x[i]
            lookformax = True

    return np.array(maxtab), np.array(mintab)


def findpeaks(data, spacing=1, limit=None):
    """Finds peaks in `data` which are of `spacing` width and >=`limit`.
    :param data: values
    :param spacing: minimum spacing to the next peak (should be 1 or more)
    :param limit: peaks should have value greater or equal
    :return:
    """
    n_points = data.size
    x = np.zeros(n_points+2*spacing)
    x[:spacing] = data[0]-1.e-6
    x[-spacing:] = data[-1]-1.e-6
    x[spacing:spacing+n_points] = data
    peak_candidate = np.zeros(n_points)
    peak_candidate[:] = True
    for s in range(spacing):
        start = spacing - s - 1
        h_b = x[start:start + n_points]  # before
        start = spacing
        h_c = x[start:start + n_points]  # central
        start = spacing + s + 1
        h_a = x[start:start + n_points]  # after
        peak_candidate = np.logical_and(
            peak_candidate,
            np.logical_and(h_c > h_b, h_c > h_a)
        )

    ind = np.argwhere(peak_candidate)
    ind = ind.reshape(ind.size)
    if limit is not None:
        ind = ind[data[ind] > limit]
    return ind


def get_registered_retrackers() -> Dict:
    """
    Get a dictionary of subclasses of `BaseRetracker`
    :return: Dictionary {class_name: class}
    """
    return {cls.__name__: cls for cls in BaseRetracker.__subclasses__()}


def get_retracker_class(name: str):
    """
    Retrieve a retracker class

    :param name: retracker class name

    :raises NotImplementedError: The retracker class cannot be retrieved. An
        appropriate error message will be displayed

    :return: Initialized retracker class
    """

    from pysiral.retracker.tfmra import CYTFMRA_OK
    from pysiral.retracker.samosa import SAMOSA_OK

    if name == "cTFMRA" and not CYTFMRA_OK:
        msg = (
            "pysiral error: cTFMRA selected but pysiral.retracker.cytfmra not compiled\n"
            + "(See documentation on compilation with cython)"
        )
        raise NotImplementedError(msg)

    if name == "SAMOSAPlus" and not SAMOSA_OK:
        msg = (
            "pysiral error: SAMOSAPlus selected but not available locally\n"
            + "(Obtain SAMOSA and install via pip first)"
        )
        raise NotImplementedError(msg)

    registered_retrackers = get_registered_retrackers()
    return registered_retrackers[name]()
