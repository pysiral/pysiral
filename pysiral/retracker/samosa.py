# -*- coding: utf-8 -*-

"""

"""

import os
import logging
import typing

import matplotlib.pyplot as plt
import multiprocessing
import numpy as np
import numpy.typing as npt
import xarray as xr
from dataclasses import dataclass
from typing import Any, Tuple, List
from loguru import logger
from pysiral.l1data import L1bdataNCFile

try:
    from samosa.help_functions import calc_sigma0, func_wind_speed
    from samosa.sampy import SAMOSA as initialize_SAMOSAlib
    from samosa.sampy import compute_ThNEcho, initialize_epoch
    SAMOSA_OK = True
except ImportError:
    SAMOSA_OK = False

from pysiral import InterceptHandler
from pysiral.retracker import BaseRetracker

# TODO: Move this to an environment variable?
SAMOSA_DEBUG_MODE = False


@dataclass
class SAMOSAConstants:
    """
    Physical constants used for the SAMOSA+ Waveform model
    """
    # speed of light in m/sec
    c0: float = 299792458.
    # Reference Ellipsoid Earth Radius in m
    R_e: float = 6378137.
    # Reference Ellipsoid Earth Flatness
    f_e: float = 1 / 298.257223563
    # Gamma Function Value at 3/4
    gamma_3_4: float = 1.2254167024651779


@dataclass
class SAMOSAFittingOptions:
    """
    Settings for the waveform model fitting (scipy.optimize)
    """
    # acronym of the minimization solver, see scipy.optimize.least_squares for details
    method: str = 'trf'
    # exit tolerance on f
    ftol: float = 1e-2
    # exit tolerance on gradient norm of f
    gtol: float = 1e-2
    # exit tolerance on x
    xtol: float = 2 * 1e-3
    # relative step size for the finite difference approximation of the Jacobian
    diff_step: Any = None
    # maximum number of function evaluations
    max_nfev: Any = None
    # loss function , see scipy.optimize.least_squares for details
    loss: str = 'linear'


@dataclass
class SAMOSALookUpTables:
    """
    Links to the lookup-table filenames in the samosa package
    """

    # filename of the F0 LUT
    F0: str = 'LUT_F0.txt'

    # filename of the F1 LUT
    F1: str = 'LUT_F1.txt'

    # filename of the alphap LUT (case no weighting)
    alphap_noweight: str = 'alphap_table_DX3000_ZP20_SWH20_10_Sept_2019(CS2_NOHAMMING).txt'

    # filename of the alphap LUT (case weighting)
    alphap_weight: str = 'alphap_table_DX3000_ZP20_SWH20_10_Sept_2019(CS2_HAMMING).txt'

    # filename of the alpha power LUT (case no weighting)
    alphapower_noweight: str = 'alphaPower_table_CONSTANT_SWH20_10_Feb_2020(CS2_NOHAMMING).txt'

    # filename of the alpha power LUT (case weighting)
    alphapower_weight: str = 'alphaPower_table_CONSTANT_SWH20_10_Feb_2020(CS2_NOHAMMING).txt'


class SAMOSARadarSpecs:
    """
    Radar altimeter specifications needed for the SAMOSA(+) Retracker
    """

    def __init__(self,
                 Np_burst: int = None,
                 Npulse: int  = None,
                 PRF_SAR: float = None,
                 BRI: float = None,
                 f_0: float = None,
                 Bs: float = None,
                 theta_3x: float = None,
                 theta_3y: float = None,
                 G_0: float = None,
                 bias_sigma0: float = None
    ) -> None:
        """
        Radar altimeter specifications for the generating the waveform model.

        :param Np_burst: number of pulses per burst
        :param Npulse: number of the range gates per pulse (without zero-padding)
        :param PRF_SAR: Pulse Repetition Frequency in SAR mode, given in Hz
        :param BRI: Burst Repetition Interval, given in sec
        :param f_0: Carrier Frequency in Hz
        :param Bs: Sampled Bandwidth in Hz
        :param theta_3x: (rad) Antenna 3 dB beam-width (along-track)
        :param theta_3y:  (rad) Antenna 3 dB beamwidth (cross-track)
        :param G_0: Boresight One-Way Antenna Power Radiation Gain (natural units)
        :param bias_sigma0: static bias in sigma0 (dB)
        """

        self.Np_burst = Np_burst
        self.Npulse = Npulse
        self.PRF_SAR = PRF_SAR
        self.BRI = BRI
        self.f_0 = f_0
        self.Bs = Bs
        self.theta_3x = theta_3x
        self.theta_3y = theta_3y
        self.G_0 = G_0
        self.bias_sigma0 = bias_sigma0

    @classmethod
    def from_preset(cls, preset: str) -> "SAMOSARadarSpecs":
        """
        Sets the radar altimeter specifications for the given presets

        :param preset: The name of the preset of type {platform}_{sensor}_{mode}

        :raise ValueError: Invalid preset

        :return: Initialized SAMOSARadarSpecs instance
        """

        if preset == "cryosat2_siral_sar":
            kwargs = {
                # number of pulses per burst
                "Np_burst": 64,
                # number of the range gates per pulse (without zero-padding)
                "Npulse": 128,
                # Pulse Repetition Frequency in SAR mode, given in Hz
                "PRF_SAR": 18181.8181818181,
                # Burst Repetition Interval, given in sec
                "BRI": 0.0117929625,
                # Carrier Frequency in Hz
                "f_0": 13.575e9,
                # Sampled Bandwidth in Hz
                "Bs": 320e6,
                # (rad) Antenna 3 dB beam-width (along-track)
                "theta_3x": np.deg2rad(1.10),
                # (rad) Antenna 3 dB beamwidth (cross-track)
                "theta_3y": np.deg2rad(1.22),
                # Boresight One-Way Antenna Power Radiation Gain (natural units)
                "G_0": 10. ** (42.6 / 10),
                # static bias in sigma0 (dB)
                "bias_sigma0": -3.04
            }
        else:
            raise ValueError(f"Unknown preset name {preset} [cryosat2_siral_sar]")

        return cls(**kwargs)


@dataclass
class SAMOSAFitResult:
    """
    Container for output of the SAMOSA+ retracker
    """
    tau: np.ndarray
    wf: np.ndarray
    wf_model: np.ndarray
    epoch_sec: float
    rng: float
    nu: float
    swh: float
    Pu: float
    misfit: float
    oceanlike_flag: bool
    sigma0: float
    pval: float
    cval: float
    rval: float
    kval: float

    def debug_plot(self):

        plt.figure(figsize=(8, 6), dpi=150)
        plt.plot(self.tau, self.wf, "-o", color="#e08214", ms=2, lw=0.5)
        plt.plot(self.tau, self.wf_model, color="#542788", lw=0.75, alpha=0.8)
        plt.axvline(self.epoch_sec, color="black", lw=0.5, ls="dashed")
        plt.annotate(f"misfit = {self.misfit:.02f}", (0.5, 0.9), xycoords="axes fraction")
        plt.annotate(f"surface height sdev (swh/4) = {self.swh/4.:.02f}m", (0.5, 0.85),
                     xycoords="axes fraction")
        plt.annotate(r"mean square slope (1/$\nu$)" + f" = {1./self.nu:.2E}", (0.5, 0.8),
                     xycoords="axes fraction")
        plt.xlabel(r"$\tau$ (sec)")
        plt.ylabel("Normed Waveform Power")
        plt.show()


class SAMOSAGeoVariables(object):

    def __init__(self,
                 lat: float,
                 lon: float,
                 height: float,
                 vs: float,
                 hrate: float,
                 pitch: float,
                 roll: float,
                 nu: float,
                 track_sign: int,
                 thn: float
                 ) -> None:
        """


        :param lat: latitude in degree for the waveform under iteration
        :param lon: longitude in degree between -180, 180 for the waveform under iteration
        :param height: Orbit Height in meter for the waveform under iteration
        :param vs:  Satellite Velocity in m/s for the waveform under iteration
        :param hrate:  Orbit Height rate in m/s for the waveform under iteration
        :param pitch: Altimeter Reference Frame Pitch in radiant
        :param roll: Altimeter Reference Frame Roll in radiant
        :param nu: Inverse of the mean square slope
        :param track_sign: if Track Ascending => -1, if Track Descending => +1, set it to zero if flag_slope=False in CONF
        :param thn: thermal noise
        """

        self.LAT = lat
        self.LON = lon
        self.Height = height
        self.Vs = vs
        self.Hrate = hrate
        self.Pitch = pitch
        self.Roll = roll
        self.nu = nu
        self.track_sign = track_sign
        self.ThN = thn

    @classmethod
    def from_l1(cls, l1b, vel, hrate, ThNEcho, index) -> "SAMOSAGeoVariables":
        """

        :param l1b:
        :param vel:
        :param ThNEcho:
        :param index:

        :return:
        """
        return cls(
            l1b.time_orbit.latitude[index],
            l1b.time_orbit.longitude[index],
            l1b.time_orbit.altitude[index],
            np.squeeze(vel)[index],
            hrate[index],
            np.radians(l1b.time_orbit.antenna_pitch[index]),
            np.radians(l1b.time_orbit.antenna_roll[index]),
            0,
            0,
            np.squeeze(ThNEcho)[index]
        )


class SAMOSAConfiguration(object):
    """
    Configuration settings for the SAMOSA+ retracker
    """

    def __init__(self,
                 cst: "SAMOSAConstants",
                 rdb: "SAMOSARadarSpecs",
                 flag_slope: bool = False,
                 beamsamp_factor: int = 1,
                 wf_weighted: bool = True,
                 N_Look_min: int = -90,
                 N_Look_max: int = 90,
                 guess_swh: float = 2.,
                 guess_pu: float =  1.,
                 guess_nu: float =  2.,
                 lb_epoch: float =  None,
                 lb_swh: float = -0.5,
                 lb_pu: float = 0.2,
                 lb_nu: float = 0.0,
                 ub_epoch: float = None,
                 ub_swh: float = 30.,
                 ub_pu: float = 1.5,
                 ub_nu: float = 1e9,
                 rtk_type: str = 'samosa+',
                 wght_factor: float = 1.4705
                 ) -> None:
        """

        :param cst: Universal constants
        :param rdb: Database of radar altimeter specification
        """

        # flag True commands to include in the model the slope of orbit and surface
        # #(this effect usually is included in LookAngles Array)
        self.flag_slope = flag_slope

        # 1 means only one beam per resolution cell is generated in the DDM, the other ones are decimated
        self.beamsamp_factor = beamsamp_factor

        # flag True if the waveform under iteration is weighted
        self.wf_weighted = wf_weighted

        # number of the first Look to generate in the DDM
        # (only used if LookAngles array is not passed in input: i.e. set to  None)
        self.N_Look_min = N_Look_min

        # number of the last Look to generate in the DDM
        # (only used if LookAngles array is not passed in input: i.e. set to  None)
        self.N_Look_max = N_Look_max

        # first-guess SWH in meter
        self.guess_swh = guess_swh

        # first-guess Pu
        self.guess_pu = guess_pu

        # first-guess nu (only used in second step of SAMOSA+)
        self.guess_nu = guess_nu

        # lower bound on epoch in sec. If set to None, lower bound will be set to the first time in input array tau
        self.lb_epoch = lb_epoch

        # lower bound on SWH in m
        self.lb_swh = lb_swh

        # lower bound on Pu
        self.lb_pu = lb_pu

        # lower bound on nu (only used in second step of SAMOSA+)
        self.lb_nu = lb_nu

        # upper bound on epoch in sec. If set to None, upper bound will be set to the last time in input array tau
        self.ub_epoch = ub_epoch

        # upper bound on SWH in m
        self.ub_swh = ub_swh

        # upper bound on Pu
        self.ub_pu = ub_pu

        # upper bound on nu (only used in second step of SAMOSA+)
        self.ub_nu = ub_nu

        # choose between 'samosa' or 'samosa+'
        self.rtk_type = rtk_type

        # widening factor of PTR main lobe after Weighting Window Application
        self.wght_factor = wght_factor

        self.Lz = cst.c0 / (2. * rdb.Bs)


class SAMOSAPlus(BaseRetracker):
    """
    Interface to the SAMOSA+ retracker by CLS.
    Retracker must be installed as a package into the environment for this class to be used.

    """

    def __init__(self):
        """
        Initialize the SAMOSA+ retracker
        """
        super(SAMOSAPlus, self).__init__()

        # Check imports
        # NOTE: The samosa package is not installed with pysiral on default.
        if SAMOSA_OK:
            logger.info("SAMOSA retracker loaded from the environment")
            logging.getLogger('samosa.sampy').addHandler(InterceptHandler())
            logging.getLogger('samosa.sampy').setLevel('INFO')
        else:
            logger.error("Unable to import the SAMOSA retracker. Has it been installed?")

        # Set a dictionary with additions retracker parameters
        self._retracker_params = {}

    def create_retracker_properties(self, n_records: int) -> None:
        """
        Initialize retracker properties with correct arrays shapes (shape = (n_records, )).
        The list of parameters depends on whether the SAMOSA_DEBUG_MODE flag is set.

        NOTE: The properties are set to an array, but can be accessed as `self.{property_name}`
        via the __getattr__ method.

        :param n_records:
        """

        parameter = (
            ["misfit", "swh", "wind_speed", "oceanlike_flag"]
            if SAMOSA_DEBUG_MODE
            else [
                "misfit",
                "swh",
                "mean_square_slope",
                "wind_speed",
                "oceanlike_flag",
                "epoch",
                "guess",
                "Pu",
                "rval",
                "kval",
                "pval",
                "cval",
            ]
        )
        for parameter_name in parameter:
            self._retracker_params[parameter_name] = np.full(n_records, np.nan, dtype=np.float32)

    def l2_retrack(self, rng, wfm, indices, radar_mode, is_valid) -> None:
        """
        API method for the retracker interface in the Level-2 processor.

        :param rng: The range per waveform samples for all waveforms in the Level-1 data object
        :param wfm: All waveforms in the Level-1 data object
        :param indices: List of waveforms for the retracker
        :param radar_mode: Flag indicating the radar mode for all waveforms
        :param is_valid: Error flag for all waveforms

        :return: None (Output is added to the instance)
        """

        # Run the retracker
        # NOTE: Output is a SAMOSAFitResult dataclass for each index in indices.
        fit_results = self._samosa_plus_retracker(rng, wfm, indices, radar_mode)

        # [fit_result.debug_plot() for fit_result in fit_results]
        # breakpoint()

        # Store retracker properties (including range)
        self._store_retracker_properties(fit_results, indices)

        # Set/compute uncertainty
        self._set_range_uncertainty()

        # Add range biases (when set in config file)
        self._set_range_bias(radar_mode)

        # Add auxiliary variables to the l2 data object
        self._register_auxiliary_variables()

        if SAMOSA_DEBUG_MODE:
            self._samosa_debug_output()

    def _store_retracker_properties(self, fit_results, indices) -> None:
        """
        Store the output of the SAMOSA+ retracker in the class.

        :param fit_results:
        :param indices:

        :return:
        """

        for index, fit_result in zip(indices, fit_results):

            self._range[index] = fit_result.rng
            self._power[index] = fit_result.sigma0

            # Store additional retracker parameters
            self.swh[index] = fit_result.swh
            self.misfit[index] = fit_result.misfit
            self.wind_speed[index] = func_wind_speed([fit_result.sigma0])
            self.oceanlike_flag[index] = fit_result.oceanlike_flag
            if fit_result.nu == 0:
                self.mean_square_slope[index] = np.nan
            else:
                self.mean_square_slope[index] = 1. / fit_result.nu
            if not SAMOSA_DEBUG_MODE:
                self.epoch[index] = fit_result.epoch_sec
                self.guess[index] = self._retracker_params["epoch0"][index]
                self.Pu[index] = 65535.0 * fit_result.Pu/np.max(fit_result.wf)
                self.pval[index] = fit_result.pval
                self.cval[index] = fit_result.cval
                self.rval[index] = fit_result.rval
                self.kval[index] = fit_result.kval

    def _set_range_bias(self, radar_mode) -> None:
        """
        Set range bias bases on radar mode (from config file)

        :param radar_mode:

        :return:
        """

        # Apply range bias from config files
        if "range_bias" in self._options:
            for radar_mode_index in np.arange(3):
                indices = np.where(radar_mode == radar_mode_index)[0]
                if len(indices) == 0:
                    continue
                range_bias = self._options.range_bias[radar_mode_index]
                self._range[indices] -= range_bias

    def _set_range_uncertainty(self) -> None:
        """
        Estimate the uncertainty of the SAMOSA+ retracker result.
        At the moment there is no other implementation than using
        a fixed uncertainty value from the Level-2 processor
        definition file.
        """

        # Set the uncertainty
        if (
            "uncertainty" in self._options
            and self._options.uncertainty.type == "fixed"
        ):
            self._uncertainty[:] = self._options.uncertainty.value
        else:
            logger.warning("No uncertainty definition for SAMOSA+ range")

    def _register_auxiliary_variables(self) -> None:
        """
        Add auxiliary variables to the L2 data object. The specific variables
        depend on surface type.
        """

        # General auxiliary variables
        self.register_auxdata_output("sammf", "samosa_misfit", self.misfit)
        self.register_auxdata_output("sammss", "samosa_mean_square_slope", self.mean_square_slope)

        # Lead and open ocean surfaces
        surface_type = self._options.get("surface_type", "undefined")
        if surface_type == "polar_ocean":
            self.register_auxdata_output("samswh", "samosa_swh", self.swh)
            self.register_auxdata_output("samwsp", "samosa_wind_speed", self.wind_speed)

        # Sea ice surface types
        elif surface_type == "sea_ice":
            self.register_auxdata_output("samshsd", "samosa_surface_height_standard_deviation", self.swh / 4.)

        else:
            logger.warning("No specific surface type set for SAMOSA+: No auxiliary variables added to L2")

    def _samosa_plus_retracker(self, rng, wfm, indices, radar_mode) -> List[SAMOSAFitResult]:
        """
        Run the SAMOSA+ retracker for a set of waveforms.

        NOTE: The entire data structure will likely change.

        :param rng:
        :param wfm:
        :param indices:
        :param radar_mode:
        :return:
        """
        # Range contains the range to each bin in each waveform

        # All retracker options need to be defined in the l2 settings file
        # How to get an option
        # opt_val = self._options.name_in_file

        # Get datastructure and initialize samosa/sampy
        CST, OPT, RDB, CONF, LUT = self._get_samosa_dataclasses()
        sampy_kwargs = self._options.get("sampy_kwargs", {})
        samlib = initialize_SAMOSAlib(CST, RDB, OPT, LUT, **sampy_kwargs)

        # Further settings (l2 processor options?)
        MaskRanges = None

        # Prepare input data
        # TODO: Potentially move to data class ?
        window_del_20_hr_ku_deuso = self._l1b.classifier.window_delay
        tau, raw_range, Raw_Elevation, wf_zp = self._get_range_array(wfm, RDB, CST)
        wf_norm = self._get_normalized_waveform(wfm)
        hrate, vel = self._get_altitude_velocity_from_l1()
        # TODO: Move n_start_noise, n_end_noise to config file
        ThNEcho = self._compute_thermal_noise(wfm.T, wf_zp)
        # initializing the epoch (first-guess epoch) from the waveform matrix
        epoch0 = initialize_epoch(wf_norm.T, tau, Raw_Elevation, CST, size_half_block=10)

        sin_index = np.where(radar_mode == 2)[0]
        if len(sin_index > 0):
            Raw_Elevation[sin_index] = self._l1b.time_orbit.altitude[sin_index] - rng[sin_index, np.shape(wfm)[1]//2]
            raw_range[sin_index] = rng[sin_index, np.shape(wfm)[1]//2]
            logger.info('Using L1b range array for {:d} SARIn mode records'.format(len(sin_index)))

        # for Write debugging file
        self._retracker_params["epoch0"] = epoch0
        if SAMOSA_DEBUG_MODE:
            self._retracker_params["raw_elevation"] = Raw_Elevation
            self._retracker_params["raw_range"] = raw_range
            self._retracker_params["vel"] = vel
            self._retracker_params["wf_norm"] = wf_norm

        # Fit SAMOSA+ waveform model and return list of results
        args = [samlib, self._l1b, wfm, tau, window_del_20_hr_ku_deuso, vel, hrate, ThNEcho,
                CONF, epoch0, MaskRanges, raw_range, CST, RDB]
        return [
            fit_samosa_waveform_model(index, *args) for index in indices
        ]

        # TODO Investigate speeding multi-processing
        # The use of a multiprocessing Pool is currently slower (!) than just using
        # a single thread. The reason might the work overhead for creating each process.
        # One idea would be to divide the indices in x chunks and not to start a process
        # for each waveform.

        # n_processes = multiprocessing.cpu_count()
        # logger.info(f"Use multi-processing with {n_processes} workers")
        # process_pool = multiprocessing.Pool(n_processes)
        # jobs = [process_pool.apply_async(fit_samosa_waveform_model, args=(index, *args)) for index in indices]
        # process_pool.close()
        # process_pool.join()
        # fit_results = [job.get() for job in jobs]

    def _get_samosa_dataclasses(self) -> Tuple[
        "SAMOSAConstants",
        "SAMOSAFittingOptions",
        "SAMOSARadarSpecs",
        "SAMOSAConfiguration",
        "SAMOSALookUpTables"
    ]:
        """
        Get the data classes requires for the SAMOSA retracker

        :return:
        """

        # Universal constants
        cst = SAMOSAConstants()

        # minimization scheme settings
        opt = SAMOSAFittingOptions()

        # Radar altimeter specifications
        rdb = SAMOSARadarSpecs.from_preset("cryosat2_siral_sar")

        # SAMOSA retracker configuration
        sampy_conf_kwargs = self._options.get("sampy_conf_kwargs", {})
        conf = SAMOSAConfiguration(cst, rdb, **sampy_conf_kwargs)

        # Lookup table for resources filenames in the samosa package
        lut = SAMOSALookUpTables()

        return cst, opt, rdb, conf, lut

    @staticmethod
    def _compute_thermal_noise(
            wfm: npt.NDArray,
            wf_zp: int,
            n_start_noise: int = 2,
            n_end_noise: int = 6
    ) -> float:
        """
        Compute thermal noise for all waveforms with samosa/sampy package from the
        early range gated in the range window.

        :param wfm: waveform power
        :param wf_zp: zero padding waveform oversampling factor
        :param n_start_noise: noise range gate counting from 1, no oversampling
        :param n_end_noise: noise range gate counting from 1, no oversampling

        :return: Thermal noise as computed by the `sampy` package
        """
        return compute_ThNEcho(wfm, n_start_noise * wf_zp, n_end_noise * wf_zp)

    @staticmethod
    def _get_normalized_waveform(wfm: npt.NDArray) -> npt.NDArray:
        """
        The samosa/SAMPy package expecte waveforms with 16bit scaling. Scale
        any waveform to this range.

        :param wfm:

        :return:
        """
        wf_norm = np.zeros_like(wfm)
        for rec in np.arange(np.shape(wfm)[0]):
            wf_norm[rec, :] = (65535.0 * wfm[rec, :] / np.nanmax(wfm[rec, :])).round().astype(np.uint16)
        return wf_norm

    def _get_range_array(self,
                         wfm,
                         rdb: "SAMOSARadarSpecs",
                         cst: "SAMOSAConstants"
                         ) -> Tuple[npt.NDArray, npt.NDArray, npt.NDArray, float]:

        # zero-padding factor of the waveform
        wf_zp = np.shape(wfm)[1] / rdb.Npulse
        logger.info('Waveform zero padding factor is {:f}'.format(wf_zp))
        Nstart = rdb.Npulse * wf_zp
        Nend = rdb.Npulse * wf_zp
        # time sampling step for the array tau, it includes the zero-padding factor
        dt = 1. / (rdb.Bs * wf_zp)

        # print(np.arange(-(4 / 2), ((4 - 1) / 2)))
        # [-2. -1.  0.  1.]
        # Zero bin as at len()//2. So use the range at this location as the base to adjust from
        tau = np.arange(-(Nstart / 2) * dt, ((Nend - 1) / 2) * dt, dt)

        # window_del_20_hr_ku_deuso = window_del_20_hr_ku * (uso_cor_20_hr_ku + 1)
        window_del_20_hr_ku_deuso = self._l1b.classifier.window_delay
        # Raw_Elevation = alt_20_hr_ku - cst.c0 / 2 * window_del_20_hr_ku_deuso
        raw_range = cst.c0 * window_del_20_hr_ku_deuso * 0.5
        raw_elevation = self._l1b.time_orbit.altitude - raw_range

        return tau, raw_range, raw_elevation, wf_zp

    def _get_altitude_velocity_from_l1(self) -> Tuple[npt.NDArray, npt.NDArray]:
        """
        Get altitude (height) rate and satellite velocity from l1 data

        :return:
        """

        hrate = np.zeros_like(self._l1b.time_orbit.altitude_rate)
        vel = np.sqrt(self.get_l1b_parameter("classifier", "satellite_velocity_x")**2
                      + self.get_l1b_parameter("classifier", "satellite_velocity_y")**2
                      + self.get_l1b_parameter("classifier", "satellite_velocity_z")**2)
        return hrate, vel

    def _samosa_debug_output(self) -> None:
        """
        Write a netCDF file with debugging parameters

        :return:
        """

        outds = xr.Dataset({'SSHunc': (['time_20_ku'], self._l1b.time_orbit.altitude-self._range),
                            'raw_elev': (['time_20_ku'], self._retracker_params["raw_elevation"]),
                            'range': (['time_20_ku'], self._range),
                            'wd_range': (['time_20_ku'], self._retracker_params["raw_range"]),
                            'epoch': (['time_20_ku'], self.epoch),
                            'guess': (['time_20_ku'], self.guess),
                            'Pu': (['time_20_ku'], self.Pu),
                            'lat': (['time_20_ku'], self._l1b.time_orbit.latitude),
                            'height': (['time_20_ku'], self._l1b.time_orbit.altitude),
                            'vel': (['time_20_ku'], self._retracker_params["vel"]),
                            'pval': (['time_20_ku'], self.pval),
                            'cval': (['time_20_ku'], self.cval),
                            'rval': (['time_20_ku'], self.rval),
                            'kval': (['time_20_ku'], self.kval),
                            'wf': (['time_20_ku', 'bins'], self._retracker_params["wf_norm"]),
                            'misfit': (['time_20_ku'], self.misfit),
                            'oceanlike_flag': (['time_20_ku'], self.oceanlike_flag),
                            'SWH': (['time_20_ku'], self.swh),
                            'sigma0': (['Sigma0_20Hz'], self._power),
                            'wind_speed': (['U10_20Hz'], self.wind_speed)},
                           coords={'time_20_ku': self._l1b.time_orbit.timestamp,
                                   'bins': np.arange(256),
                                   'lon_20_ku': (['time_20_ku'], self._l1b.time_orbit.longitude),
                                   'lat_20_ku': (['time_20_ku'], self._l1b.time_orbit.latitude)},
                           attrs={'description': "Parameters from SAMOSA+ retracker"})
        outds.time_20_ku.encoding = {'calendar': 'gregorian',
                                     'units': 'seconds since 2000-01-01 0:0:0'}
        if os.path.isfile('samosa_debug.nc'):
            os.remove('samosa_debug.nc')
        outds.to_netcdf('samosa_debug.nc')

    def _filter_results(self) -> None:
        """
        Nothing here yet.
        :return:
        """
        pass

    def __getattr__(self, item: str) -> Any:
        """
        Direct attribute access to the retracker properties dictionary

        :param item: parameter name

        :raise AttributeError: item not in self._retracker_params (see self.create_retracker_properties)

        :return:
        """
        if item in self._retracker_params:
            return self._retracker_params[item]
        else:
            raise AttributeError(f"{self.__class__.__name__} has no attribute {item}")


def fit_samosa_waveform_model(
        index, samlib, l1b, wfm, tau, window_del_20_hr_ku_deuso, vel, hrate,
        ThNEcho, CONF, epoch0, MaskRanges, raw_range, CST, RDB
):
    """
    Fitting procedure for one waveform as function

    :param index:
    :param samlib:
    :param wfm:
    :param tau:
    :param window_del_20_hr_ku_deuso:
    :param vel:
    :param hrate:
    :param ThNEcho:
    :param CONF:
    :param epoch0:
    :param MaskRanges:
    :param raw_range:
    :param CST:
    :param RDB:
    :return:
    """

    look_angles = get_look_angles(l1b, index)

    # Create the GEO structure needed for the samosa pacakge
    GEO = SAMOSAGeoVariables.from_l1(l1b, vel, hrate, ThNEcho, index)

    wf = np.array(wfm[index, :]).astype("float64")

    CONF.guess_epoch = epoch0[index]

    # Do retrack in units of Watts as we don't have the scaling factors available for calc_sigma0
    epoch_sec, swh, Pu, misfit, oceanlike_flag = samlib.Retrack_Samosa(
        tau, wf, look_angles, MaskRanges, GEO, CONF
    )

    # Extract additional parameters from samlib instance
    wf_model = samlib.last_wfm_model
    wf_ref = samlib.last_wfm_ref
    nu = samlib.last_nu

    # SAMOSA returns a dR based upon the retracker chosen bin sampled from tau
    rng = raw_range[index] + epoch_sec * CST.c0 * 0.5

    # Compute sigma0
    sigma0, pval, cval, rval, kval = calc_sigma0(
        None, Pu, CST, RDB, GEO, epoch_sec, window_del_20_hr_ku_deuso[index],
        GEO.LAT, GEO.Height, GEO.Vs, l1b.classifier.transmit_power[index]
    )

    var = [tau, wf_ref, wf_model, epoch_sec, rng, nu, swh, Pu, misfit, oceanlike_flag, sigma0, pval, cval, rval, kval]
    return SAMOSAFitResult(*var)


def get_look_angles(l1b: L1bdataNCFile, index: int) -> npt.NDArray:
    """
    Compute the lookangles based on l1b stack information
    # TODO: This functions raises a ValueError for NaN values in the classifiers (LRM?)

    :param l1b: The level 1 data object
    :param index: Waveform index

    :return:
    """
    return 90. - np.linspace(
        np.rad2deg(l1b.classifier.look_angle_start[index]),
        np.rad2deg(l1b.classifier.look_angle_stop[index]),
        num=int(l1b.classifier.stack_beams[index]),
        endpoint=True
    )
