# -*- coding: utf-8 -*-

"""

"""

import os
import logging

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import xarray as xr
from dataclasses import dataclass
from typing import Any, Tuple
from loguru import logger

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

    def __init__(self, cst: "SAMOSAConstants", rdb: "SAMOSARadarSpecs") -> None:
        """

        :param cst: Universal constants
        :param rdb: Database of radar altimeter specification
        """

        # flag True commands to include in the model the slope of orbit and surface
        # #(this effect usually is included in LookAngles Array)
        self.flag_slope = False

        # 1 means only one beam per resolution cell is generated in the DDM, the other ones are decimated
        self.beamsamp_factor = 1

        # flag True if the waveform under iteration is weighted
        self.wf_weighted = True

        # number of the first Look to generate in the DDM
        # (only used if LookAngles array is not passed in input: i.e. set to  None)
        self.N_Look_min = -90

        # number of the last Look to generate in the DDM
        # (only used if LookAngles array is not passed in input: i.e. set to  None)
        self.N_Look_max = 90

        # first-guess SWH in meter
        self.guess_swh = 2

        # first-guess Pu
        self.guess_pu = 1

        # first-guess nu (only used in second step of SAMOSA+)
        self.guess_nu = 2

        # lower bound on epoch in sec. If set to None, lower bound will be set to the first time in input array tau
        self.lb_epoch = None

        # lower bound on SWH in m
        self.lb_swh = -0.5

        # lower bound on Pu
        self.lb_pu = 0.2

        # lower bound on nu (only used in second step of SAMOSA+)
        self.lb_nu = 0

        # upper bound on epoch in sec. If set to None, upper bound will be set to the last time in input array tau
        self.ub_epoch = None

        # upper bound on SWH in m
        self.ub_swh = 30

        # upper bound on Pu
        self.ub_pu = 1.5

        # upper bound on nu (only used in second step of SAMOSA+)
        self.ub_nu = 1e9

        # choose between 'samosa' or 'samosa+'
        self.rtk_type = 'samosa+'

        # widening factor of PTR main lobe after Weighting Window Application
        self.wght_factor = 1.4705

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
        API method for the retracker interface in the Level-2 processor

        :param rng:
        :param wfm:
        :param indices:
        :param radar_mode:
        :param is_valid:

        :return:
        """

        # Run the retracker
        self._samosa_plus_retracker(rng, wfm, indices, radar_mode)

        # Extract ocean/lead properties

        # Extract sea ice properties

        # Filter the results
        # Needs a filter option in the config file
        # if self._options.filter.use_filter:
        #    self._filter_results()

    def _samosa_plus_retracker(self, range, wfm, indices, radar_mode) -> None:
        # Range contains the range to each bin in each waveform

        # All retracker options need to be defined in the l2 settings file
        # How to get an option
        # opt_val = self._options.name_in_file

        # Get datastructure and initialize samosa/sampy
        CST, OPT, RDB, CONF, LUT = self._get_samosa_dataclasses()
        samlib = initialize_SAMOSAlib(CST, RDB, OPT, LUT)

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
            Raw_Elevation[sin_index] = self._l1b.time_orbit.altitude[sin_index] - range[sin_index, np.shape(wfm)[1]//2]
            raw_range[sin_index] = range[sin_index, np.shape(wfm)[1]//2]
            logger.info('Using L1b range array for {:d} SARIn mode records'.format(len(sin_index)))

        for index in indices:

            LookAngles = self._get_look_angles(index)

            # Create the GEO structure needed for the samosa pacakge
            GEO = SAMOSAGeoVariables.from_l1(self._l1b, vel, hrate, ThNEcho, index)

            wf = np.array(wfm[index, :]).astype("float64")

            CONF.guess_epoch = epoch0[index]

            # Do retrack in units of Watts as we don't have the scaling factors available for calc_sigma0
            epoch_sec, swh, Pu, misfit, oceanlike_flag = samlib.Retrack_Samosa(
                tau, wf, LookAngles, MaskRanges, GEO, CONF
            )

            # SAMOSA returns a dR based upon the retracker chosen bin sampled from tau
            self._range[index] = raw_range[index] + epoch_sec * CST.c0 * 0.5

            # Compute sigma0
            sigma0, pval, cval, rval, kval = calc_sigma0(
                None, Pu, CST, RDB, GEO, epoch_sec, window_del_20_hr_ku_deuso[index],
                GEO.LAT, GEO.Height, GEO.Vs, self._l1b.classifier.transmit_power[index]
            )

            # fig, axs = plt.subplots(2)
            # axs[0].plot(samlib.last_wfm_ref, label="Waveform")
            # axs[0].plot(samlib.last_wfm_model, label="Model")
            # axs[0].legend()
            # axs[1].plot(samlib.last_wfm_ref-samlib.last_wfm_model)
            # axs[1].set_ylim(-0.25, 0.25)
            # plt.savefig(r"D:\temp\samosa\wfm_fit"+f"{index:06g}.jpg", dpi=300)
            # plt.close(fig)

            self._power[index] = sigma0

            # Store additional retracker parameters
            self.swh[index] = swh
            self.misfit[index] = misfit
            self.wind_speed[index] = func_wind_speed([sigma0])
            self.oceanlike_flag[index] = oceanlike_flag
            if not SAMOSA_DEBUG_MODE:
                self.epoch[index] = epoch_sec
                self.guess[index] = epoch0[index]
                self.Pu[index] = 65535.0 * Pu/np.max(wf)
                self.pval[index] = pval
                self.cval[index] = cval
                self.rval[index] = rval
                self.kval[index] = kval

        # Register additional auxiliary variables
        # TODO: Potentially register a lot more variables (waveforms mode, mean square slope, ...)
        self.register_auxdata_output("samswh", "samosa_swh", self.swh)
        self.register_auxdata_output("samwsp", "samosa_wind_speed", self.wind_speed)

        # Apply range bias from config files
        if "range_bias" in self._options:
            for radar_mode_index in np.arange(3):
                indices = np.where(radar_mode == radar_mode_index)[0]
                if len(indices) == 0:
                    continue
                range_bias = self._options.range_bias[radar_mode_index]
                self._range[indices] -= range_bias

        # Set the uncertainty
        if (
            "uncertainty" in self._options
            and self._options.uncertainty.type == "fixed"
        ):
            self._uncertainty[:] = self._options.uncertainty.value

        # Write debugging file
        if SAMOSA_DEBUG_MODE:
            self._samosa_debug_output(Raw_Elevation, raw_range, vel, wf_norm)

    @staticmethod
    def _get_samosa_dataclasses() -> Tuple[
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
        conf = SAMOSAConfiguration(cst, rdb)

        # Lookup table for resources filenames in the samosa package
        lut = SAMOSALookUpTables()

        return cst, opt, rdb, conf, lut

    @staticmethod
    def _compute_thermal_noise(wfm, wf_zp, n_start_noise: int = 2, n_end_noise: int = 6):
        """
        Compute thermal noise for all waveforms with samosa package

        :param wfm:
        :param wf_zp:
        :param n_start_noise: noise range gate counting from 1, no oversampling
        :param n_end_noise: noise range gate counting from 1, no oversampling
        :return:
        """
        return compute_ThNEcho(wfm, n_start_noise * wf_zp, n_end_noise * wf_zp)

    def _get_look_angles(self, index: int):
        """
        Compute the lookangles based on l1b stack information

        :param index: Waveform index

        :return:
        """

        return 90. - np.linspace(
            np.rad2deg(self._l1b.classifier.look_angle_start[index]),
            np.rad2deg(self._l1b.classifier.look_angle_stop[index]),
            num=int(self._l1b.classifier.stack_beams[index]),
            endpoint=True
        )

    @staticmethod
    def _get_normalized_waveform(wfm: npt.NDArray) -> npt.NDArray:
        """
        The samosa/SAMPy package expectes waveforms with 16bit scaling

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
                         ) -> Tuple[npt.NDArray, npt.NDArray, npt.NDArray, npt.NDArray]:

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

    def _samosa_debug_output(self, Raw_Elevation, raw_range, vel, wf_norm) -> None:
        """
        Write a netCDF file with debugging parameters

        :param Raw_Elevation:
        :param raw_range:
        :param vel:
        :param wf_norm:

        :return:
        """

        outds = xr.Dataset({'SSHunc': (['time_20_ku'], self._l1b.time_orbit.altitude-self._range),
                            'raw_elev': (['time_20_ku'], Raw_Elevation),
                            'range': (['time_20_ku'], self._range),
                            'wd_range': (['time_20_ku'], raw_range),
                            'epoch': (['time_20_ku'], self.epoch),
                            'guess': (['time_20_ku'], self.guess),
                            'Pu': (['time_20_ku'], self.Pu),
                            'lat': (['time_20_ku'], self._l1b.time_orbit.latitude),
                            'height': (['time_20_ku'], self._l1b.time_orbit.altitude),
                            'vel': (['time_20_ku'], vel),
                            'pval': (['time_20_ku'], self.pval),
                            'cval': (['time_20_ku'], self.cval),
                            'rval': (['time_20_ku'], self.rval),
                            'kval': (['time_20_ku'], self.kval),
                            'wf': (['time_20_ku', 'bins'], wf_norm),
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
