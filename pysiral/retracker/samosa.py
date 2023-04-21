# -*- coding: utf-8 -*-

"""

"""

import logging

import numpy as np
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


class SAMOSAPlus(BaseRetracker):
    """
    Interface to the SAMOSA+ retracker by CLS.
    Retracker must be installed as a package into the environment for this class to be used.

    """

    def __init__(self):
        super(SAMOSAPlus, self).__init__()
        if SAMOSA_OK:
            logger.info("SAMOSA retracker loaded from the environment")
            logging.getLogger('samosa.sampy').addHandler(InterceptHandler())
            logging.getLogger('samosa.sampy').setLevel('INFO')
        else:
            logger.error("Unable to import the SAMOSA retracker. Has it been installed?")

    def create_retracker_properties(self, n_records):
        # False branches here and below were used for debugging
        if True:
            parameter = ["misfit", "swh", "wind_speed", "oceanlike_flag", "epoch", "guess", "Pu", "rval", "kval", "pval", "cval"]
        else:
            parameter = ["misfit", "swh", "wind_speed", "oceanlike_flag"]
        for parameter_name in parameter:
            setattr(self, parameter_name,
                    np.ndarray(shape=(n_records), dtype=np.float32) * np.nan)

    def l2_retrack(self, range, wfm, indices, radar_mode, is_valid):
        # Run the retracker
        self._samosa_plus_retracker(range, wfm, indices, radar_mode)
        # Filter the results
        # Needs a filter option in the config file
        #if self._options.filter.use_filter:
        #    self._filter_results()

    def _samosa_plus_retracker(self, range, wfm, indices, radar_mode):
        # Range contains the range to each bin in each waveform

        # All retracker options need to be defined in the l2 settings file
        # How to get an option
        # opt_val = self._options.name_in_file
        ### CST is a structure collecting universal constants

        CST = type('', (), {})()

        CST.c0 = 299792458.  ## speed of light in m/sec
        CST.R_e = 6378137.  ## Reference Ellipsoid Earh Radius in m
        CST.f_e = 1 / 298.257223563  ## Reference Ellipsoid Earth Flatness
        CST.gamma_3_4 = 1.2254167024651779  ## Gamma Function Value at 3/4

        ### OPT is a structure collecting parameters relative to the minimization scheme settings

        OPT = type('', (), {})()

        OPT.method = 'trf'  ## acronym of the minimization solver, see scipy.optimize.least_squares for details
        OPT.ftol = 1e-2  ## exit tolerance on f
        OPT.gtol = 1e-2  ## exit tolerance on gradient norm of f
        OPT.xtol = 2 * 1e-3  ## exit tolerance on x
        OPT.diff_step = None  ## relative step size for the finite difference approximation of the Jacobian
        OPT.max_nfev = None  ## maximum number of function evaluations
        OPT.loss = 'linear'  ## loss function , see scipy.optimize.least_squares for details

        ### RDB is a structure collecting parameters relative to the sensor radar database

        RDB = type('', (), {})()

        RDB.Np_burst = 64  # number of pulses per burst
        RDB.Npulse = 128  # number of the range gates per pulse (without zero-padding)
        RDB.PRF_SAR = 18181.8181818181  # Pulse Repetition Frequency in SAR mode , given in Hz
        RDB.BRI = 0.0117929625  # Burst Repetition Interval, given in sec
        RDB.f_0 = 13.575e9  # Carrier Frequency in Hz
        RDB.Bs = 320e6  # Sampled Bandwidth in Hz
        RDB.theta_3x = np.deg2rad(1.10)  # (rad) Antenna 3 dB beamwidth (along-track)
        RDB.theta_3y = np.deg2rad(1.22)  # (rad) Antenna 3 dB beamwidth (cross-track)
        RDB.G_0 = 10. ** (42.6 / 10)  # Boresight One-Way Antenna Power Radiation Gain (natural units)
        RDB.bias_sigma0 = -3.04  # static bias in sigma0 (dB)

        ### LUT is a structure collecting filenames of the all SAMOSA LUT
        ### All the LUT files must be in a folder named auxi and located in the same folder as sampy.py

        # FIXME: Need to add a configuration parameter that gives a path to these files. Per surface type. Possibly
        # per mode as well?
        LUT = type('', (), {})()

        LUT.F0 = 'LUT_F0.txt'  ## filename of the F0 LUT
        LUT.F1 = 'LUT_F1.txt'  ## filename of the F1 LUT
        LUT.alphap_noweight = 'alphap_table_DX3000_ZP20_SWH20_10_Sept_2019(CS2_NOHAMMING).txt'  ## filename of the alphap LUT ( case no weighting)
        LUT.alphap_weight = 'alphap_table_DX3000_ZP20_SWH20_10_Sept_2019(CS2_HAMMING).txt'  ## filename of the alphap LUT ( case weighting)
        LUT.alphapower_noweight = 'alphaPower_table_CONSTANT_SWH20_10_Feb_2020(CS2_NOHAMMING).txt'  ## filename of the alpha power LUT ( case no weighting)
        LUT.alphapower_weight = 'alphaPower_table_CONSTANT_SWH20_10_Feb_2020(CS2_NOHAMMING).txt'  ## filename of the alpha power LUT ( case weighting)

        ### time array tau : it gives the relative time of each range gate of the radar waveform with respect a time zero
        ### time zero corresponds at the time of the reference gate

        wf_zp = np.shape(wfm)[1] / RDB.Npulse  #### zero-padding factor of the waveform
        logger.info('Waveform zero padding factor is {:f}'.format(wf_zp))
        Nstart = RDB.Npulse * wf_zp
        Nend = RDB.Npulse * wf_zp
        dt = 1. / (RDB.Bs * wf_zp)  #### time sampling step for the array tau, it includes the zero-padding factor

        # print(np.arange(-(4 / 2), ((4 - 1) / 2)))
        # [-2. -1.  0.  1.]
        # Zero bin as at len()//2. So use the range at this location as the base to adjust from
        tau = np.arange(-(Nstart / 2) * dt, ((Nend - 1) / 2) * dt, dt)

        NstartNoise = 2  ## noise range gate counting from 1, no oversampling
        NendNoise = 6  ## noise range gate counting from 1, no oversampling

        #window_del_20_hr_ku_deuso = window_del_20_hr_ku * (uso_cor_20_hr_ku + 1)
        window_del_20_hr_ku_deuso = self._l1b.classifier.window_delay
        #Raw_Elevation = alt_20_hr_ku - CST.c0 / 2 * window_del_20_hr_ku_deuso
        raw_range = CST.c0 * window_del_20_hr_ku_deuso * 0.5
        Raw_Elevation = self._l1b.time_orbit.altitude - raw_range

        sin_index = np.where(radar_mode == 2)[0]
        if len(sin_index > 0):
            Raw_Elevation[sin_index] = self._l1b.time_orbit.altitude[sin_index] - range[sin_index,np.shape(wfm)[1]//2]
            raw_range[sin_index] = range[sin_index,np.shape(wfm)[1]//2]
            logger.info('Using L1b range array for {:d} SARIn mode records'.format(len(sin_index)))

        ThNEcho = compute_ThNEcho(wfm.T, NstartNoise * wf_zp,
                                  NendNoise * wf_zp)  ### computing Thermal Noise from the waveform matric

        # initialize_epoch relies on the waveform being in counts, not watts, so revert
        wf_norm = np.zeros_like(wfm)
        for rec in np.arange(np.shape(wfm)[0]):
            wf_norm[rec,:] = (65535.0 * wfm[rec,:] / np.nanmax(wfm[rec,:])).round().astype(np.uint16)

        epoch0 = initialize_epoch(wf_norm.T, tau, Raw_Elevation, CST,
                                  size_half_block=10)  ### initializing the epoch (first-guess epoch) from the waveform matrix

        samlib = initialize_SAMOSAlib(CST, RDB, OPT,
                                      LUT)  #### initializing the SAMOSA library sampy, it's a mandatory step

        n = np.shape(wfm)[0]

        GEO = type('', (), {})()

        CONF = type('', (), {})()

        CONF.flag_slope = False                    ### flag True commands to include in the model the slope of orbit and surface (this effect usually is included in LookAngles Array)
        CONF.beamsamp_factor = 1                   ### 1 means only one beam per resolution cell is generated in the DDM, the other ones are decimated
        CONF.wf_weighted = True                   ### flag True if the waveform under iteration is weighted
        CONF.N_Look_min = -90                      ### number of the first Look to generate in the DDM (only used if LookAngles array is not passed in input: i.e. set to  None)
        CONF.N_Look_max = 90                       ### number of the last Look to generate in the DDM (only used if LookAngles array is not passed in input: i.e. set to  None)
        CONF.guess_swh = 2                         ### first-guess SWH in meter
        CONF.guess_pu = 1                          ### first-guess Pu
        CONF.guess_nu = 2                          ### first-guess nu (only used in second step of SAMOSA+)
        CONF.lb_epoch = None                       ### lower bound on epoch in sec. If set to None, lower bound will be set to the first time in input array tau
        CONF.lb_swh = -0.5                         ### lower bound on SWH in m
        CONF.lb_pu = 0.2                           ### lower bound on Pu
        CONF.lb_nu = 0                             ### lower bound on nu (only used in second step of SAMOSA+)
        CONF.ub_epoch = None                       ### upper bound on epoch in sec. If set to None, upper bound will be set to the last time in input array tau
        CONF.ub_swh = 30                           ### upper bound on SWH in m
        CONF.ub_pu = 1.5                           ### upper bound on Pu
        CONF.ub_nu = 1e9                           ### upper bound on nu (only used in second step of SAMOSA+)
        CONF.rtk_type = 'samosa+'                  ### choose between 'samosa' or 'samosa+'
        CONF.wght_factor= 1.4705                   ### widening factor of PTR main lobe after Weighting Window Application
        CONF.Lz = CST.c0 / (2. * RDB.Bs)

        # dummy array for interpolation of actual retracker range window
        x = np.arange(wfm.shape[1])

        # Make an altitude rate. Currently zero as L1b rate is nan. Not used anyway with flag_slope false.
        hrate = np.zeros_like(self._l1b.time_orbit.altitude_rate)

        vel = np.sqrt(self.get_l1b_parameter("classifier", "satellite_velocity_x")**2
                      + self.get_l1b_parameter("classifier", "satellite_velocity_y")**2
                      + self.get_l1b_parameter("classifier", "satellite_velocity_z")**2)

        # Loop over waveform indices marked as surface type
        for index in indices:

            LookAngles = 90 - np.linspace(np.rad2deg(self._l1b.classifier.look_angle_start[index]),
                                          np.rad2deg(self._l1b.classifier.look_angle_stop[index]),
                                          num=int(self._l1b.classifier.stack_beams[index]), endpoint=True)
            MaskRanges = None
            GEO.LAT=self._l1b.time_orbit.latitude[index]                              ### latitude in degree for the waveform under iteration
            GEO.LON=self._l1b.time_orbit.longitude[index]                              ### longitude in degree between -180, 180 for the waveform under iteration
            GEO.Height=self._l1b.time_orbit.altitude[index]                            ### Orbit Height in meter for the waveform under iteration
            GEO.Vs=np.squeeze(vel)[index]                       ### Satellite Velocity in m/s for the waveform under iteration
            GEO.Hrate=hrate[index]                   ### Orbit Height rate in m/s for the waveform under iteration
            GEO.Pitch=np.radians(self._l1b.time_orbit.antenna_pitch[index])     ### Altimeter Reference Frame Pitch in radiant
            GEO.Roll=np.radians(self._l1b.time_orbit.antenna_roll[index])       ### Altimeter Reference Frame Roll in radiant
            GEO.nu=0                                                         ### Inverse of the mean square slope
            GEO.track_sign=0                                                 ### if Track Ascending => -1, if Track Descending => +1, set it to zero if flag_slope=False in CONF
            GEO.ThN=np.squeeze(ThNEcho)[index]                                   ### Thermal Noise

            wf = np.array(wfm[index, :]).astype("float64")

            CONF.guess_epoch = epoch0[index]

            # Do retrack in units of Watts as we don't have the scaling factors available for calc_sigma0
            epoch_sec,swh,Pu,misfit,oceanlike_flag=samlib.Retrack_Samosa(tau,wf,LookAngles,MaskRanges,GEO,CONF)

            # SAMOSA returns a dR based upon the retracker chosen bin sampled from tau
            self._range[index] = raw_range[index] + epoch_sec * CST.c0 * 0.5
            sigma0,pval,cval,rval,kval = calc_sigma0(None, Pu, CST, RDB, GEO, epoch_sec, window_del_20_hr_ku_deuso[index],
                                 GEO.LAT, GEO.Height, GEO.Vs,
                                 self._l1b.classifier.transmit_power[index])
            wind_speed = func_wind_speed([sigma0])

            self._power[index] = sigma0

            # Store additional retracker parameters
            self.swh[index] = swh
            self.misfit[index] = misfit
            self.wind_speed[index] = wind_speed
            self.oceanlike_flag[index] = oceanlike_flag
            if True:
                self.epoch[index] = epoch_sec
                self.guess[index] = epoch0[index]
                self.Pu[index] = 65535.0 * Pu/np.max(wf)
                self.pval[index] = pval
                self.cval[index] = cval
                self.rval[index] = rval
                self.kval[index] = kval

        self.register_auxdata_output("samswh", "samosa_swh", self.swh)
        self.register_auxdata_output("samwsp", "samosa_wind_speed", self.wind_speed)

        if "range_bias" in self._options:
            for radar_mode_index in np.arange(3):
                indices = np.where(radar_mode == radar_mode_index)[0]
                if len(indices) == 0:
                    continue
                range_bias = self._options.range_bias[radar_mode_index]
                self._range[indices] -= range_bias

        if "uncertainty" in self._options:
            if self._options.uncertainty.type == "fixed":
                self._uncertainty[:] = self._options.uncertainty.value

        if False:
            import xarray as xr
            outds = xr.Dataset({'SSHunc': (['time_20_ku'], self._l1b.time_orbit.altitude-self._range),
                                'raw_elev' : (['time_20_ku'], Raw_Elevation),
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
                                'wf': (['time_20_ku','bins'], wf_norm),
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


    def _filter_results(self):
        pass
        """ These threshold are based on the SICCI code"""
        #thrs = self._options.filter

        #valid = ANDCondition()
        #valid.add(self.leading_edge_width < thrs.maximum_leading_edge_width)

        # Error flag is also computed for other surface types, do not
        # overide those
        #error_flag = self._flag
        #error_flag[self.indices] = np.logical_not(valid.flag[self.indices])
        #self._flag = error_flag



