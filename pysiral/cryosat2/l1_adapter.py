
import re
from pathlib import Path

import numpy as np
import xarray
from astropy.time import Time
from cftime import num2pydate
from loguru import logger
from scipy import interpolate

from pysiral import __version__ as pysiral_version
from pysiral.classifier import CS2OCOGParameter, CS2PulsePeakiness
from core.clocks import StopWatch
from pysiral.core.flags import ESA_SURFACE_TYPE_DICT
from pysiral.cryosat2 import cs2_procstage2timeliness
from pysiral.helper import parse_datetime_str
from pysiral.l1bdata import Level1bData
from pysiral.l1preproc import Level1PInputHandlerBase


class ESACryoSat2PDSBaselineD(Level1PInputHandlerBase):

    def __init__(self, cfg, raise_on_error=False):

        cls_name = self.__class__.__name__
        super(ESACryoSat2PDSBaselineD, self).__init__(cfg, raise_on_error, cls_name)

        # Init main class variables
        self.nc = None
        self.filepath = None
        self.l1 = None

    @staticmethod
    def translate_opmode2radar_mode(op_mode):
        """ Converts the ESA operation mode str in the pysiral compliant version """
        translate_dict = {"sar": "sar", "lrm": "lrm", "sarin": "sin"}
        return translate_dict.get(op_mode, None)

    def get_l1(self, filepath, polar_ocean_check=None):
        """
        Main entry point to the CryoSat-2 Baseline-D Input Adapter
        :param filepath:
        :param polar_ocean_check:
        :return:
        """

        timer = StopWatch()
        timer.start()

        # Save filepath
        self.filepath = filepath

        # Create an empty Level-1 data object
        self.l1 = Level1bData()

        # Input Validation
        if not Path(filepath).is_file():
            msg = f"Not a valid file: {filepath}"
            logger.warning(msg)
            self.error.add_error("invalid-filepath", msg)
            return self.empty

        # Parse the input file
        self._read_input_netcdf(filepath, attributes_only=True)
        if self.nc is None:
            return self.empty

        # CAVEAT: An issue has been identified with baseline-D L1b data when the orbit solution
        # is based on predicted orbits and not the DORIS solution (Nov 2020).
        # The source of the orbit data can be identified by the `vector_source` global attribute
        # in the L1b source files. This can take/should take the following values:
        #
        #     nrt:  "fos predicted" (predicted orbit)
        #           "doris_navigator" (DORIS Nav solution)
        #
        #     rep:  "doris_precise" (final and precise DORIS solution)
        #
        # To prevent l1 data with erroneous orbit solution entering the processing chain, l1 data
        # with the predicted orbit can be excluded here. The process of exclusion requires to set
        # a flag in the l1 processor definition for the input handler:
        #
        #   exclude_predicted_orbits: True
        #
        exclude_predicted_orbits = self.cfg.get("exclude_predicted_orbits", False)
        is_predicted_orbit = self.nc.vector_source.lower().strip() == "fos predicted"
        if is_predicted_orbit and exclude_predicted_orbits:
            logger.warning("Predicted orbit solution detected -> skip file")
            return self.empty

        # Get metadata
        self._set_input_file_metadata()
        if polar_ocean_check is not None:
            has_polar_ocean_data = polar_ocean_check.has_polar_ocean_segments(self.l1.info)
            if not has_polar_ocean_data:
                timer.stop()
                return self.empty

        # Polar ocean check passed, now fill the rest of the l1 data groups
        self._set_l1_data_groups()

        timer.stop()
        logger.info("- Created L1 object in %.3f seconds" % timer.get_seconds())

        # Return the l1 object
        return self.l1

    @staticmethod
    def get_wfm_range(window_delay, n_range_bins):
        """
        Returns the range for each waveform bin based on the window delay and the number of range bins
        :param window_delay: The two-way delay to the center of the range window in seconds
        :param n_range_bins: The number of range bins (256: sar, 512: sin)
        :return: The range for each waveform bin as array (time, ns)
        """
        lightspeed = 299792458.0
        bandwidth = 320000000.0
        # The two-way delay time give the distance to the central bin
        central_window_range = window_delay * lightspeed / 2.0
        # Calculate the offset from the center to the first range bin
        window_size = (n_range_bins * lightspeed) / (4.0 * bandwidth)
        first_bin_offset = window_size / 2.0
        # Calculate the range increment for each bin
        range_increment = np.arange(n_range_bins) * lightspeed / (4.0 * bandwidth)

        # Reshape the arrays
        range_offset = np.tile(range_increment, (window_delay.shape[0], 1)) - first_bin_offset
        window_range = np.tile(central_window_range, (n_range_bins, 1)).transpose()

        return window_range + range_offset

    @staticmethod
    def interp_1hz_to_20hz(variable_1hz, time_1hz, time_20hz, **kwargs):
        """
        Computes a simple linear interpolation to transform a 1Hz into a 20Hz variable
        :param variable_1hz: an 1Hz variable array
        :param time_1hz: 1Hz reference time
        :param time_20hz: 20 Hz reference time
        :return: the interpolated 20Hz variable
        """
        error_status = False
        try:
            f = interpolate.interp1d(time_1hz, variable_1hz, bounds_error=False, **kwargs)
            variable_20hz = f(time_20hz)
        except ValueError:
            fill_value = np.nan
            variable_20hz = np.full(time_20hz.shape, fill_value)
            error_status = True
        return variable_20hz, error_status

    def _read_input_netcdf(self, filepath, **kwargs):
        """ Read the netCDF file via xarray """
        try:
            self.nc = xarray.open_dataset(filepath, decode_times=False, mask_and_scale=True)
        except:
            msg = "Error encountered by xarray parsing: %s" % filepath
            self.error.add_error("xarray-parse-error", msg)
            self.nc = None
            logger.warning(msg)
            return

    def _set_input_file_metadata(self):
        """ Fill the product info """

        # Short cuts
        metadata = self.nc.attrs
        info = self.l1.info

        # Processing environment metadata
        info.set_attribute("pysiral_version", pysiral_version)

        # General product metadata
        info.set_attribute("mission", "cryosat2")
        info.set_attribute("mission_sensor", "siral")
        info.set_attribute("mission_data_version", "D")
        info.set_attribute("orbit", metadata["abs_orbit_start"])
        info.set_attribute("rel_orbit", metadata["rel_orbit_number"])
        info.set_attribute("cycle", metadata["cycle_number"])
        info.set_attribute("mission_data_source", Path(self.filepath).name)
        info.set_attribute("timeliness", cs2_procstage2timeliness(metadata["processing_stage"]))

        # Time-Orbit Metadata
        tcs_tai = parse_datetime_str(metadata["first_record_time"][4:])
        tce_tai = parse_datetime_str(metadata["last_record_time"][4:])
        tcs_utc, tce_utc = Time([tcs_tai, tce_tai], scale="tai").utc.datetime

        lats = [float(metadata["first_record_lat"])*1e-6, float(metadata["last_record_lat"])*1e-6]
        lons = [float(metadata["first_record_lon"])*1e-6, float(metadata["last_record_lon"])*1e-6]
        info.set_attribute("start_time", tcs_utc)
        info.set_attribute("stop_time", tce_utc)
        info.set_attribute("lat_min", np.amin(lats))
        info.set_attribute("lat_max", np.amax(lats))
        info.set_attribute("lon_min", np.amin(lons))
        info.set_attribute("lon_max", np.amax(lons))

        # Product Content Metadata
        for mode in ["sar", "sin", "lrm"]:
            percent_value = 0.0
            if metadata["sir_op_mode"].strip().lower() == mode:
                percent_value = 100.
            info.set_attribute(f"{mode}_mode_percent", percent_value)
        info.set_attribute("open_ocean_percent", float(metadata["open_ocean_percent"])*0.01)

    def _set_l1_data_groups(self):
        """
        Fill all data groups of the Level-1 data object with the content of the netCDF file. This is just the
        overview method, see specific sub-methods below
        :return: None
        """
        self._set_time_orbit_data_group()
        self._set_waveform_data_group()
        self._set_range_correction_group()
        self._set_surface_type_group()
        self._set_classifier_group()

    def _set_time_orbit_data_group(self):
        """
        Transfer the time orbit parameter from the netcdf to l1 data object
        :return: None
        """

        # Transfer the timestamp
        # NOTE: Here it is critical that the xarray does not automatically decodes time since it is
        #       difficult to work with the numpy datetime64 date format. Better to compute datetimes using
        #       a know num2pydate conversion
        tai_datetime = num2pydate(self.nc.time_20_ku.values, units=self.nc.time_20_ku.units)
        self.l1.time_orbit.timestamp = Time(tai_datetime, scale="tai").utc.datetime

        # Set the geolocation
        self.l1.time_orbit.set_position(
            self.nc.lon_20_ku.values,
            self.nc.lat_20_ku.values,
            self.nc.alt_20_ku.values,
            self.nc.orb_alt_rate_20_ku.values)

        # Set antenna attitude
        self.l1.time_orbit.set_antenna_attitude(
            self.nc.off_nadir_pitch_angle_str_20_ku.values,
            self.nc.off_nadir_roll_angle_str_20_ku.values,
            self.nc.off_nadir_yaw_angle_str_20_ku.values)

    def _set_waveform_data_group(self):
        """
        Transfer of the waveform group to the Level-1 object. This includes
          1. the computation of waveform power in Watts
          2. the computation of the window delay in meter for each waveform bin
          3. extraction of the waveform valid flag
        :return: None
        """

        # Get the waveform
        # NOTE: Convert the waveform units to Watts. From the documentation:is applied as follows:
        #       pwr_waveform_20_ku(time, ns) * echo_scale_factor_20_ku(time, ns) * 2 ^ echo_scale_pwr_20_ku(time)
        wfm_linear = self.nc.pwr_waveform_20_ku.values

        # Get the shape of the waveform array
        dim_time, dim_ns = wfm_linear.shape

        # Scaling parameter are 1D -> Replicate to same shape as waveform array
        echo_scale_factor = self.nc.echo_scale_factor_20_ku.values
        echo_scale_pwr = self.nc.echo_scale_pwr_20_ku.values
        echo_scale_factor = np.tile(echo_scale_factor, (dim_ns, 1)).transpose()
        echo_scale_pwr = np.tile(echo_scale_pwr, (dim_ns, 1)).transpose()

        # Convert the waveform from linear counts to Watts
        wfm_power = wfm_linear*echo_scale_factor * 2.0**echo_scale_pwr

        # Get the window delay
        # From the documentation:
        #   Calibrated 2-way window delay: distance from CoM to middle range window (at sample ns/2 from 0).
        #   It includes all the range corrections given in the variable instr_cor_range and in the
        #   variable uso_cor_20_ku. This is a 2-way time and 2-way corrections are applied.
        window_delay = self.nc.window_del_20_ku.values

        # Convert window delay to range for each waveform range bin
        wfm_range = self.get_wfm_range(window_delay, dim_ns)

        # Set the waveform
        op_mode = str(self.nc.attrs["sir_op_mode"].strip().lower())
        radar_mode = self.translate_opmode2radar_mode(op_mode)
        self.l1.waveform.set_waveform_data(wfm_power, wfm_range, radar_mode)

        # --- Get the valid flag ---
        #
        # From the documentation
        # :comment = "Measurement confidence flags. Generally the MCD flags indicate problems when set.
        #             If the whole MCD is 0 then no problems or non-nominal conditions were detected.
        #             Serious errors are indicated by setting the most significant bit, i.e. block_degraded,
        #             in which case the block must not be processed. Other error settings can be regarded
        #             as warnings.";
        #
        # :flag_masks = -2147483648, block_degraded        <- most severe error
        #                1073741824, blank_block
        #                536870912, datation_degraded
        #                268435456, orbit_prop_error
        #                134217728, orbit_file_change
        #                67108864, orbit_gap
        #                33554432, echo_saturated
        #                16777216, other_echo_error
        #                8388608, sarin_rx1_error
        #                4194304, sarin_rx2_error
        #                2097152, window_delay_error
        #                1048576, agc_error
        #                524288, cal1_missing
        #                262144, cal1_default
        #                131072, doris_uso_missing
        #                65536, ccal1_default
        #                32768, trk_echo_error
        #                16384, echo_rx1_error
        #                8192, echo_rx2_error
        #                4096, npm_error                   <- Defined as maximum permissible error level
        #                2048, cal1_pwr_corr_type
        #                128, phase_pert_cor_missing       <- Seems to be always set for SARin
        #                64, cal2_missing
        #                32, cal2_default
        #                16, power_scale_error
        #                8, attitude_cor_missing
        #                1, phase_pert_cor_default
        measurement_confident_flag = self.nc.flag_mcd_20_ku.values
        valid_flag = (measurement_confident_flag >= 0) & (measurement_confident_flag <= 4096)
        self.l1.waveform.set_valid_flag(valid_flag)

    def _set_range_correction_group(self):
        """
        Transfer the range corrections defined in the l1p config file to the Level-1 object
        NOTE: The range corrections are all in 1 Hz and must be interpolated to 20Hz
        :return: None
        """

        # Get the reference times for interpolating the range corrections from 1Hz -> 20Hz
        time_1hz = self.nc.time_cor_01.values
        time_20hz = self.nc.time_20_ku.values

        # Loop over all range correction variables defined in the processor definition file
        for key in self.cfg.range_correction_targets.keys():
            pds_var_name = self.cfg.range_correction_targets[key]
            variable_1hz = getattr(self.nc, pds_var_name)
            variable_20hz, error_status = self.interp_1hz_to_20hz(variable_1hz.values, time_1hz, time_20hz)
            if error_status:
                msg = f"- Error in 20Hz interpolation for variable `{pds_var_name}` -> set only dummy"
                logger.warning(msg)
            self.l1.correction.set_parameter(key, variable_20hz)

    def _set_surface_type_group(self):
        """
        Transfer of the surface type flag to the Level-1 object
        NOTE: In the current state (TEST dataset), the surface type flag is only 1 Hz. A nearest neighbour
              interpolation is used to get the 20Hz surface type flag.
        :return: None
        """

        # Get the reference times for interpolating the flag from 1Hz -> 20Hz
        time_1hz = self.nc.time_cor_01.values
        time_20hz = self.nc.time_20_ku.values

        # Interpolate 1Hz surface type flag to 20 Hz
        surface_type_1hz = self.nc.surf_type_01.values
        surface_type_20hz, error_status = self.interp_1hz_to_20hz(surface_type_1hz, time_1hz, time_20hz, kind="nearest")
        if error_status:
            msg = "- Error in 20Hz interpolation for variable `surf_type_01` -> set only dummy"
            logger.warning(msg)

        # Set the flag
        for key in ESA_SURFACE_TYPE_DICT.keys():
            flag = surface_type_20hz == ESA_SURFACE_TYPE_DICT[key]
            self.l1.surface_type.add_flag(flag, key)

    def _set_classifier_group(self):
        """
        Transfer the classifiers defined in the l1p config file to the Level-1 object.
        NOTE: It is assumed that all classifiers are 20Hz

        In addition, a few legacy parameter are computed based on the waveform counts that is only available at
        this stage. Computation of other parameter such as sigma_0, leading_edge_width, ... are moved to the
        post-processing

        :return: None
        """
        # Loop over all classifier variables defined in the processor definition file
        for key in self.cfg.classifier_targets.keys():
            variable_20hz = getattr(self.nc, self.cfg.classifier_targets[key])
            self.l1.classifier.add(variable_20hz, key)

        # Calculate the OCOG Parameter (CryoSat-2 notation)
        ocog = CS2OCOGParameter(self.l1.waveform.power)
        self.l1.classifier.add(ocog.width, "ocog_width")
        self.l1.classifier.add(ocog.amplitude, "ocog_amplitude")

        # Get satellite velocity vector (classifier needs to be vector -> manual extraction needed)
        satellite_velocity_vector = self.nc.sat_vel_vec_20_ku.values
        self.l1.classifier.add(satellite_velocity_vector[:, 0], "satellite_velocity_x")
        self.l1.classifier.add(satellite_velocity_vector[:, 1], "satellite_velocity_y")
        self.l1.classifier.add(satellite_velocity_vector[:, 2], "satellite_velocity_z")

    @property
    def empty(self):
        return None


class ESACryoSat2PDSBaselineDPatchFES(ESACryoSat2PDSBaselineD):
    def __init__(self, cfg, raise_on_error=False):
        ESACryoSat2PDSBaselineD.__init__(self, cfg, raise_on_error)

    def _set_l1_data_groups(self):
        ESACryoSat2PDSBaselineD._set_l1_data_groups(self)
        fespath = self._get_fes_path(self.filepath)
        if not Path(fespath).is_file():
            msg = f"Not a valid file: {fespath}"
            logger.warning(msg)
            self.error.add_error("invalid-filepath", msg)
            raise FileNotFoundError
        try:
            nc_fes = xarray.open_dataset(fespath, decode_times=False, mask_and_scale=True)

            # time_1hz = self.nc.time_cor_01.values
            # time_20hz = self.nc.time_20_ku.values

            msg = f"Patching FES2014b tide data from: {fespath}"
            logger.info(msg)

            # ocean_tide_elastic: ocean_tide_01
            variable_20hz = getattr(nc_fes, 'ocean_tide_20')
            # variable_20hz, error_status = self.interp_1hz_to_20hz(variable_1hz.values, time_1hz, time_20hz)
            # if error_status:
            #    msg = "- Error in 20Hz interpolation for variable `%s` -> set only dummy" % 'ocean_tide_01'
            #    logger.warning(msg)
            #    raise FileNotFoundError
            self.l1.correction.set_parameter('ocean_tide_elastic', variable_20hz)

            # ocean_tide_long_period: ocean_tide_eq_01
            variable_20hz = getattr(nc_fes, 'ocean_tide_eq_20')
            # variable_20hz, error_status = self.interp_1hz_to_20hz(variable_1hz.values, time_1hz, time_20hz)
            # if error_status:
            #    msg = "- Error in 20Hz interpolation for variable `%s` -> set only dummy" % 'ocean_tide_eq_01'
            #    logger.warning(msg)
            #    raise FileNotFoundError
            self.l1.correction.set_parameter('ocean_tide_long_period', variable_20hz)

            # ocean_loading_tide: load_tide_01
            variable_20hz = getattr(nc_fes, 'load_tide_20')
            # variable_20hz, error_status = self.interp_1hz_to_20hz(variable_1hz.values, time_1hz, time_20hz)
            # if error_status:
            #     msg = "- Error in 20Hz interpolation for variable `%s` -> set only dummy" % 'load_tide_01'
            #     logger.warning(msg)
            #     raise FileNotFoundError
            self.l1.correction.set_parameter('ocean_loading_tide', variable_20hz)
        except:
            msg = f"Error encountered by xarray parsing: {fespath}"
            self.error.add_error("xarray-parse-error", msg)
            self.nc = None
            logger.warning(msg)
            raise FileNotFoundError

    def _get_fes_path(self, filepath):
        # TODO: get the substitutions to make from config file. Get a list of pairs of sub 'this' to 'that'.
        # pathsubs = [ ( 'L1B', 'L1B/FES2014' ), ( 'nc', 'fes2014b.nc' ) ]
        newpath = str(filepath)
        p = re.compile('L1B')
        newpath = p.sub('L1B/FES2014', newpath)
        p = re.compile('nc')
        newpath = p.sub('fes2014b.nc', newpath)
        p = re.compile('TEST')
        newpath = p.sub('LTA_', newpath)
        return newpath


class ESACryoSat2PDSBaselineDPatchFESArctide(ESACryoSat2PDSBaselineDPatchFES):
    def __init__(self, cfg, raise_on_error=False):
        ESACryoSat2PDSBaselineDPatchFES.__init__(self, cfg, raise_on_error)

    def _set_l1_data_groups(self):
        ESACryoSat2PDSBaselineDPatchFES._set_l1_data_groups(self)
        arcpath = self._get_arctide_path(self.filepath)
        if not Path(arcpath).is_file():
            msg = f"Not a valid file: {arcpath}"
            logger.warning(msg)
            self.error.add_error("invalid-filepath", msg)
            # The handling of missing files here is different so that we can still process
            # south files even though we don't have Arctide for them
            self.l1.correction.set_parameter('ocean_tide_elastic_2',
                                             self.l1.correction.get_parameter_by_name('ocean_tide_elastic'))
        else:
            nc_arc = xarray.open_dataset(arcpath, decode_times=False, mask_and_scale=True)

            # time_1hz = self.nc.time_cor_01.values
            # time_20hz = self.nc.time_20_ku.values

            msg = f"Patching ARCTIDE tide data from: {arcpath}"
            logger.info(msg)

            # ocean_tide_elastic: ocean_tide_01
            variable_20hz = getattr(nc_arc, 'tide_Arctic')
            # variable_20hz, error_status = self.interp_1hz_to_20hz(variable_1hz.values, time_1hz, time_20hz)
            # if error_status:
            #    msg = "- Error in 20Hz interpolation for variable `%s` -> set only dummy" % 'ocean_tide_01'
            #    logger.warning(msg)
            #    raise FileNotFoundError
            nans_indices = np.where(np.isnan(variable_20hz))[0]
            if len(nans_indices) > 0:
                msg = 'Arctide file had {numnan} NaN values of {numval}. These have been replaced with FES2014b data'.format(numnan=len(nans_indices), numval=len(variable_20hz))
                logger.warning(msg)
                variable_20hz[nans_indices] = self.l1.correction.get_parameter_by_name('ocean_tide_elastic')[nans_indices]
            self.l1.correction.set_parameter('ocean_tide_elastic_2', self.l1.correction.get_parameter_by_name('ocean_tide_elastic'))
            self.l1.correction.set_parameter('ocean_tide_elastic', variable_20hz)


    def _get_arctide_path(self, filepath):
        # TODO: get the substitutions to make from config file. Get a list of pairs of sub 'this' to 'that'.
        # pathsubs = [ ( 'L1B', 'L1B/FES2014' ), ( 'nc', 'fes2014b.nc' ) ]
        newpath = str(filepath)
        p = re.compile('L1B')
        newpath = p.sub('L1B/ARCTIDE', newpath)
        p = re.compile('nc')
        newpath = p.sub('RegAT_Arctic_tides_v1.2.nc', newpath)
        p = re.compile('TEST')
        newpath = p.sub('LTA_', newpath)
        return newpath