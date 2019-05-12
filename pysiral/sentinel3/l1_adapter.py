
import os
import sys

import xarray
import numpy as np
from scipy import interpolate
from netCDF4 import num2date

from pysiral import __version__ as pysiral_version
from pysiral.classifier import CS2OCOGParameter, CS2LTPP, CS2PulsePeakiness
from pysiral.clocks import StopWatch, UTCTAIConverter
from pysiral.cryosat2 import cs2_procstage2timeliness
from pysiral.errorhandler import ErrorStatus
from pysiral.helper import parse_datetime_str
from pysiral.l1bdata import Level1bData
from pysiral.logging import DefaultLoggingClass
from pysiral.path import filename_from_path
from pysiral.surface_type import ESA_SURFACE_TYPE_DICT


class Sentinel3CODAL2Wat(DefaultLoggingClass):

    def __init__(self, cfg, raise_on_error=False):

        cls_name = self.__class__.__name__
        super(Sentinel3CODAL2Wat, self).__init__(cls_name)
        self.error = ErrorStatus(caller_id=cls_name)

        # Store arguments
        self.raise_on_error = raise_on_error
        self.cfg = cfg

        # Init main class variables
        self.nc = None

    @staticmethod
    def translate_opmode2radar_mode(op_mode):
        """ Converts the ESA operation mode str in the pysiral compliant version """
        translate_dict = {"sar": "sar", "lrm": "lrm", "sarin": "sin"}
        return translate_dict.get(op_mode, None)

    def get_l1(self, filepath, polar_ocean_check=None):
        """
        Main entry point to the CryoSat-2 Baseline-D Input Adapter
        :param filepath:
        :return:
        """

        timer = StopWatch()
        timer.start()

        # Save filepath
        self.filepath = filepath

        # Create an empty Level-1 data object
        self.l1 = Level1bData()

        # Input Validation
        if not os.path.isfile(filepath):
            msg = "Not a valid file: %s" % filepath
            self.log.warning(msg)
            self.error.add_error("invalid-filepath", msg)
            return self.empty

        # Parse the input file
        self._read_input_netcdf(filepath, attributes_only=True)

        if self.error.status:
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
        self.log.info("- Created L1 object in %.3f seconds" % timer.get_seconds())

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
        # The two way delay time give the distance to the central bin
        central_window_range = window_delay * lightspeed / 2.0
        # Calculate the offset from the center to the first range bin
        window_size = (n_range_bins * lightspeed) / (4.0 * bandwidth)
        first_bin_offset = window_size / 2.0
        # Calculate the range increment for each bin
        range_increment = np.arange(n_range_bins) * lightspeed / (4.0 * bandwidth)

        # Reshape the arrays
        range_offset = np.tile(range_increment, (window_delay.shape[0], 1)) - first_bin_offset
        window_range = np.tile(central_window_range, (n_range_bins, 1)).transpose()

        # Compute the range for each bin and return
        wfm_range = window_range + range_offset
        return wfm_range

    @staticmethod
    def interp_1Hz_to_20Hz(variable_1Hz, time_1Hz, time_20Hz, **kwargs):
        """
        Computes a simple linear interpolation to transform a 1Hz into a 20Hz variable
        :param variable_1Hz: an 1Hz variable array
        :param time_1Hz: 1Hz reference time
        :param time_20Hz: 20 Hz reference time
        :return: the interpolated 20Hz variable
        """
        error_status = False
        try:
            f = interpolate.interp1d(time_1Hz, variable_1Hz, bounds_error=False, **kwargs)
            variable_20Hz = f(time_20Hz)
        except ValueError:
            fill_value = np.nan
            variable_20Hz = np.full(time_20Hz.shape, fill_value)
            error_status = True
        return variable_20Hz, error_status

    def _read_input_netcdf(self, filepath, attributes_only=False):
        """ Read the netCDF file via xarray """
        try:
            self.nc = xarray.open_dataset(filepath, decode_times=False, mask_and_scale=True)
        except:
            msg = "Error encountered by xarray parsing: %s" % filepath
            self.error.add_error("xarray-parse-error", msg)
            self.log.warning(msg)
            return

    def _set_input_file_metadata(self):
        """ Fill the product info """

        # Short cuts
        metadata =  self.nc.attrs
        info = self.l1.info

        # Processing environment metadata
        info.set_attribute("pysiral_version", pysiral_version)

        # General product metadata
        mission = metadata["mission_name"].lower().replace(" ", "")
        info.set_attribute("mission", mission)
        info.set_attribute("mission_sensor", "sral")
        info.set_attribute("mission_data_version", metadata["source"])
        info.set_attribute("orbit", metadata["absolute_rev_number"])
        info.set_attribute("cycle", metadata["cycle_number"])

        product_name = metadata["product_name"]
        info.set_attribute("mission_data_source", product_name)
        timeliness_tag = product_name.split("_")[-2]
        info.set_attribute("timeliness", self.cfg.timeliness_dict[timeliness_tag])

        # Time-Orbit Metadata
        lats = [float(metadata["first_meas_lat"]), float(metadata["last_meas_lat"])]
        lons = [float(metadata["first_meas_lon"]), float(metadata["last_meas_lon"])]
        info.set_attribute("start_time", parse_datetime_str(metadata["first_meas_time"][4:]))
        info.set_attribute("stop_time", parse_datetime_str(metadata["last_meas_time"][4:]))
        info.set_attribute("lat_min", np.amin(lats))
        info.set_attribute("lat_max", np.amax(lats))
        info.set_attribute("lon_min", np.amin(lons))
        info.set_attribute("lon_max", np.amax(lons))

        # Product Content Metadata
        for mode in ["sar", "sin", "lrm"]:
            percent_value = 0.0
            if metadata["sir_op_mode"].strip().lower() == mode:
                percent_value = 100.
            info.set_attribute("{}_mode_percent".format(mode), percent_value)
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
        #       a know num2date conversion
        tai_datetime = num2date(self.nc.time_20_ku.values, units=self.nc.time_20_ku.units)
        converter = UTCTAIConverter()
        utc_timestamp = converter.tai2utc(tai_datetime, check_all=False)
        self.l1.time_orbit.timestamp = utc_timestamp

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

        # Make sure that parameter are float and not double
        # -> Import for cythonized algorithm parts (ctfrma specifically uses floats)
        wfm_power = wfm_power.astype(np.float32)
        wfm_range = wfm_range.astype(np.float32)

        # Set the waveform
        op_mode = str(self.nc.attrs["sir_op_mode"].strip().lower())
        radar_mode = self.translate_opmode2radar_mode(op_mode)
        self.l1.waveform.set_waveform_data(wfm_power, wfm_range, radar_mode)

        # Get the valid flags
        measurement_confident_flag = self.nc.flag_mcd_20_ku.values
        valid_flag = measurement_confident_flag == 0
        self.l1.waveform.set_valid_flag(valid_flag)

    def _set_range_correction_group(self):
        """
        Transfer the range corrections defined in the l1p config file to the Level-1 object
        NOTE: The range corrections are all in 1 Hz and must be interpolated to 20Hz
        :return: None
        """

        # Get the reference times for interpolating the range corrections from 1Hz -> 20Hz
        time_1Hz =  self.nc.time_cor_01.values
        time_20Hz = self.nc.time_20_ku.values

        # Loop over all range correction variables defined in the processor definition file
        for key in self.cfg.range_correction_targets.keys():
            pds_var_name = self.cfg.range_correction_targets[key]
            variable_1Hz = getattr(self.nc, pds_var_name)
            variable_20Hz, error_status = self.interp_1Hz_to_20Hz(variable_1Hz.values, time_1Hz, time_20Hz)
            if error_status:
                msg = "- Error in 20Hz interpolation for variable `%s` -> set only dummy" % pds_var_name
                self.log.warning(msg)
            self.l1.correction.set_parameter(key, variable_20Hz)

    def _set_surface_type_group(self):
        """
        Transfer of the surface type flag to the Level-1 object
        NOTE: In the current state (TEST dataset), the surface type flag is only 1 Hz. A nearest neighbour
              interpolation is used to get the 20Hz surface type flag.
        :return: None
        """

        # Get the reference times for interpolating the flag from 1Hz -> 20Hz
        time_1Hz =  self.nc.time_cor_01.values
        time_20Hz = self.nc.time_20_ku.values

        # Interpolate 1Hz surface type flag to 20 Hz
        surface_type_1Hz = self.nc.surf_type_01.values
        surface_type_20Hz, error_status = self.interp_1Hz_to_20Hz(surface_type_1Hz, time_1Hz, time_20Hz, kind="nearest")
        if error_status:
            msg = "- Error in 20Hz interpolation for variable `surf_type_01` -> set only dummy"
            self.log.warning(msg)

        # Set the flag
        for key in ESA_SURFACE_TYPE_DICT.keys():
            flag = surface_type_20Hz == ESA_SURFACE_TYPE_DICT[key]
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
            variable_20Hz = getattr(self.nc, self.cfg.classifier_targets[key])
            self.l1.classifier.add(variable_20Hz, key)

        # Calculate Parameters from waveform counts
        # XXX: This is a legacy of the CS2AWI IDL processor
        #      Threshold defined for waveform counts not power in dB
        wfm_counts = self.nc.pwr_waveform_20_ku.values

        # Calculate the OCOG Parameter (CryoSat-2 notation)
        ocog = CS2OCOGParameter(wfm_counts)
        self.l1.classifier.add(ocog.width, "ocog_width")
        self.l1.classifier.add(ocog.amplitude, "ocog_amplitude")

        # Calculate the Peakiness (CryoSat-2 notation)
        pulse = CS2PulsePeakiness(wfm_counts)
        self.l1.classifier.add(pulse.peakiness, "peakiness")
        self.l1.classifier.add(pulse.peakiness_r, "peakiness_r")
        self.l1.classifier.add(pulse.peakiness_l, "peakiness_l")

        # fmi version: Calculate the LTPP
        ltpp = CS2LTPP(wfm_counts)
        self.l1.classifier.add(ltpp.ltpp, "late_tail_to_peak_power")

        # Get satellite velocity vector (classifier needs to be vector -> manual extraction needed)
        satellite_velocity_vector = self.nc.sat_vel_vec_20_ku.values
        self.l1.classifier.add(satellite_velocity_vector[:, 0], "satellite_velocity_x")
        self.l1.classifier.add(satellite_velocity_vector[:, 1], "satellite_velocity_y")
        self.l1.classifier.add(satellite_velocity_vector[:, 2], "satellite_velocity_z")



    @property
    def empty(self):
        return None