
import os
import sys
import pandas
import xarray
import numpy as np
from scipy import interpolate


from pysiral import __version__ as pysiral_version
from pysiral.clocks import StopWatch, UTCTAIConverter
from pysiral.cryosat2 import cs2_procstage2timeliness
from pysiral.errorhandler import ErrorStatus
from pysiral.helper import parse_datetime_str
from pysiral.l1bdata import Level1bData
from pysiral.logging import DefaultLoggingClass
from pysiral.path import filename_from_path


class ESAPDSBaselineD(DefaultLoggingClass):

    def __init__(self, cfg, raise_on_error=False):

        cls_name = self.__class__.__name__
        super(ESAPDSBaselineD, self).__init__(cls_name)
        self.error = ErrorStatus(caller_id=cls_name)

        # Store arguments
        self.raise_on_error = raise_on_error
        self.cfg = cfg

        # Init main class variables
        self.nc = None

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
        self.log.info("Created L1 object in %.3f seconds" % timer.get_seconds())

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
        f = interpolate.interp1d(time_1Hz, variable_1Hz, bounds_error=False, **kwargs)
        variable_20Hz = f(time_20Hz)
        return variable_20Hz

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
        info.set_attribute("pysiral_version", "pysiral_version")

        # General product metadata
        info.set_attribute("mission", "cryosat2")
        info.set_attribute("mission_sensor", "siral")
        info.set_attribute("mission_data_version", "D")
        info.set_attribute("orbit", metadata["abs_orbit_start"])
        info.set_attribute("cycle", metadata["cycle_number"])
        info.set_attribute("mission_data_source", filename_from_path(self.filepath))
        info.set_attribute("timeliness", cs2_procstage2timeliness(metadata["processing_stage"]))

        # Time-Orbit Metadata
        lats = [float(metadata["first_record_lat"])*1e-6, float(metadata["last_record_lat"])*1e-6]
        lons = [float(metadata["first_record_lon"])*1e-6, float(metadata["last_record_lon"])*1e-6]
        info.set_attribute("start_time", parse_datetime_str(metadata["first_record_time"][4:]))   # TAI=....
        info.set_attribute("stop_time", parse_datetime_str(metadata["last_record_time"][4:]))     # TAI=....
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
        # NOTE: Unfortunately is seem to be necessary to first convert the np.datetime64 from xarray to
        #       Pandas TimeStamp and then to datetime.datetime as there seems to be no reliable way
        #       to convert directly from np.datetime64 to datetime.datetime
        tai_datetime = pandas.to_datetime(self.nc.time_20_ku.values).to_pydatetime()

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

        # Set the waveform
        radar_mode = str(self.nc.attrs["sir_op_mode"].strip().lower())
        self.l1.waveform.set_waveform_data(wfm_power, wfm_range, radar_mode)

        # Get the valid flags
        measurement_confident_flag = self.nc.flag_mcd_20_ku.values
        valid_flag = measurement_confident_flag == 0
        self.l1.waveform.set_valid_flag(valid_flag)

    def _set_range_correction_group(self):
        pass

    def _set_surface_type_group(self):
        pass

    def _set_classifier_group(self):
        pass

    @property
    def empty(self):
        return None