
import os
import sys
import xarray
import numpy as np

from pysiral import __version__ as pysiral_version
from pysiral.clocks import StopWatch
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
                return self.empty

        # Polar ocean check passed, now fill the rest of the l1 data groups
        self._set_l1_data_groups()

        # Return the l1 object
        return self.l1

    def _read_input_netcdf(self, filepath, attributes_only=False):
        """ Read the netCDF file via xarray """
        timer = StopWatch()
        timer.start()
        try:
            self.nc = xarray.open_dataset(filepath)
        except:
            timer.stop()
            msg = "Error encountered by xarray parsing: %s" % filepath
            self.error.add_error("xarray-parse-error", msg)
            self.log.warning(msg)
            return

        timer.stop()
        self.log.info("Read netCDF file in %s" % timer.get_duration())

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

        # Set the geolocation
        self.l1.time_orbit.set_position(
            self.nc.lon_20_ku[:],
            self.nc.lat_20_ku[:],
            self.nc.alt_20_ku[:],
            self.nc.orb_alt_rate_20_ku[:])

        # Set antenna attitude
        self.l1.time_orbit.set_antenna_attitude(
            self.nc.off_nadir_pitch_angle_str_20_ku[:],
            self.nc.off_nadir_roll_angle_str_20_ku[:],
            self.nc.off_nadir_yaw_angle_str_20_ku[:])

    def _set_waveform_data_group(self):
        stop
        pass

    def _set_range_correction_group(self):
        pass

    def _set_surface_type_group(self):
        pass

    def _set_classifier_group(self):
        pass

    @property
    def empty(self):
        return None