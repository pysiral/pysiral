
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

        # Get metadata
        self._get_input_file_metadata()

        if polar_ocean_check is not None:
            has_polar_ocean_data = polar_ocean_check.has_polar_ocean_segments(self.l1.info)
            if not has_polar_ocean_data:
                return None

    def _read_input_netcdf(self, filepath, attributes_only=False):
        timer = StopWatch()
        timer.start()
        self.nc = xarray.open_dataset(filepath)
        timer.stop()
        self.log.info("Read netCDF file in %s" % timer.get_duration())


    def _get_input_file_metadata(self):
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

    @property
    def empty(self):
        return None