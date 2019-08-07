
import os
import xmltodict

import xarray
import numpy as np
from scipy import interpolate
from netCDF4 import num2date

from pysiral import __version__ as pysiral_version
from pysiral.classifier import CS2OCOGParameter, CS2LTPP, CS2PulsePeakiness
from pysiral.clocks import StopWatch
from pysiral.errorhandler import ErrorStatus
from pysiral.helper import parse_datetime_str
from pysiral.l1bdata import Level1bData
from pysiral.logging import DefaultLoggingClass
from pysiral.path import folder_from_filename
from pysiral.surface_type import ESA_SURFACE_TYPE_DICT


class Sentinel3CODAL2Wat(DefaultLoggingClass):

    def __init__(self, cfg, raise_on_error=False):
        """
        Input handler for Sentinel-3 L2WAT netCDF files from the CODA.
        :param cfg: A treedict object (root.input_handler.options) from the corresponding Level-1 pre-processor
                    config file
        :param raise_on_error: Boolean value if the class should raise an exception upon an error (default: False)
        """

        cls_name = self.__class__.__name__
        super(Sentinel3CODAL2Wat, self).__init__(cls_name)
        self.error = ErrorStatus(caller_id=cls_name)

        # Store arguments
        self.raise_on_error = raise_on_error
        self.cfg = cfg

        # Init main class variables
        self.nc = None

        # Debug variables
        self.timer = None

    def get_l1(self, filepath, polar_ocean_check=None):
        """
        Create a Level-1 data container from Sentinel-3 CODA L2WAT files
        :param filepath: The full file path to the netCDF file
        :param polar_ocean_check:
        :return: The parsed (or empty) Level-1 data container
        """

        #  for debug purposes
        self.timer = StopWatch()
        self.timer.start()

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

        # Parse xml header file
        self._parse_xml_manifest(filepath)

        # Parse the input netCDF file
        self._read_input_netcdf(filepath)
        if self.error.status:
            return self.empty

        # Get metadata
        self._set_input_file_metadata()

        # Test if input file contains data over polar oceans (optional)
        if polar_ocean_check is not None:
            has_polar_ocean_data = polar_ocean_check.has_polar_ocean_segments(self.l1.info)
            if not has_polar_ocean_data:
                self.timer.stop()
                return self.empty

        # Polar ocean check passed, now fill the rest of the l1 data groups
        self._set_l1_data_groups()

        self.timer.stop()
        self.log.info("- Created L1 object in %.3f seconds" % self.timer.get_seconds())

        # Return the l1 object
        return self.l1

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

    @staticmethod
    def parse_sentinel3_l1b_xml_header(filename):
        """
        Reads the XML header file of a Sentinel 3 L1b Data set
        and returns the contents as an OrderedDict
        """
        with open(filename) as fd:
            content_odereddict = xmltodict.parse(fd.read())
        return content_odereddict[u'xfdu:XFDU']

    def _parse_xml_manifest(self, filepath):
        """
        Parse the Sentinel-3 XML header file and extract key attributes for filtering
        :param filepath: the filepath for the netcdf
        :return: None
        """
        # Retrieve header information from mission settings
        xml_header_file = self.cfg.xml_manifest
        dataset_folder = folder_from_filename(filepath)
        filename_header = os.path.join(dataset_folder, xml_header_file)
        self._xmlh = self.parse_sentinel3_l1b_xml_header(filename_header)

    def _get_xml_content(self, section_name, tag):
        """ Returns the generalProductInformation content of the xml manifest
        :return: dictionary
        """

        # Extract Metadata
        metadata = self._xmlh["metadataSection"]["metadataObject"]

        # Extract General Product Info
        index = self.cfg.xml_metadata_object_index[section_name]
        product_info = metadata[index]["metadataWrap"]["xmlData"]
        product_info = product_info[tag]

        return product_info

    def _read_input_netcdf(self, filepath):
        """
        Read the netCDF file via xarray
        :param filepath: The full filepath to the netCDF file
        :return: none
        """
        try:
            self.nc = xarray.open_dataset(filepath, decode_times=False, mask_and_scale=True)
        except:
            msg = "Error encountered by xarray parsing: %s" % filepath
            self.error.add_error("xarray-parse-error", msg)
            self.log.warning(msg)
            return

    def _set_input_file_metadata(self):
        """
        Populates the product info segment of the Level1Data object with information from
        the global attributes of the netCDF and content of the xml manifest
        :return: None
        """

        # Short cuts
        metadata =  self.nc.attrs
        info = self.l1.info

        # Get xml manifest content
        product_info = self._get_xml_content("generalProductInformation", "sentinel3:generalProductInformation")
        sral_info = self._get_xml_content("sralProductInformation", "sralProductInformation")

        # Processing environment metadata
        info.set_attribute("pysiral_version", pysiral_version)

        # General product metadata
        mission = metadata["mission_name"].lower().replace(" ", "")
        info.set_attribute("mission", str(mission))
        info.set_attribute("mission_sensor", "sral")
        info.set_attribute("mission_data_version", metadata["source"])
        info.set_attribute("orbit", metadata["absolute_rev_number"])
        info.set_attribute("cycle", metadata["cycle_number"])
        info.set_attribute("mission_data_source",  metadata["product_name"])
        info.set_attribute("timeliness", self.cfg.timeliness_dict[str(product_info["sentinel3:timeliness"])])

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
            if mode == "sar":
                percent_value = 100.
            info.set_attribute("{}_mode_percent".format(mode), percent_value)
        info.set_attribute("open_ocean_percent", float(sral_info["sral:openOceanPercentage"]))

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
        utc_timestamp = num2date(self.nc.time_20_ku.values, units=self.nc.time_20_ku.units)
        self.l1.time_orbit.timestamp = utc_timestamp

        # Set the geolocation
        self.l1.time_orbit.set_position(
            self.nc.lon_20_ku.values,
            self.nc.lat_20_ku.values,
            self.nc.alt_20_ku.values,
            self.nc.orb_alt_rate_20_ku.values)

        # Set antenna attitude
        # NOTE: This are only available in 1Hz and need to be interpolated
        time_01, time_20 = self.nc.time_01.values, self.nc.time_20_ku.values
        pitch_angle_20, stat = self.interp_1Hz_to_20Hz(self.nc.off_nadir_pitch_angle_pf_01.values, time_01, time_20)
        roll_angle_20, stat = self.interp_1Hz_to_20Hz(self.nc.off_nadir_roll_angle_pf_01.values, time_01, time_20)
        yaw_angle_20, stat = self.interp_1Hz_to_20Hz(self.nc.off_nadir_yaw_angle_pf_01.values, time_01, time_20)
        self.l1.time_orbit.set_antenna_attitude(pitch_angle_20, roll_angle_20, yaw_angle_20)

    def _set_waveform_data_group(self):
        """
        Transfer of the waveform group to the Level-1 object. This includes
          1. the computation of waveform power in Watts
          2. the computation of the window delay in meter for each waveform bin
          3. extraction of the waveform valid flag
        :return: None
        """

        # Get the waveform
        # NOTE: The waveform is given in counts
        wfm_counts = self.nc.waveform_20_ku.values
        n_records, n_range_bins = wfm_counts.shape

        # Convert the waveform to power
        # TODO: This needs to be verified. Currently using the scale factor and documentation in netcdf unclear
        # From the documentation:
        # "This scaling factor represents the backscatter coefficient for a waveform amplitude equal to 1.
        #  It is corrected for AGC instrumental errors (agc_cor_20_ku) and internal calibration (sig0_cal_20_ku)"
        # NOTE: Make sure type of waveform is float and not double
        #       (double will cause issues with cythonized retrackers)
        wfm_power = np.ndarray(shape=wfm_counts.shape, dtype=np.float32)
        waveform_scale_factor = self.nc.scale_factor_20_ku.values
        for record in np.arange(n_records):
            wfm_power[record, :] = waveform_scale_factor[record] * wfm_counts[record, :].astype(float)

        # Get the window delay
        # "The tracker_range_20hz is the range measured by the onboard tracker
        #  as the window delay, corrected for instrumental effects and
        #  CoG offset"
        tracker_range_20hz = self.nc.tracker_range_20_ku.values
        wfm_range = np.ndarray(shape=wfm_counts.shape, dtype=np.float32)
        range_bin_index = np.arange(n_range_bins)
        for record in np.arange(n_records):
            wfm_range[record, :] = tracker_range_20hz[record] + \
                (range_bin_index*self.cfg.range_bin_width) - \
                (self.cfg.nominal_tracking_bin*self.cfg.range_bin_width)

        # Set the operation mode
        op_mode = self.nc.instr_op_mode_20_ku.values
        op_mode_translator = self.cfg.instr_op_mode_list
        radar_mode = np.array([op_mode_translator[int(val)] for val in op_mode]).astype("int8")

        # Set the waveform
        self.l1.waveform.set_waveform_data(wfm_power, wfm_range, radar_mode)

        # Get the valid flags
        # TODO: Find a way to get a valid flag
        # measurement_confident_flag = self.nc.flag_mcd_20_ku.values
        # valid_flag = measurement_confident_flag == 0
        # self.l1.waveform.set_valid_flag(valid_flag)

    def _set_range_correction_group(self):
        """
        Transfer the range corrections defined in the l1p config file to the Level-1 object
        NOTE: The range corrections are all in 1 Hz and must be interpolated to 20Hz
        :return: None
        """

        # Get the reference times for interpolating the range corrections from 1Hz -> 20Hz
        time_1Hz = self.nc.time_01.values
        time_20Hz = self.nc.time_20_ku.values

        # Loop over all range correction variables defined in the processor definition file
        for key in self.cfg.range_correction_targets.keys():
            var_name = self.cfg.range_correction_targets[key]
            variable_1Hz = getattr(self.nc, var_name)
            variable_20Hz, error_status = self.interp_1Hz_to_20Hz(variable_1Hz.values, time_1Hz, time_20Hz)
            if error_status:
                msg = "- Error in 20Hz interpolation for variable `%s` -> set only dummy" % var_name
                self.log.warning(msg)
            self.l1.correction.set_parameter(key, variable_20Hz)

    def _set_surface_type_group(self):
        """
        Transfer of the surface type flag to the Level-1 object
        NOTE: In the current state (TEST dataset), the surface type flag is only 1 Hz. A nearest neighbour
              interpolation is used to get the 20Hz surface type flag.
        :return: None
        """

        # Set the flag
        for key in ESA_SURFACE_TYPE_DICT.keys():
            flag = self.nc.surf_type_20_ku.values == ESA_SURFACE_TYPE_DICT[key]
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

    @property
    def empty(self):
        """
        Default return object, if nodata should be returned
        :return: Representation of an empty object (None)
        """
        return None