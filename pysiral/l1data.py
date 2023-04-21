# -*- coding: utf-8 -*-
"""
Created on Tue Jul 07 14:10:34 2015

@author: Stefan

L1bdata is a data container that unifies radar altimeter L1b orbit data
from different missions. It allows subsetting and merging of adjacent
orbit segments. L1bdata can be stored as a netCDF file, thus allowing
faster access to pre-processed subsets of RA orbit data for L2 processing.

The scope of the L1bdata container comprises:

---------
Metadata:
---------

- descriptors of RA source data
- period and geographical location
- processing history (subsetted, merged)
- software version

-------------
Waveform Data
-------------

- waveform echo power
  dimension: (n_records, n_range_bins)
- range for each range bin to the satellite in meters
  dimension: (n_records, n_range_bins)
- radar mode flag for each waveform:
    0: LRM
    1: SAR
    2: SIN
  (this is necessary for merging CryoSat-2 SAR and SIN adjacent orbit segments)
- summarizing flag from source data
    0: invalid
    1: valid
- optional: Additional named flags


----------------------
Time-Orbit Information
----------------------

- timestamp in UTC
- longitude, latitude (of satellite/nadir point)
- altitude (of satellite above WGS84 reference ellipsoid)

All parameter are of dimension (n_records).


-----------------
Range Corrections
-----------------

A list of range corrections (usually from RA source data files). The list
of correction is not predefined, but usally contains range corrections for:

- dry troposphere
- wet troposphere
- ionosphere
- inverse barometric / dynamic atmosphere
- ocean tide
- solid earth tide
- long period tide
- pole tide
- tidal loading

All parameter are of dimension (n_records) and of unit meter


----------
Classifier
----------

A list of optional named parameters that can be used for waveform
classification in the L2 processor. (e.g. stack parameter from the
CryoSat-2 l1b files)

All parameter are of dimension (n_records)


------------
Surface Type
------------



"""

import copy
from collections import OrderedDict
from typing import Any, List, Union

import numpy as np
import numpy.typing as npt
from cftime import num2pydate as cn2pyd
from pysiral.core.class_template import DefaultLoggingClass
from loguru import logger
from netCDF4 import Dataset, date2num
from scipy.spatial.transform import Rotation

from pysiral.core.config import RadarModes
from pysiral.core.flags import SurfaceType
from pysiral.core.output import NCDateNumDef

DATE2NUM_UNIT = "seconds since 1970-01-01 00:00:00.0"


class Level1bData(DefaultLoggingClass):
    """
    Unified L1b Data Class
    """
    data_groups = ["time_orbit", "correction", "classifier", "waveform", "surface_type"]

    def __init__(self):

        super(Level1bData, self).__init__(self.__class__.__name__)

        self.info = L1bMetaData()
        self.waveform = L1bWaveforms(self.info)
        self.time_orbit = L1bTimeOrbit(self.info)
        self.correction = L1bRangeCorrections(self.info)
        self.classifier = L1bClassifiers(self.info)
        self.surface_type = SurfaceType()

    def append(self,
               l1b_annex: "Level1bData",
               remove_overlap: bool = False,
               warn_if_temporal_offset_seconds: int = 10,
               raise_on_error: bool = False
               ) -> None:
        """
        Appends another l1b object to this one. The `l1b_annex` object is expected
        to have data after the self instance.

        :param l1b_annex:

        :return:
        """

        # Validity Checks
        timedelta_seconds = (l1b_annex.tcs - self.tce).total_seconds()
        if np.abs(timedelta_seconds) >= warn_if_temporal_offset_seconds:
            logger.warning(f"Appending l1 segment with time offset of {timedelta_seconds} seconds")

        # 1. The appending segment must not be a full subset of the current one
        time_coverage_str = f"self=[{l1b_annex.tcs}-{l1b_annex.tce}] annex=[{self.tcs}-{self.tce}]"
        if l1b_annex.tce < self.tce and l1b_annex.tcs > self.tcs:
            msg = f"l1.append: l1 segment to append is full subset of base segment: ({time_coverage_str})"
            logger.error(msg)
            if raise_on_error:
                raise ValueError(msg)
            logger.error(" -> l1 segment will not be appended")
            return

        # 2. There must be data after this one
        if l1b_annex.tce <= self.tce:
            msg = f"l1.append: No data after base segment: ({time_coverage_str})"
            logger.error(msg)
            if raise_on_error:
                raise ValueError(msg)
            logger.error(" -> l1 segment will not be appended")
            return

        # 3. There should not be an overlap
        if l1b_annex.tcs < self.tce:
            logger.warning(f"l1.append: Partial overlap between segments: ({time_coverage_str})")
            if remove_overlap:
                annex_non_overlap_idx = np.where(l1b_annex.time_orbit.timestamp > self.tce)[0]
                logger.warning("-> trimming appending segment")
                l1b_annex.trim_to_subset(annex_non_overlap_idx)

        # Append data in each data group
        for data_group in self.data_groups:
            this_data_group = getattr(self, data_group)
            annex_data_group = getattr(l1b_annex, data_group)
            this_data_group.append(annex_data_group)

        # Update the statistics
        self.info.set_attribute("is_merged_orbit", True)
        self.info.set_attribute("n_records", len(self.time_orbit.timestamp))
        mission_data_source = ";".join([self.info.mission_data_source,
                                        l1b_annex.info.mission_data_source])
        self.info.set_attribute("mission_data_source", mission_data_source)
        self.update_l1b_metadata()

    def trim_to_subset(self, subset_list: Union[List, npt.NDArray]) -> None:
        """ Create a subset from an index list """

        if len(subset_list) == 0:
            raise ValueError("subset list is emtpy")

        # Trim all data groups
        for data_group in self.data_groups:
            content = getattr(self, data_group)
            content.set_subset(subset_list)

        # Update metadata
        self.info.set_attribute("is_orbit_subset", True)
        self.info.set_attribute("n_records", len(subset_list))
        self.update_l1b_metadata()

    def apply_range_correction(self, correction):
        """  Apply range correction """
        # TODO: This method has no place here
        range_delta = self.correction.get_parameter_by_name(correction)
        if range_delta is None:
            # TODO: raise warning
            return

        # Check if NaN's, set values to zero and provide warning
        nans_indices = np.where(np.isnan(range_delta))[0]
        if len(nans_indices) > 0:
            range_delta[nans_indices] = 0.0
            logger.warning(f"NaNs encountered in range correction parameter: {correction}")

        self.waveform.add_range_delta(range_delta)

    def extract_subset(self, subset_list):
        """ Same as trim_to_subset, except returns a new l1bdata instance """
        if len(subset_list) == 0:
            return None
        l1b = copy.deepcopy(self)
        l1b.trim_to_subset(subset_list)
        return l1b

    def extract_region_of_interest(self, roi):
        """ Extracts data for a given region of interest definition """
        subset_list = roi.get_roi_list(self.time_orbit.longitude,
                                       self.time_orbit.latitude)
        if len(subset_list) > 0:
            l1b = copy.copy(self)
            l1b.trim_to_subset(subset_list)
        else:
            l1b = Level1bData()
        return l1b

    def detect_and_fill_gaps(self, gap_tolerance=0.02):
        """ Some radar altimeter input products are not provided with a
        regular time/sample spacing, e.g. by omitting short sections of the
        along-track data. This method is supposed to rectify that by
        1. Compute the nominal data repitition rate (in seconds)
        2. Detect gaps in the timestamp value
        3. Fills all arrays in all data groups with empty values
           (no interpolation) """

        # Get the time stamp and the time increment in seconds
        time = self.time_orbit.timestamp
        timedelta = np.ediff1d(time)
        timedelta_secs = [td.seconds + td.microseconds / 1e6 for td in timedelta]

        # Compute thresholds
        median_timedelta_secs = np.nanmedian(timedelta_secs)
        threshold_tolerance = median_timedelta_secs * gap_tolerance
        gap_sec_threshold = median_timedelta_secs + threshold_tolerance

        # Detect thresholds
        gap_start_indices = np.where(timedelta_secs > gap_sec_threshold)[0]

        # No gaps -> nothing to do
        if len(gap_start_indices) == 0:
            return

        # timedelta is [1:] because of np.ediff1d
        gap_start_indices += 1

        # indices_map' maps from the current indices to the corrected
        # indices and gap_indices contains the indices of the corrected
        # gap-filled arrays that need to be interpalated / filles with nodata
        # values.
        #
        # Example usage:
        #   timestamp_corrected[indices_map] = timestamp_old
        #   timestamp_corrected[gap_indices] = interpolated timestamp_old ...
        indices_map = np.arange(self.n_records)
        gap_indices = []
        for gap_start_index in gap_start_indices:
            gap_seconds = timedelta_secs[gap_start_index - 1]
            gap_width = int(np.round(gap_seconds / median_timedelta_secs))
            gap_indices.extend(
                np.arange(gap_width) + gap_start_index + len(gap_indices))
            indices_map[gap_start_index:] += gap_width

        # Get corrected n_records
        corrected_n_records = indices_map[-1] + 1

        # Update metadata
        self.info.set_attribute("n_records", corrected_n_records)

        # Loop over all data groups and fill gaps
        for data_group in self.data_groups:
            content = getattr(self, data_group)
            content.fill_gaps(corrected_n_records, gap_indices, indices_map)

    def update_l1b_metadata(self):
        self.update_data_limit_attributes()
        self.update_waveform_statistics()
        self.update_surface_type_statistics()
        self.update_region_name()

    def update_data_limit_attributes(self):
        """
        Set latitude/longitude and timestamp limits in the metadata container
        """

        info = self.info

        # Check if timestamp is monotonically increasing
        tdelta_dt = self.time_orbit.timestamp[1:]-self.time_orbit.timestamp[:-1]
        tdelta_secs = np.array([t.total_seconds() for t in tdelta_dt])
        if np.any(np.logical_and(tdelta_secs < 0.0, tdelta_secs >= -1.0)):
            logger.warning("- Found anomaly (small negative time step < -0.1 sec)")
        elif np.any(tdelta_secs < -1.0):
            logger.error("- Found anomaly (large negative time step > -1.0 sec)")
            import matplotlib.pyplot as plt
            plt.figure(dpi=150)
            plt.plot(tdelta_secs)
            plt.show()
            breakpoint()

        # time orbit group infos
        info.set_attribute("lat_min", np.nanmin(self.time_orbit.latitude))
        info.set_attribute("lat_max", np.nanmax(self.time_orbit.latitude))
        info.set_attribute("lon_min", np.nanmin(self.time_orbit.longitude))
        info.set_attribute("lon_max", np.nanmax(self.time_orbit.longitude))
        info.set_attribute("start_time", self.time_orbit.timestamp[0])
        info.set_attribute("stop_time", self.time_orbit.timestamp[-1])

    def update_waveform_statistics(self):
        """ Compute waveform metadata attributes """

        # waveform property infos (lrm, sar, sarin)
        radar_modes = RadarModes()
        radar_mode = self.waveform.radar_mode

        # Check if radar mode is none
        # (e.g. if only header information is parsed at this stage)
        if radar_mode is None:
            return

        # Compute the percentage of each radar mode in the l1b object
        nrecs_fl = float(self.n_records)
        for flag in range(radar_modes.num):
            is_this_radar_mode = np.where(radar_mode == flag)[0]
            radar_mode_percent = 100. * float(len(is_this_radar_mode)) / nrecs_fl
            attribute_name = f"{radar_modes.name(flag)}_mode_percent"
            self.info.set_attribute(attribute_name, radar_mode_percent)

    def update_surface_type_statistics(self):
        """ Re-calculate the open ocean percent """
        n_ocean_records = self.surface_type.get_by_name("ocean").num
        open_ocean_percent = 100. * float(n_ocean_records) / float(self.n_records)
        self.info.set_attribute("open_ocean_percent", open_ocean_percent)

    def update_region_name(self):
        """ Estimate the region (north/south/global) for metatdata class """

        lat_range = np.array([self.info.lat_min, self.info.lat_max])

        if np.amin(lat_range) > 0 and np.amax(lat_range) > 0:
            region_name = "north"
        elif np.amin(lat_range) < 0 and np.amax(lat_range) < 0:
            region_name = "south"
        else:
            region_name = "global"
        self.info.set_attribute("region_name", region_name)

    # TODO: Move to waveform class
    def reduce_waveform_bin_count(self, target_count: int, maxloc: float = 0.4) -> None:
        """
        Reduce the bin count of waveform power and range arrays.
        (e.g. for merging CryoSat-2 SAR [256 bins] and SIN [1024 bins])

        Creates a subset and updates the l1b.waveform container

        :param target_count: target number of waveform bins
          (needs to be smaller than full waveform bin count)
        :param maxloc:  preferred location of the maximum of the waveform in the subset

        :raises None:

        :return: None
        """

        # Extract original waveform
        orig_power, orig_range = self.waveform.power, self.waveform.range
        n_records, n_bins = orig_power.shape

        # Get the bin with the waveform maximum
        max_index = np.argmax(orig_power, axis=1)

        # Compute number of leading and trailing bins
        lead_bins = int(maxloc * target_count)
        trail_bins = target_count - lead_bins

        # Get the start/stop indices for each waveform
        start, stop = max_index - lead_bins, max_index + trail_bins

        # Create new arrays
        rebin_shape = (n_records, target_count)
        power = np.ndarray(shape=rebin_shape, dtype=orig_power.dtype)
        range_ = np.ndarray(shape=rebin_shape, dtype=orig_range.dtype)

        # Validity check
        overflow = np.where(stop > n_bins)[0]
        if len(overflow) > 0:
            offset = n_bins - stop[overflow]
            stop[overflow] += offset
            start[overflow] += offset

        underflow = np.where(start < 0)[0]
        if len(underflow) > 0:
            offset = start[underflow]
            stop[underflow] -= offset
            start[underflow] -= offset

        # Extract the waveform with reduced bin count
        for i in np.arange(n_records):
            power[i, :] = orig_power[i, start[i]:stop[i]]
            range_[i, :] = orig_range[i, start[i]:stop[i]]

        # Push to waveform container
        self.waveform.set_waveform_data(power, range_, self.radar_modes)

    # TODO: Move to waveform class
    def increase_waveform_bin_count(self, target_count: int) -> None:
        """
        Increase the bin count of waveform power and range arrays.
        (e.g. for merging CryoSat-2 LRM [128 bins] and SAR [256 bins])

        Creates a subset and updates the l1b.waveform container

        :param target_count: target number of waveform bins
                             (needs to be bigger than full waveform bin count)

        :raises None:

        :return: None
        """
        # Extract original waveform
        orig_power, orig_range = self.waveform.power, self.waveform.range
        n_records, n_bins = orig_power.shape

        # Add the zero bins at the beginning of the range window
        pwr = np.full((n_records, target_count), 0.0)
        pwr[:, 0:n_bins] = orig_power

        # Extend the range
        rng = np.full((n_records, target_count), 0.0)
        rng[:, 0:n_bins] = orig_range

        n_new_bins = target_count - n_bins
        approx_bins_size = rng[0, 1] - rng[0, 0]
        artificial_range = np.arange(1, n_new_bins + 1) * approx_bins_size
        rng[:, n_bins:] = np.tile(artificial_range, (n_records, 1))
        rng[:, n_bins:] += np.tile(rng[:, n_bins - 1], (n_new_bins, 1)).transpose()

        # Push to waveform container
        self.waveform.set_waveform_data(pwr, rng, self.radar_modes)

    def get_parameter_by_name(self, data_group: str, parameter_name: str) -> Union[None, np.ndarray]:
        """
        API method to retrieve any parameter from any data group

        :param data_group:
        :param parameter_name:

        :return:
        """
        try:
            data_group = getattr(self, data_group)
            return getattr(data_group, parameter_name)
        except AttributeError:
            return None

    def set_parameter_by_name(self, data_group_name: str, parameter_name: str, value: np.ndarray) -> None:
        """ API method to set any parameter in any data group """
        # Sanity check
        if len(value) != self.n_records:
            raise ValueError(f"value for {data_group_name}.{parameter_name} has wrong shape: {str(value.shape)}")
        try:
            # Get data group
            data_group = getattr(self, data_group_name)
            # Update data group
            setattr(data_group, parameter_name, value)
            # Update l1b container
            setattr(self, data_group_name, data_group)
        except AttributeError as e:
            raise ValueError(f"Could not set value for {data_group_name}.{parameter_name}") from e

    @property
    def n_records(self):
        return self.info.n_records

    @property
    def tcs(self):
        return self.info.start_time

    @property
    def tce(self):
        return self.info.stop_time

    @property
    def radar_modes(self):
        radar_modes = RadarModes()
        radar_mode_flag_list = np.unique(self.waveform.radar_mode)
        radar_mode_list = [radar_modes.name(radar_mode_flag) for radar_mode_flag in radar_mode_flag_list]
        return ";".join(radar_mode_list)


class L1bdataNCFile(Level1bData):

    def __init__(self, filename):

        super(L1bdataNCFile, self).__init__()
        self.filename = filename
        self.nc = None
        self.time_def = NCDateNumDef()
        self.ncattrs_ignore_list = ['_NCProperties']

    def parse(self):
        """ populated the L1b data container from the l1bdata netcdf file """
        self.nc = Dataset(self.filename, "r")
        self.nc.set_auto_scale(False)
        self._import_metadata()
        self._import_timeorbit()
        self._import_waveforms()
        self._import_corrections()
        self._import_surface_type()
        self._import_classifier()
        self.nc.close()

    def _import_metadata(self):
        """
        transfers l1b metadata attributes
        (stored as global attributes in l1bdata netCDF files)
        """
        for attribute_name in self.nc.ncattrs():
            if attribute_name in self.ncattrs_ignore_list:
                continue
            attribute_value = getattr(self.nc, attribute_name)
            # Convert timestamps back to datetime objects
            if attribute_name in ["start_time", "stop_time"]:
                attribute_value = cn2pyd(attribute_value, self.time_def.units, calendar=self.time_def.calendar)
            # Convert flags (integers back to bool)
            if attribute_name in ["is_orbit_subset", "is_merged_orbit"]:
                attribute_value = bool(attribute_value)
            self.info.set_attribute(attribute_name, attribute_value)

    def _import_timeorbit(self):
        """
        transfers l1b timeorbit group
        (timeorbit datagroup in l1bdata netCDF files)
        """
        # Get the datagroup
        datagroup = self.nc.groups["time_orbit"]

        # Set satellite position data (measurement is nadir)
        self.time_orbit.set_position(
            datagroup.variables["longitude"][:],
            datagroup.variables["latitude"][:],
            datagroup.variables["altitude"][:])

        antenna_angles = {}
        for angle in ["pitch", "roll", "yaw"]:
            try:
                value = datagroup.variables[f"antenna_{angle}"][:]
            except KeyError:
                value = np.full(self.time_orbit.longitude.shape, 0.0)
            antenna_angles[angle] = value

        # Set satellite position data (measurement is nadir)
        self.time_orbit.set_antenna_attitude(
            antenna_angles["pitch"],
            antenna_angles["roll"],
            antenna_angles["yaw"])

        # Convert the timestamp to datetimes
        self.time_orbit.timestamp = cn2pyd(
            datagroup.variables["timestamp"][:],
            self.time_def.units,
            calendar=self.time_def.calendar)

    def _import_waveforms(self):
        """
        transfers l1b waveform group
        (waveform datagroup in l1bdata netCDF files)
        """
        # Get the datagroup
        datagroup = self.nc.groups["waveform"]

        # Set waveform (measurement is nadir)
        self.waveform.set_waveform_data(
            datagroup.variables["power"][:],
            datagroup.variables["range"][:],
            datagroup.variables["radar_mode"][:])
        # Set the valid flag
        is_valid = datagroup.variables["is_valid"][:].astype(bool)
        self.waveform.set_valid_flag(is_valid)

    def _import_corrections(self):
        """
        transfers l1b corrections group
        (waveform corrections in l1bdata netCDF files)
        """
        # Get the datagroup
        datagroup = self.nc.groups["correction"]
        # Loop over parameters
        for key in datagroup.variables.keys():
            variable = np.array(datagroup.variables[key][:])
            self.correction.set_parameter(key, variable)

    def _import_surface_type(self):
        """
        transfers l1b surface_type group
        (waveform corrections in l1bdata netCDF files)
        """
        # Get the datagroup
        datagroup = self.nc.groups["surface_type"]
        self.surface_type.set_flag(datagroup.variables["flag"][:])

    def _import_classifier(self):
        """
        transfers l1b corrections group
        (waveform corrections in l1bdata netCDF files)
        """
        # Get the data group
        datagroup = self.nc.groups["classifier"]
        # Loop over parameters
        for key in datagroup.variables.keys():
            variable = np.array(datagroup.variables[key][:])
            self.classifier.add(variable, key)


class L1bMetaData(object):
    """
    Container for L1B Metadata information
    (see property attribute_list for a list of attributes)
    """

    def __init__(self):
        """
        Class containing a specific set of metadata attributes for l1b/l1p data
        """

        # Init all fields
        self._attribute_list = [
            "pysiral_version", "mission", "mission_data_version",
            "mission_sensor", "mission_data_source", "n_records", "orbit", "rel_orbit",
            "cycle", "sar_mode_percent", "lrm_mode_percent", "sin_mode_percent",
            "is_orbit_subset", "is_merged_orbit", "start_time", "stop_time",
            "region_name", "lat_min", "lat_max", "lon_min", "lon_max",
            "open_ocean_percent", "timeliness"]

        self._attrs = {attr_name: None for attr_name in self._attribute_list}

        # Set some fields to False (instead of none)
        self._attrs["orbit"] = 999999
        self._attrs["is_orbit_subset"] = False
        self._attrs["is_merged_orbit"] = False
        self._attrs["n_records"] = -1

    def __repr__(self):
        output = "pysiral.L1bdata object:\n"
        for field in self._attribute_list:
            output += "%22s: %s" % (field, getattr(self, field))
            output += "\n"
        return output

    def __getattr__(self, item):
        """
        Modify the attribute getter to provide a shortcut to the data content
        :param item: Name of the parameter
        :return:
        """
        if item == "__setstate__":
            raise AttributeError(item)
        if item in self._attrs:
            return self._attrs[item]
        else:
            raise AttributeError(f"L1BMetadata does not have the attribute {item}")

    @property
    def attribute_list(self):
        return self._attribute_list

    @property
    def attdict(self):
        """ Return attributes as dictionary (e.g. for netCDF export) """
        return {field: getattr(self, field) for field in self.attribute_list}

    @property
    def hemisphere(self):
        hemisphere = "global"
        if self.lat_min > 0. and self.lat_max > 0.:
            hemisphere = "north"
        if self.lat_min < 0. and self.lat_max < 0.0:
            hemisphere = "south"
        return hemisphere

    @property
    def year(self):
        return "%04g" % self.start_time.year

    @property
    def month(self):
        return "%02g" % self.start_time.month

    def set_attribute(self, tag, value):
        if tag not in self.attribute_list:
            raise ValueError("Unknown attribute: ", tag)
        self._attrs[tag] = value

    def check_n_records(self, n_records: int) -> None:
        """
        First time a data set is set: Store number of records as reference

        :param n_records: Number of records

        :return: None

        :raises: ValueError
        """

        if self._attrs["n_records"] == -1:
            self._attrs["n_records"] = n_records

        elif n_records != self.n_records:
            msg = f"n_records mismatch, len must be: {self.n_records} (was {n_records})"
            raise ValueError(msg)


class L1bTimeOrbit(object):
    """ Container for Time and Orbit Information of L1b Data """

    def __init__(self, info, is_evenly_spaced=True):
        self._info = info  # Pointer to metadata container
        self._timestamp = None
        self._longitude = None
        self._latitude = None
        self._altitude = None
        self._altitude_rate = None
        self._antenna_pitch = None
        self._antenna_roll = None
        self._antenna_yaw = None
        self._antenna_mispointing = None
        self._orbit_flag = None
        self._is_evenly_spaced = is_evenly_spaced

    @property
    def longitude(self):
        return np.array(self._longitude)

    @property
    def latitude(self):
        return np.array(self._latitude)

    @property
    def altitude(self):
        return np.array(self._altitude)

    @property
    def altitude_rate(self):
        return np.array(self._altitude_rate)

    @property
    def antenna_pitch(self):
        return np.array(self._antenna_pitch)

    @property
    def antenna_roll(self):
        return np.array(self._antenna_roll)

    @property
    def antenna_yaw(self):
        return np.array(self._antenna_yaw)

    @property
    def antenna_mispointing(self):
        return np.array(self._antenna_yaw)

    @property
    def orbit_flag(self):
        return np.array(self._orbit_flag)

    @property
    def timestamp(self):
        return np.array(self._timestamp)

    @timestamp.setter
    def timestamp(self, value):
        if self._info is not None:
            self._info.check_n_records(len(value))
        self._timestamp = value

    @property
    def parameter_list(self):
        return ["timestamp", "longitude", "latitude", "altitude", "altitude_rate",
                "antenna_pitch", "antenna_roll", "antenna_yaw", "antenna_mispointing", "orbit_flag"]

    @property
    def geolocation_parameter_list(self):
        return ["longitude", "latitude", "altitude", "altitude_rate",
                "antenna_pitch", "antenna_roll", "antenna_yaw", "antenna_mispointing"]

    @property
    def dimdict(self):
        """ Returns dictionary with dimensions"""
        return OrderedDict([("n_records", len(self._timestamp))])

    @property
    def is_evenly_spaced(self):
        return self._is_evenly_spaced

    def set_position(self, longitude, latitude, altitude, altitude_rate=None):
        # Check dimensions
        if self._info is not None:
            self._info.check_n_records(len(longitude))
            self._info.check_n_records(len(latitude))
            self._info.check_n_records(len(altitude))
            if altitude_rate is not None:
                self._info.check_n_records(len(altitude_rate))

        # All fine => set values
        self._longitude = longitude
        self._latitude = latitude
        self._altitude = altitude

        # Parameter that were added later
        dummy_val = np.full(self.longitude.shape, np.nan)
        self._altitude_rate = altitude_rate if altitude_rate is not None else dummy_val

        # Set a dummy value for pitch, roll & yaw for backward compability
        if self.antenna_pitch is None:
            self.set_antenna_attitude(dummy_val, dummy_val, dummy_val)

        # Compute orbit flag (0: ascending, 1: descending)
        latitude_rate = latitude[1:] - latitude[:-1]
        try:
            latitude_rate = np.insert(latitude_rate, 0, latitude_rate[0])
            self._orbit_flag = (latitude_rate < 0).astype(int)
        except IndexError:
            self._orbit_flag = np.full(latitude.shape, -1)

    def set_antenna_attitude(self, pitch, roll, yaw, mispointing=None):
        # Check dimensions
        if self._info is not None:
            self._info.check_n_records(len(pitch))
            self._info.check_n_records(len(roll))
            self._info.check_n_records(len(yaw))
            if self._antenna_mispointing is not None:
                self._info.check_n_records(len(mispointing))

        # All fine => set values
        self._antenna_pitch = pitch
        self._antenna_roll = roll
        self._antenna_yaw = yaw
        self._antenna_mispointing = mispointing if mispointing is not None else self.mispointing_from_angles(pitch, roll, yaw)

    def append(self, annex):
        for parameter in self.parameter_list:
            this_data = getattr(self, f"_{parameter}")
            annex_data = getattr(annex, parameter)
            this_data = np.append(this_data, annex_data)
            setattr(self, f"_{parameter}", this_data)

    def set_subset(self, subset_list):
        for parameter in self.parameter_list:
            data = getattr(self, f"_{parameter}")
            data = data[subset_list]
            setattr(self, f"_{parameter}", data)

    def fill_gaps(self, corrected_n_records, gap_indices, indices_map):
        """ API gap filler method. Note: It is assumed that this method is
        only evoked for filling small gaps. Therefore, we use simple linear
        interpolation for the parameters of the time orbit group """

        # Set the geolocation parameters first (lon, lat, alt)
        geoloc_parameters = []
        corrected_indices = np.arange(corrected_n_records)
        for parameter_name in self.geolocation_parameter_list:
            data_old = getattr(self, parameter_name)
            data_corr = np.interp(corrected_indices, indices_map, data_old)
            geoloc_parameters.append(data_corr)
        self.set_position(*geoloc_parameters)

        # Update the timestamp
        time_old_num = date2num(self.timestamp, DATE2NUM_UNIT)
        time_num = np.interp(corrected_indices, indices_map, time_old_num)
        self.timestamp = cn2pyd(time_num, DATE2NUM_UNIT)

    def get_parameter_by_name(self, name):
        try:
            return getattr(self, name)
        except AttributeError:
            return None

    @staticmethod
    def mispointing_from_angles(pitch_deg: npt.NDArray,
                                roll_deg: npt.NDArray,
                                heading_deg: npt.NDArray
                                ) -> npt.NDArray:
        """
        Compute the mispointing (angle between -z direction in spacecraft frame and true nadir)
        from pitch/roll/heading, assuming rotation of spacecraft is around antenna.

        :param pitch_deg: pitch angles (rotation around y-axis)
        :param roll_deg: roll angles (rotation around x-axis)
        :param heading_deg: true heading (rotation around z-axis)

        :return: mispointing angle in degrees
        """

        # Array of [0, 0, -1] vector (down in satellite reference frame)
        x_arr = np.repeat(np.array([[0, 0, -1]]), pitch_deg.shape[0], axis=0)

        # Rotate the down vector with pitch, roll, heading
        r = Rotation.from_rotvec(np.c_[roll_deg, pitch_deg, heading_deg], degrees=True)
        y_arr = r.apply(x_arr)

        # Mis-pointing is the angle between spacecraft down and nadir
        return np.rad2deg(np.arccos([np.dot(x, y) for x, y in zip(x_arr, y_arr)]))

    def __getstate__(self):
        return self.__dict__

    def __setstate__(self, d):
        self.__dict__.update(d)


class L1bRangeCorrections(object):
    """ Container for Range Correction Information """

    def __init__(self, info):
        self._info = info  # Pointer to Metadata object
        self._parameter_list = []

    def set_parameter(self, tag, value):
        self._info.check_n_records(len(value))
        setattr(self, tag, value)
        if tag not in self._parameter_list:
            self._parameter_list.append(tag)

    @property
    def parameter_list(self):
        return self._parameter_list

    @property
    def n_records(self):
        parameter, name = self.get_parameter_by_index(0)
        return len(parameter)

    @property
    def dimdict(self):
        """ Returns dictionary with dimensions"""
        return OrderedDict([("n_records", self.n_records)])

    def get_parameter_by_index(self, index):
        name = self._parameter_list[index]
        return getattr(self, name), name

    def get_parameter_by_name(self, name):
        try:
            return getattr(self, name)
        except AttributeError:
            return None

    def append(self, annex):
        for parameter in self.parameter_list:
            this_data = getattr(self, parameter)
            annex_data = getattr(annex, parameter)
            this_data = np.append(this_data, annex_data)
            setattr(self, parameter, this_data)

    def set_subset(self, subset_list):
        for parameter in self.parameter_list:
            data = getattr(self, parameter)
            data = data[subset_list]
            setattr(self, parameter, data)

    def fill_gaps(self, corrected_n_records, gap_indices, indices_map):
        """ API gap filler method. Note: Gaps will be filled with
        the nodata=0.0 value"""

        for parameter_name in self.parameter_list:
            data_corr = np.full(corrected_n_records, 0.0)
            data_old = self.get_parameter_by_name(parameter_name)
            data_corr[indices_map] = data_old
            self.set_parameter(parameter_name, data_corr)


class L1bClassifiers(object):
    """ Containier for parameters that can be used as classifiers """

    def __init__(self, info):
        self._info = info  # Pointer to Metadate object
        # Make a pre-selection of different classifier types
        self._list = {
            "surface_type": [],
            "warning": [],
            "error": []}

    def add(self, value, name, classifier_type="surface_type"):
        """ Add a parameter for a given classifier type """
        setattr(self, name, np.array(value))
        if name not in self._list[classifier_type]:
            self._list[classifier_type].append(name)

    @property
    def parameter_list(self):
        parameter_list = []
        for key in self._list.keys():
            parameter_list.extend(self._list[key])
        return parameter_list

    @property
    def n_records(self):
        parameter_list = self.parameter_list
        return 0 if len(parameter_list) == 0 else len(getattr(self, parameter_list[0]))

    @property
    def dimdict(self):
        """ Returns dictionary with dimensions"""
        return OrderedDict([("n_records", self.n_records)])

    def has_parameter(self, parameter_name):
        return parameter_name in self.parameter_list

    def get_parameter(self, parameter_name: str, raise_on_error: bool = False) -> Union[np.ndarray, None]:
        if raise_on_error:
            return getattr(self, parameter_name, None)
        else:
            return getattr(self, parameter_name)

    def append(self, annex):
        for parameter in self.parameter_list:
            this_data = getattr(self, parameter)
            annex_data = getattr(annex, parameter)
            this_data = np.append(this_data, annex_data)
            setattr(self, parameter, this_data)

    def set_subset(self, subset_list):
        for parameter in self.parameter_list:
            data = getattr(self, parameter)
            data = data[subset_list]
            setattr(self, parameter, data)

    def fill_gaps(self, corrected_n_records, gap_indices, indices_map):
        """ API gap filler method. Note: Gaps will be filled with
        the nodata=nan value """

        for parameter_name in self.parameter_list:
            data_corr = np.full(corrected_n_records, np.nan)
            data_old = self.get_parameter(parameter_name)
            data_corr[indices_map] = data_old
            self.add(data_corr, parameter_name)

    def __getattr__(self, item: str) -> Any:
        """
        Direct attribute access to the cfg dictionary

        :param item:
        :return:
        """
        if self.has_parameter(item):
            return self.get_parameter(item)
        else:
            raise AttributeError(f"attribute {item} not found in classifier container")

    def __getstate__(self):
        return self.__dict__

    def __setstate__(self, d):
        self.__dict__.update(d)


class L1bWaveforms(object):
    """ Container for Echo Power Waveforms """

    _valid_radar_modes = ["lrm", "sar", "sin"]
    _parameter_list = ["power", "range", "radar_mode", "is_valid", "classification_flag"]
    _attribute_list = ["echo_power_unit"]

    def __init__(self, info):
        self._info = info  # Pointer to Metadate object
        # Attributes
        self.echo_power_unit = None
        self.radar_mode_def = RadarModes()
        # Parameter
        self._power = None
        self._range = None
        self._radar_mode = None
        self._is_valid = None
        self._classification_flag = None

    @property
    def power(self):
        return np.copy(self._power)

    @property
    def classification_flag(self):
        if self._classification_flag is None:
            return np.full(self.n_records, -1, dtype=int)
        return np.copy(self._classification_flag)

    @property
    def range(self):
        return np.copy(self._range)

    @property
    def radar_mode(self):
        return np.copy(self._radar_mode)

    @property
    def is_valid(self):
        return np.copy(self._is_valid)

    @property
    def parameter_list(self):
        return list(self._parameter_list)

    @property
    def n_range_bins(self):
        return self._get_wfm_shape(1)

    @property
    def n_records(self):
        return self._get_wfm_shape(0)

    @property
    def radar_modes(self):
        if self._radar_mode is None:
            return "none"
        flags = np.unique(self._radar_mode)
        return [self.radar_mode_def.get_name(flag) for flag in flags]

    @property
    def dimdict(self):
        """ Returns dictionary with dimensions"""
        shape = np.shape(self._power)
        return OrderedDict([("n_records", shape[0]), ("n_bins", shape[1])])

    def set_waveform_data(self, power, range, radar_mode, classification_flag=None):
        """
        Set the waveform data
        :param power:
        :param range:
        :param radar_mode:
        :param classification_flag:
        :return:
        """
        # Validate input
        if power.shape != range.shape:
            raise ValueError("power and range must be of same shape", power.shape, range.shape)
        if len(power.shape) != 2:
            raise ValueError("power and range arrays must be of dimension (n_records, n_bins)")

            # Validate number of records
        self._info.check_n_records(power.shape[0])

        # Assign values
        self._power = power
        self._range = range

        # Create radar mode arrays
        if type(radar_mode) is str and radar_mode in self._valid_radar_modes:
            mode_flag = self.radar_mode_def.get_flag(radar_mode)
            self._radar_mode = np.repeat(mode_flag, self.n_records).astype(np.byte)
        elif len(radar_mode) == self._info.n_records:
            self._radar_mode = radar_mode.astype(np.int8)
        else:
            raise ValueError("Invalid radar_mode: ", radar_mode)

        # Set valid flag (assumed to be valid for all waveforms)
        # Flag can be set separately using the set_valid_flag method
        if self._is_valid is None:
            self._is_valid = np.ones(shape=self.n_records, dtype=bool)

    def set_valid_flag(self, valid_flag):
        # Validate number of records
        self._info.check_n_records(len(valid_flag))
        self._is_valid = valid_flag

    def set_classification_flag(self, classification_flag):
        """
        Add or update the waveform classification flag
        :param classification_flag: intarray with shape (n_records, n_range_bins)
        :return:
        """
        # Validate number of records
        if classification_flag.shape != self.power.shape:
            raise ValueError(f"Invalid dimensions: {classification_flag.shape} [{self.power.shape}]")
        self._classification_flag = classification_flag

    def append(self, annex):
        self._power = np.concatenate((self._power, annex.power), axis=0)
        self._range = np.concatenate((self._range, annex.range), axis=0)
        self._radar_mode = np.append(self._radar_mode, annex.radar_mode)
        self._is_valid = np.append(self._is_valid, annex.is_valid)

    def set_subset(self, subset_list):
        self._power = self._power[subset_list, :]
        self._range = self._range[subset_list, :]
        self._radar_mode = self._radar_mode[subset_list]
        self._is_valid = self._is_valid[subset_list]

    def add_range_delta(self, range_delta):
        """
        Add a range delta to all range bins
        :param range_delta:
        :return:
        """
        range_delta_reshaped = np.repeat(range_delta, self.n_range_bins)
        range_delta_reshaped = range_delta_reshaped.reshape(self.n_records, self.n_range_bins)
        self._range += range_delta_reshaped

    def fill_gaps(self, corrected_n_records, gap_indices, indices_map):
        """ API gap filler method. Note: Gaps will be filled with
        custom values for each parameter, see below"""

        # is_valid flag: False for gaps
        is_valid = np.full(corrected_n_records, False)
        is_valid[indices_map] = self.is_valid
        self.set_valid_flag(is_valid)

        # Power/range: set gaps to nan
        power = np.full((corrected_n_records, self.n_range_bins), np.nan)
        power[indices_map, :] = self.power
        range_ = np.full((corrected_n_records, self.n_range_bins), np.nan)
        range_[indices_map, :] = self.range

        # Radar map: set gaps to lrm
        radar_mode = np.full(corrected_n_records, 1, dtype=self.radar_mode.dtype)
        radar_mode[indices_map] = self.radar_mode

        # And set new values
        self.set_waveform_data(power, range_, radar_mode)

    def _get_wfm_shape(self, index):
        shape = np.shape(self._power)
        return shape[index]
