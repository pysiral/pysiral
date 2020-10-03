# -*- coding: utf-8 -*-

"""
@author: Stefan Hendricks

pysiral module for estimating sea surface height (ssh) respectively sea level anomaly (sla) from along-track
radar altimeter data. The classes are designed to work with the Level-2 Processor, therefore need to be
child classes of pysiral.l2proc.procsteps.Level2ProcessorStep.

NOTES:

    1. The convention of the pysiral Level-2 processor is to compute and sla (with ssh = sla + mss)
       to the Level-2 data object. Thus, the auxiliary data set mean sea surface (mss) is a mandatory
       auxiliary data set for the functionality of all classes in this module.

"""


import numpy as np

from pysiral.l2proc.procsteps import Level2ProcessorStep
from pysiral.filter import (fill_nan, idl_smooth)


class SLASmoothedLinear(Level2ProcessorStep):
    """
    Default implemetation of a smoothed interpolation of ssh tiepoints (with various filter options).
    This class will compute the sea level anomaly (sla)
    """

    def __init__(self, *args, **kwargs):
        """
        Init the class. Options will be passed to parent class
        (pysiral.l2proc.procsteps.Level2ProcessorStep)
        :param args:
        :param kwargs:
        """
        super(SLASmoothedLinear, self).__init__(*args, **kwargs)

        # Init Properties
        self.ssh_tiepoints = None      # A boolean flag which elevation obs is a valid ssh tiepoint
        self.sla_raw = None            # The raw (unfiltered) sea level anomaly
        self.sla = None                # The final (filtered) sea level anomaly
        self.sla_uncertainty = None    # Uncertainty of the sea level anomaly

    def execute_procstep(self, l1b, l2):
        """
        Mandatory Level-2 processor method that will execute the processing step
        and modify the L2 data object in-place.
        This method will interpolate ssh tiepoints given by lead (+ ocean) elevations
        and smooth the result. Filter options are available to filter unreasonable
        sla values.
        :param l1b:
        :param l2:
        :return:
        """

        # Step 1: A linear interpolation between lead elevations
        # -> will add properties `sla_raw`, `ssh_tiepoints` and `sla` to the instance
        self.smoothed_linear_interpolation_between_tiepoints(l2)

        # Step 2: Compute sea level anomaly uncertainty
        # -> will add properties `sla_uncertainty`
        self.calculate_sla_uncertainty(l2)

        # Step 3 (optional): Filter small marine segments surrounded by land
        # Note: This intends to remove small segments in fjords/channels for which
        #       the SLA computation is very likely not trustworthy
        if "marine_segment_filter" in self.cfg:
            self.marine_segment_filter(l2)

        # Step 4 (optional): Filter SLA segments that are far away from the next SSH
        #    tie point
        if "tiepoint_maxdist_filter" in self.cfg:
            self.tiepoint_maxdist_filter(l2)

        # Step 5: Modify the Level-2 data container with the result in-place
        l2.sla.set_value(self.sla)
        l2.sla.set_uncertainty(self.sla_uncertainty)

        # Return the error status
        error_status = np.isnan(l2.sla[:])
        return error_status

    def smoothed_linear_interpolation_between_tiepoints(self, l2):
        """
        The main SLA computation method in this class
        :param l2: Level-2 data container
        :return: None
        """

        # Collect the first estimate of ssh tie points
        # NOTE: The use of ocean waveforms is optional and needs to activated
        #       in the options dictionary of the Level-2 processor definition file
        self.ssh_tiepoints = l2.surface_type.lead.indices
        if self.cfg.options.use_ocean_wfm:
            self.ssh_tiepoints.append(l2.surface_type.ocean.indices)
            self.ssh_tiepoints = np.sort(self.ssh_tiepoints)

        # Get initial elevation at tie point locations
        mss_frb = l2.elev - l2.mss

        # Remove ssh tie points from the list if their elevation
        # corrected by the median offset of all tie points from the mss
        # exceeds a certain threshold
        if self.cfg.options.pre_filtering:

            # Startup
            index_dict = np.arange(l2.surface_type.lead.num)
            threshold = self.cfg.options.pre_filter_maximum_mss_median_offset

            # Compute the mean distance to the mss
            tiepoint_mss_distance = mss_frb[self.ssh_tiepoints]
            valid_points = np.where(np.isfinite(tiepoint_mss_distance))[0]
            tiepoint_mss_distance = tiepoint_mss_distance[valid_points]
            median_mss_offset = np.median(tiepoint_mss_distance)

            # Filter points for outliers
            offset = np.abs(mss_frb[self.ssh_tiepoints] - median_mss_offset)
            valid = np.where(offset < threshold)
            self.ssh_tiepoints = self.ssh_tiepoints[index_dict[valid]]

        # Compute the first SLA estimate
        self.sla_raw = np.ndarray(shape=l2.n_records)*np.nan
        self.sla_raw[self.ssh_tiepoints] = mss_frb[self.ssh_tiepoints]
        non_tiepoints = np.where(np.isnan(self.sla_raw))

        # Filtered raw values (python implementation of the CS2AWI IDL code)
        # Use custom implementation of IDL SMOOTH:
        # idl_smooth(x, w) equivalent to SMOOTH(x, w, /edge_truncate, /nan)
        sla_filter1 = idl_smooth(self.sla_raw, self.filter_width)

        # Leave only the original ssh tie points
        sla_filter1[non_tiepoints] = np.nan

        # Fill nans with linear interpolation and constant values at borders
        # python: fill_nan(x) = IDL: FILL_NAN(x, /NEIGHBOUR)
        sla_filter2 = fill_nan(sla_filter1)

        # Final smoothing
        sla = idl_smooth(sla_filter2, self.filter_width)
        self.sla = sla

    def calculate_sla_uncertainty(self, l2):
        """
        Components that add to sea surface anomaly uncertainty
            - mss uncertainty (if known)
            - uncertainty of lead elevations
            - distance to next lead tiepoint
        :param l2: Level-2 data container
        :return: None
        """

        # short cuts to options
        max_distance = self.cfg.options.uncertainty_tiepoints_distance_max
        sla_unc_min = self.cfg.options.uncertainty_minimum
        sla_unc_max = self.cfg.options.uncertainty_maximum

        # get tie point distance
        tiepoint_distance = self.get_tiepoint_distance(l2)

        # Compute the influence of distance to next tie points
        # in the range of 0: minimum influence to 1: maximum influence
        # It is assumed that the uncertainty has a quadratic dependence
        # on tie point distance
        tiepoint_distance_scalefact = tiepoint_distance / max_distance
        above_distance_limit = np.where(tiepoint_distance_scalefact > 1.)[0]
        tiepoint_distance_scalefact[above_distance_limit] = 1.
        tiepoint_distance_scalefact = tiepoint_distance_scalefact**2.

        # Compute the sla uncertainty based on a min/max approach scaled by factor
        sla_unc_range = sla_unc_max - sla_unc_min
        sla_unc = sla_unc_min + sla_unc_range * tiepoint_distance_scalefact

        # Save result to instance
        self.sla_uncertainty = sla_unc

    def marine_segment_filter(self, l2):
        """ Check all sections divided by land masses for reliable
        information content """

        filter_options = self.cfg.options.marine_segment_filter
        minimum_lead_number = filter_options.minimum_lead_number
        footprint_size = self.cfg.options.smooth_filter_width_footprint_size
        section_prop = {"i0": 0.0, "i1": 0.0,
                        "width": 0.0, "n_tiepoints": 0,
                        "land_before": (9999.0, 0),
                        "land_after": (9999.0, 0)}

        # Find sea ice clusters
        land = l2.surface_type.land

        # No land -> nothing to do
        if land.num == 0:
            return

        # Get indices for land sections
        lead_flag = l2.surface_type.lead.flag
        land_flag = land.flag.astype(int)
        land_start = np.where(np.ediff1d(land_flag) > 0)[0]
        land_stop = np.where(np.ediff1d(land_flag) < 0)[0]

        # It is assumed here, that the l1b orbit segment never starts
        # or end with land. Thus the number of land start and land stop
        # events need to be identical.
        if len(land_start) != len(land_stop):
            code = "l2-crop-error"
            msg = "l2 segments either starts or ends with land"
            self.error.add_error(code, msg)
            return

        # Add artificial large land sections on beginning and end of profile
        n_marine_segments = len(land_start) + 1
        n = l2.n_records
        land_start = np.concatenate(([-1000], land_start, [n-1]))
        land_stop = np.concatenate(([-1], land_stop, [n+1000]))

        # Loop over marine segments and collect information
        marine_segments = []
        for i in np.arange(n_marine_segments):

            marine_segment = section_prop.copy()

            # Get the start stop indices for marine section
            i0 = land_stop[i]+1
            i1 = land_start[i+1]
            marine_segment["i0"] = i0
            marine_segment["i1"] = i1
            marine_segment["width"] = (i1-i0) * footprint_size

            # get the number of leads
            marine_section_indices = np.arange(i0, i1+1)
            n_tiepoints = np.where(lead_flag[marine_section_indices])[0].size
            marine_segment["n_tiepoints"] = n_tiepoints

            if marine_segment["n_tiepoints"] < minimum_lead_number:
                self.sla[marine_section_indices] = np.nan

            marine_segments.append(marine_segment)

    def tiepoint_maxdist_filter(self, l2):
        """
        A filter that does not removes sla values which distance to
        the next ssh tiepoint exceeds a defined threshold
        :param l2: Level-2 data container
        :return: None
        """

        # Get options
        filter_options = self.cfg.options.tiepoint_maxdist_filter
        edges_only = filter_options.edges_only
        distance_threshold = filter_options.maximum_distance_to_tiepoint

        # Compute distance to next tie point
        tiepoint_distance = self.get_tiepoint_distance(l2)

        # Get indices
        invalid_indices = np.where(tiepoint_distance > distance_threshold)[0]

        # Only remove sla values at the edges (if edges_only:True)
        if edges_only:
            lead_indices = l2.surface_type.lead.indices
            before_first_lead = invalid_indices < lead_indices[0]
            after_last_lead = invalid_indices > lead_indices[-1]
            edge_condition = np.logical_or(before_first_lead, after_last_lead)
            invalid_indices = invalid_indices[np.where(edge_condition)[0]]

        # Validity check
        if len(invalid_indices) == 0:
            return

        # Set the values
        self.sla[invalid_indices] = np.nan

    def get_tiepoint_distance(self, l2):
        """
        Returns the distance in meter to the next ssh tiepoint for each record
        :param l2: Level-2 data container
        :return: array(float32, shape=l2.n_records)
        """

        # prepare parameter arrays
        lead_indices = l2.surface_type.lead.indices
        lead_elevation = np.full(l2.n_records, np.nan, dtype=np.float32)
        lead_elevation[lead_indices] = self.sla[lead_indices]

        # Compute distance to next lead tiepoint in meter
        tiepoint_distance = get_tiepoint_distance(np.isfinite(lead_elevation))
        tiepoint_distance = tiepoint_distance.astype(np.float32)
        tiepoint_distance *= self.cfg.options.smooth_filter_width_footprint_size
        return tiepoint_distance

    @property
    def filter_width(self):
        """
        Compute the filter width in points
        :return:
        """
        filter_width = self.cfg.options.smooth_filter_width_m / self.cfg.options.smooth_filter_width_footprint_size
        # Make sure filter width is odd integer
        filter_width = np.floor(filter_width) // 2 * 2 + 1
        filter_width = filter_width.astype(int)
        return filter_width

    @property
    def l2_input_vars(self):
        """
        Mandatory property for Level2ProcessorStep children
        :return: list (str)
        """
        return ["surface_type", "elev", "mss"]

    @property
    def l2_output_vars(self):
        """
        Mandatory property for Level2ProcessorStep children
        :return: list (str)
        """
        return ["sla"]

    @property
    def error_bit(self):
        return self.error_flag_bit_dict["sla"]


def get_tiepoint_distance(is_tiepoint):
    """
    Calculates the distance to the next tiepoints (boolean array)
    """
    distance_forward = get_tiepoints_oneway_distance(is_tiepoint)
    distance_reverse = get_tiepoints_oneway_distance(is_tiepoint, reverse=True)
    return np.minimum(distance_forward, distance_reverse)


def get_tiepoints_oneway_distance(a, reverse=False):
    """ loops through array and determines distance to latest flag=true """
    n = len(a)
    distance = np.full(a.shape, n+1, dtype=np.int32)
    if reverse:
        a = a[::-1]
    dist = n
    for i in np.arange(n):
        if a[i]:
            dist = 0
        elif dist == n:
            pass
        else:
            dist += 1
        distance[i] = dist
    if reverse:
        distance = distance[::-1]
    return distance


def gaussian_process(sla_raw):

    from sklearn import gaussian_process
    from sklearn.gaussian_process.kernels import Matern, WhiteKernel

    n = len(sla_raw)
    x = np.arange(n)
    y = np.copy(sla_raw)
    ssh_indices = np.where(np.isfinite(y))[0]
    mean_val = np.nanmean(sla_raw)
    y -= mean_val
    x_fit = x[ssh_indices].reshape(-1, 1)
    y_fit = y[ssh_indices].reshape(-1, 1)
    matern_kernel_props = dict(length_scale=200.0, length_scale_bounds=(10, 1000))
    white_noise_kernel_props = dict(noise_level=1, noise_level_bounds=(0.5, 5))
    kernel = Matern(**matern_kernel_props) + WhiteKernel(**white_noise_kernel_props)
    gp = gaussian_process.GaussianProcessRegressor(kernel=kernel)
    gp.fit(x_fit, y_fit)
    x_pred = x.reshape(-1, 1)
    sla, sigma = gp.predict(x_pred, return_std=True)
    return sla + mean_val
