from pysiral.errorhandler import ErrorStatus
from pysiral.auxdata import AuxdataBaseClass
from pysiral.iotools import ReadNC
from pysiral.filter import (fill_nan, idl_smooth)

from attrdict import AttrDict
import scipy.ndimage as ndimage
import numpy as np


class SSAInterpolator(object):
    """
    Parent class for sea surface anomaly retrieval and interpolation
    """
    def __init__(self, *args, **kwargs):
        self.error = ErrorStatus()

    def set_options(self, **opt_dict):
        self._options = AttrDict(opt_dict)

    def interpolate(self, l2):
        self._value = np.ndarray(shape=(l2.n_records))
        self._uncertainty = np.ndarray(shape=(l2.n_records))
        self._interpolate(l2)

    @property
    def value(self):
        return self._value

    @property
    def uncertainty(self):
        return self._uncertainty


class SSASmoothedLinear(SSAInterpolator):
    """ Default CS2AWI Method """

    def __init__(self, *args, **kwargs):
        super(SSASmoothedLinear, self).__init__()
        self.error.caller_id = self.__class__.__name__

    def _interpolate(self, l2):
        self._linear_smoothed_interpolation_between_tiepoints(l2)
        self._calculate_uncertainty(l2)
        if "marine_segment_filter" in self._options:
            self._marine_segment_filter(l2)
        if "tiepoint_maxdist_filter" in self._options:
            self._tiepoint_maxdist_filter(l2)

    def _linear_smoothed_interpolation_between_tiepoints(self, l2):
        """ Based in cs2awi code from Robert Ricker """

        # Use ocean and lead elevations
        self._ssh_tiepoints = l2.surface_type.lead.indices
        if self._options.use_ocean_wfm:
            self._ssh_tiepoints.append(l2.surface_type.ocean.indices)
            self._ssh_tiepoints = np.sort(self._ssh_tiepoints)

        # Get initial elevation at tiepoint locations
        mss_frb = l2.elev - l2.mss

        # Remove ssh tiepoints from the list if their elevation
        # corrected by the median offset of all tiepoints from the mss
        # exceeds a certain threshold
        if self._options.pre_filtering:

            # Startup
            index_dict = np.arange(l2.surface_type.lead.num)
            threshold = self._options.pre_filter_maximum_mss_median_offset

            # Compute the mean distance to the mss
            tiepoint_mss_distance = mss_frb[self._ssh_tiepoints]
            valid_points = np.where(np.isfinite(tiepoint_mss_distance))[0]
            tiepoint_mss_distance = tiepoint_mss_distance[valid_points]
            median_mss_offset = np.median(tiepoint_mss_distance)

            # Filter points for outliers
            offset = np.abs(mss_frb[self._ssh_tiepoints] - median_mss_offset)
            valid = np.where(offset < threshold)
            self._ssh_tiepoints = self._ssh_tiepoints[index_dict[valid]]

        self.ssa_raw = np.ndarray(shape=(l2.n_records))*np.nan
        self.ssa_raw[self._ssh_tiepoints] = mss_frb[self._ssh_tiepoints]
        non_tiepoints = np.where(np.isnan(self.ssa_raw))

        # Filtered raw values
        # Use custom implementation of IDL SMOOTH:
        # idl_smooth(x, w) equivalent to SMOOTH(x, w, /edge_truncate, /nan)
        ssa_filter1 = idl_smooth(self.ssa_raw, self.filter_width)

        # Leave only the original ssh tie points
        ssa_filter1[non_tiepoints] = np.nan

        # Fill nans with linear interpolation and contant values at borders
        # python: fill_nan(x) = IDL: FILL_NAN(x, /NEIGHBOUR)
        ssa_filter2 = fill_nan(ssa_filter1)

        # Final smoothing
        ssa = idl_smooth(ssa_filter2, self.filter_width)
        self._value = ssa

##        # TODO: Make example plot of individual filter steps
#        import matplotlib.pyplot as plt
#        x = np.arange(l2.n_records)
#
##        plt.figure("land")
##        plt.plot(l2.surface_type.land.flag)
#
#        lead = l2.surface_type.lead
#
#        plt.figure()
#        plt.plot(x, l2.mss, label="mss")
#        plt.plot(x, l2.mss+ssa, label="mss+ssa")
#        plt.scatter(x, l2.elev, marker="+", alpha=0.5, label="elev")
#        plt.scatter(x[lead.indices], l2.elev[lead.indices],
#                    marker="o", alpha=0.5, label="leads")
#        plt.legend()
#
#        plt.figure()
#        plt.scatter(x, l2.surface_type.flag)
#
##        plt.figure(facecolor="white")
##        # plt.scatter(x, ssa_raw, color="black")
##        plt.scatter(x, ssa_filter1, color="red", alpha=0.5)
##        plt.plot(x, ssa_filter2, color="blue", lw=2, alpha=0.5)
##        plt.plot(x, ssa, color="orange", lw=3)
#        plt.show()
#        stop

    def _calculate_uncertainty(self, l2):
        """
        Components that add to sea surface anomaly uncertainty
        - mss uncertainty (if known)
        - uncertainty of lead elevations
        - distance to next lead tiepoint
        """

        # short cuts to options
        max_distance = self._options.uncertainty_tiepoints_distance_max
        ssa_unc_min = self._options.uncertainty_minimum
        ssa_unc_max = self._options.uncertainty_maximum

        # get tiepoint distance
        tiepoint_distance = self.get_tiepoint_distance(l2)

        # Compute the influence of distance to next tiepoints
        # in the range of 0: minimum influence to 1: maximum influence
        # It is assumed that the uncertainty has a quadratic dependance
        # on tiepoint distance
        tiepoint_distance_scalefact = tiepoint_distance / max_distance
        above_distance_limit = np.where(tiepoint_distance_scalefact > 1.)[0]
        tiepoint_distance_scalefact[above_distance_limit] = 1.
        tiepoint_distance_scalefact = tiepoint_distance_scalefact**2.

        # Compute the ssa uncertainty based on a min/max approach
        # scaled by
        ssa_unc_range = ssa_unc_max - ssa_unc_min
        ssa_unc = ssa_unc_min + ssa_unc_range * tiepoint_distance_scalefact

        # Save in class
        self._uncertainty = ssa_unc

    def _marine_segment_filter(self, l2):
        """ Check all sections divided by land masses for reliable
        information content """

        filter_options = self._options.marine_segment_filter
        minimum_lead_number = filter_options.minimum_lead_number
        footprint_size = self._options.smooth_filter_width_footprint_size
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
                self._value[marine_section_indices] = np.nan

            marine_segments.append(marine_segment)

    def _tiepoint_maxdist_filter(self, l2):
        """  A filter that does not removes ssa values which distance to
        the next ssh tiepoint exceeds a defined threshold """

        # Get options
        filter_options = self._options.tiepoint_maxdist_filter
        edges_only = filter_options.edges_only
        distance_threshold = filter_options.maximum_distance_to_tiepoint

        # Compute distance to next tiepoint
        tiepoint_distance = self.get_tiepoint_distance(l2)

        # Get indices
        invalid_indices = np.where(tiepoint_distance > distance_threshold)[0]

        # Only remove ssa values at the edges (if edges_only:True)
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
        self._value[invalid_indices] = np.nan

    def get_tiepoint_distance(self, l2):
        """ Returns the distance in meter to the next ssh tiepoint for
        each record """

        # prepare parameter arrays
        lead_indices = l2.surface_type.lead.indices
        lead_elevation = np.full((l2.n_records), np.nan, dtype=np.float32)
        lead_elevation[lead_indices] = self._value[lead_indices]

        # Compute distance to next lead tiepoint in meter
        tiepoint_distance = get_tiepoint_distance(np.isfinite(lead_elevation))
        tiepoint_distance = tiepoint_distance.astype(np.float32)
        tiepoint_distance *= self._options.smooth_filter_width_footprint_size
        return tiepoint_distance

    @property
    def filter_width(self):
        filter_width = self._options.smooth_filter_width_m / \
            self._options.smooth_filter_width_footprint_size
        # Make sure filter width is odd integer
        filter_width = np.floor(filter_width) // 2 * 2 + 1
        filter_width = filter_width.astype(int)
        return filter_width


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


def get_l2_ssh_class(name):
    pyclass = globals().get(name, None)
    if pyclass is not None:
        return pyclass()
    else:
        return pyclass
