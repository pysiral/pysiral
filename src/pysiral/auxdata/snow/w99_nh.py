# -*- coding: utf-8 -*-

import numpy as np
from pyproj import Proj

from pysiral.auxdata import AuxdataBaseClass
from pysiral.auxdata.snow._data import SnowParameterContainer
from pysiral.filter import idl_smooth


class Warren99(AuxdataBaseClass):

    # Snow depth Coefficients
    sd_coefs = np.array([
        [28.01, 0.1270, -1.1833, -0.1164, -0.0051, 0.0243, 7.6, -0.06, 0.07, 4.6],
        [30.28, 0.1056, -0.5908, -0.0263, -0.0049, 0.0044, 7.9, -0.06, 0.08, 5.5],
        [33.89, 0.5486, -0.1996, 0.0280, 0.0216, -0.0176, 9.4, -0.04, 0.10, 6.2],
        [36.80, 0.4046, -0.4005, 0.0256, 0.0024, -0.0641, 9.4, -0.09, 0.09, 6.1],
        [36.93, 0.0214, -1.1795, -0.1076, -0.0244, -0.0142, 10.6, -0.21, 0.09, 6.3],
        [36.59, 0.7021, -1.4819, -0.1195, -0.0009, -0.0603, 14.1, -0.16, 0.12, 8.1],
        [11.02, 0.3008, -1.2591, -0.0811, -0.0043, -0.0959, 9.5, 0.02, 0.10, 6.7],
        [4.64, 0.3100, -0.6350, -0.0655, 0.0059, -0.0005, 4.6, -0.01, 0.05, 3.3],
        [15.81, 0.2119, -1.0292, -0.0868, -0.0177, -0.0723, 7.8, -0.03, 0.06, 3.8],
        [22.66, 0.3594, -1.3483, -0.1063, 0.0051, -0.0577, 8.0, -0.08, 0.06, 4.0],
        [25.57, 0.1496, -1.4643, -0.1409, -0.0079, -0.0258, 7.9, -0.05, 0.07, 4.3],
        [26.67, -0.1876, -1.4229, -0.1413, -0.0316, -0.0029, 8.2, -0.06, 0.07, 4.8]])

    swe_coefs = np.array([
        [8.37, -0.0270, -0.3400, -0.0319, -0.0056, -0.0005, 2.5, -0.005, 0.024, 1.6],
        [9.43, 0.0058, -0.1309, 0.0017, -0.0021, -0.0072, 2.6, -0.007, 0.028, 1.8],
        [10.74, 0.1618, 0.0276, 0.0213, 0.0076, -0.0125, 3.1, 0.007, 0.032, 2.1],
        [11.67, 0.0841, -0.1328, 0.0081, -0.0003, -0.0301, 3.2, -0.013, 0.032, 2.1],
        [11.80, -0.0043, -0.4284, -0.0380, -0.0071, -0.0063, 3.5, -0.047, 0.033, 2.2],
        [12.48, 0.2084, -0.5739, -0.0468, -0.0023, -0.0253, 4.9, -0.030, 0.044, 2.9],
        [4.01, 0.0970, -0.4930, -0.0333, -0.0026, -0.0343, 3.5, 0.008, 0.037, 2.4],
        [1.08, 0.0712, -0.1450, -0.0155, 0.0014, -0.0000, 1.1, -0.001, 0.012, 0.8],
        [3.84, 0.0393, -0.2107, -0.0182, -0.0053, -0.0190, 2.0, -0.003, 0.016, 1.0],
        [6.24, 0.1158, -0.2803, -0.0215, 0.0015, -0.0176, 2.3, -0.005, 0.021, 1.4],
        [7.54, 0.0567, -0.3201, -0.0284, -0.0032, -0.0129, 2.4, -0.000, 0.023, 1.5],
        [8.00, -0.0540, -0.3650, -0.0362, -0.0112, -0.0035, 2.5, -0.003, 0.024, 1.5]])

    earth_radius = 6371000.8
    water_density = 1024.0

    def __init__(self, *args, **kwargs):
        super(Warren99, self).__init__(*args, **kwargs)
        self.p = Proj(proj="stere", lat_0=90, lon_0=-90, lat_ts=70)

    def evaluate(self, lons, lats, month_num):
        """ Return the result of the Warren Climatology for a
        given set of lons, lats and the month number (1-12) """

        # Compute coordinates in cartesian reference system of climatology
        x, y = self.p(lons, lats)

        # convert to degrees of arc
        x = x / (self.earth_radius * np.pi / 180.0)
        y = y / (self.earth_radius * np.pi / 180.0)

        # Get W99 snow depth & uncertainty
        sd = self._get_snow_depth(month_num, x, y)

        # Get W99 snow density
        sdens = self._get_snow_density(sd, month_num, x, y)

        # Get the uncertainties
        sd_unc, sdens_unc = self._get_warren_uncertainty(month_num, sd)

        # Put everything in a container
        snow = SnowParameterContainer()
        snow.depth = sd
        snow.density = sdens
        snow.depth_uncertainty = sd_unc
        snow.density_uncertainty = sdens_unc

        return snow

    def get_l2_track_vars(self, l2):
        """ Get the snow depth, density and their uncertainties for the track in the l2 data object
        including the potential modification of the original climatology and filters """

        # Validate hemisphere
        if l2.hemisphere == "south":
            snow = SnowParameterContainer()
            snow.set_dummy(l2.n_records)
            msg = "Warren99 not valid for southern hemisphere, returning 0"
            self.error.add_error("warren99-invalid-hemisphere", msg)
            return snow

        # Get the original warren climatology values
        # NOTE: snow is a class with properties depth, depth_uncertainty, density & density_uncertainty
        snow = self._get_warren99_fit_from_l2(l2)

        # Filter invalid values
        valid_min, valid_max = self.cfg.options.valid_snow_depth_range
        invalid = np.logical_or(snow.depth < valid_min, snow.depth > valid_max)
        invalid_records = np.where(invalid)[0]
        snow.set_invalid(invalid_records)

        # Apply ice_type (myi_fraction correction)
        scale_factor = (1.0 - l2.sitype) * self.cfg.options.fyi_correction_factor

        # The scaling factor affects the snow depth ...
        snow.depth = snow.depth - scale_factor * snow.depth

        # ... and the uncertainty. Here it is assumed that the uncertainty
        # is similar affected by the scaling factor.
        snow.depth_uncertainty = snow.depth_uncertainty - scale_factor * snow.depth_uncertainty

        # the uncertainty of the myi fraction is acknowledged by adding
        # an additional term that depends on snow depth, the magnitude of
        # scaling and the sea ice type uncertainty
        scaling_uncertainty = snow.depth * scale_factor * l2.sitype.uncertainty
        snow.depth_uncertainty = snow.depth_uncertainty + scaling_uncertainty

        # Smooth snow depth (if applicable)
        if self.cfg.options.smooth_snow_depth:
            filter_width = self.cfg.options.smooth_filter_width_m
            # Convert filter width to index
            filter_width /= l2.footprint_spacing
            # Round to odd number
            filter_width = np.floor(filter_width) // 2 * 2 + 1
            snow.depth = idl_smooth(snow.depth, filter_width)

        # Register Variables
        self.register_auxvar("sd", "snow_depth", snow.depth, snow.depth_uncertainty)
        self.register_auxvar("sdens", "snow_density", snow.density, snow.density_uncertainty)

    def _get_warren99_fit_from_l2(self, l2):
        """ This convinience function translates the information from the l2 object
        for the evaluate method """
        # get projection coordinates
        month = l2.track.timestamp[0].month
        snow = self.evaluate(l2.track.longitude, l2.track.latitude, month)
        return snow

    def _get_sd_coefs(self, month):
        return self.sd_coefs[month-1, 0:6]

    def _get_swe_coefs(self, month):
        return self.swe_coefs[month-1, 0:6]

    def _get_snow_depth(self, month, l2x, l2y):
        sd = self._get_sd_coefs(month)
        snow_depth = sd[0] + sd[1]*l2x + sd[2]*l2y + sd[3]*l2x*l2y + sd[4]*l2x*l2x + sd[5]*l2y*l2y
        snow_depth *= 0.01
        return snow_depth

    def _get_warren_uncertainty(self, month, sd):
        """
        Get the uncertainty from the Warren climatology for
        snow depth and density

        snow depth:
            sum of fit rms and interannual variability

        snow density
            fit rms of snow water equivalent
        """
        # get w99 coefficients
        sd_coef = self.sd_coefs[month-1, :]
        swe_coef = self.swe_coefs[month-1, :]

        # Snow depth uncertainties
        sd_rms_fit_error = np.full(sd.shape, sd_coef[6]*0.01)
        sd_interannual_var = np.full(sd.shape, sd_coef[9]*0.01)
        sd_unc = sd_rms_fit_error + sd_interannual_var

        # Snow density uncertainty
        sdens_rms_fit_error = (swe_coef[6]*0.01)/sd*self.water_density
        sdens_interannual_var = (swe_coef[9]*0.01)/sd*self.water_density
        sdens_unc = sdens_rms_fit_error + sdens_interannual_var

        return sd_unc, sdens_unc

    def _get_snow_density(self, snow_depth, month, l2x, l2y):
        """ Extract along-track snow density """

        # get snow water equivalent coefs
        swe = self._get_swe_coefs(month)
        snow_water_equivalent = swe[0] + swe[1]*l2x + swe[2]*l2y + swe[3]*l2x*l2y + swe[4]*l2x*l2x + swe[5]*l2y*l2y
        snow_water_equivalent *= 0.01

        # Convert sd and swe to snow density
        snow_density = snow_water_equivalent/snow_depth*self.water_density

        return snow_density
