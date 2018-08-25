# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 13:57:56 2016

@author: Stefan
"""

from pysiral.auxdata import AuxdataBaseClass, GridTrackInterpol
from pysiral.filter import idl_smooth
from pysiral.iotools import ReadNC

import scipy.ndimage as ndimage

from pyproj import Proj
import numpy as np
import os


class SnowBaseClass(AuxdataBaseClass):

    def __init__(self):
        super(SnowBaseClass, self).__init__()

    def get_along_track_snow(self, l2):
        """ Standard API method """
        # TODO: This might be moved into a standard get along-track data method
        snow, msg = self._get_along_track_snow(l2)
        return snow, msg


class NoneHandler(SnowBaseClass):

    def __init__(self):
        super(NoneHandler, self).__init__()

    def _get_along_track_snow(self, l2):
        dummy = np.full((l2.n_records), np.nan)
        snow = SnowParameterContainer()
        snow.depth = dummy
        snow.density = dummy
        snow.depth_uncertainty = dummy
        snow.density_uncertainty = dummy
        return snow, ""


class Warren99(SnowBaseClass):

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
    p = Proj(proj="stere", lat_0=90, lon_0=-90, lat_ts=70)

    def __init__(self):
        super(Warren99, self).__init__()

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

    def _get_along_track_snow(self, l2):
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
        valid_min, valid_max = self._options.valid_snow_depth_range
        invalid = np.logical_or(snow.depth < valid_min, snow.depth > valid_max)
        invalid_records = np.where(invalid)[0]
        snow.set_invalid(invalid_records)

        # Apply ice_type (myi_fraction correction)
        scale_factor = (1.0 - l2.sitype) * self._options.fyi_correction_factor

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
        if self._options.smooth_snow_depth:
            filter_width = self._options.smooth_filter_width_m
            # Convert filter width to index
            filter_width /= l2.footprint_spacing
            # Round to odd number
            filter_width = np.floor(filter_width) // 2 * 2 + 1
            snow.depth = idl_smooth(snow.depth, filter_width)

        return snow, ""

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
        snow_depth = sd[0] + sd[1]*l2x + sd[2]*l2y + \
            sd[3]*l2x*l2y + sd[4]*l2x*l2x + sd[5]*l2y*l2y
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


class Warren99AMSR2Clim(SnowBaseClass):
    """ Class for monthly snow depth & density climatology based on merged Warren99 climatology and
     monthly AMSR2 snow depth composite (source: IUP) """

    def __init__(self):
        super(Warren99AMSR2Clim, self).__init__()
        self._data = None
        self.error.caller_id = self.__class__.__name__

    def _get_along_track_snow(self, l2):
        """ This is the method that will be evoked by the Level-2 processor """

        # Set the requested date
        self.set_requested_date_from_l2(l2)

        # Update the external data
        # NOTE: This will only be done once as the climatology has the same period as the Level-2 processor
        self.update_external_data()

        # Check if error with file I/O
        if self.error.status:
            return None, self._msg

        # Extract along track snow depth and density
        snow = self._get_snow_track(l2)

        return snow, self._msg

    def load_requested_auxdata(self):
        """ Required subclass method: Load the data file necessary to satisfy condition for requested date"""

        # Retrieve the file path for the requested date from a property of the auxdata parent class
        path = self.requested_filepath

        # Validation
        if not os.path.isfile(path):
            self._msg = self.__class__.__name__+": File not found: %s " % path
            self.error.add_error("auxdata_missing_snow", self._msg)
            return

        # Read the data
        self._data = ReadNC(path)

        # This step is important for calculation of image coordinates
        self._msg = self.__class__.__name__+": Loaded snow file: %s" % path

    def _get_snow_track(self, l2):
        """ Get the along-track data from the loaded data"""

        # Extract track data from grid
        griddef = self._options[l2.hemisphere]
        grid_lons, grid_lats = self._data.longitude, self._data.latitude
        grid2track = GridTrackInterpol(l2.track.longitude, l2.track.latitude, grid_lons, grid_lats, griddef)

        # Extract data (Map the extracted tracks directly on the snow parameter container)
        var_map = self.options.variable_map
        snow = SnowParameterContainer()
        for var_name in var_map.keys():
            source_name = var_map[var_name]
            sdgrid = getattr(self._data, source_name)
            setattr(snow, var_name, grid2track.get_from_grid_variable(sdgrid))

        # Extract the W99 weight
        w99_weight = grid2track.get_from_grid_variable(self._data.w99_weight)

        # Apply the same modification as the Warren climatology
        # Apply ice_type (myi_fraction correction) but this time modified by the regional weight
        # of the Warren climatology. The weight ranges from 0 to 1 and make sure no fyi scaling is
        # applied over the AMSR2 region data
        scale_factor = (1.0 - l2.sitype) * self._options.fyi_correction_factor * w99_weight

        # The scaling factor affects the snow depth ...
        snow.depth = snow.depth - scale_factor * snow.depth

        # ... and the uncertainty. Here it is assumed that the uncertainty
        # is similar affected by the scaling factor.
        snow.depth_uncertainty = snow.depth_uncertainty - scale_factor * snow.depth_uncertainty

        # the uncertainty of the myi fraction is acknowledged by adding
        # an additional term that depends on snow depth, the magnitude of
        # scaling and the sea ice type uncertainty
        scaling_uncertainty = snow.depth * scale_factor * l2.sitype.uncertainty * w99_weight
        snow.depth_uncertainty = snow.depth_uncertainty + scaling_uncertainty

        return snow


class FixedSnowDepthDensity(SnowBaseClass):
    """ Always returns zero snow depth """

    def __init__(self):
        super(FixedSnowDepthDensity, self).__init__()

    def _get_along_track_snow(self, l2):
        snow_depth = np.ones(shape=(l2.n_records), dtype=np.float32)
        snow_depth *= self._options.fixed_snow_depth
        snow_density = np.ones(shape=(l2.n_records), dtype=np.float32)
        snow_density *= self._options.fixed_snow_density

        snow = SnowParameterContainer()
        snow.depth = snow_depth
        snow.density = snow_density

        return snow, ""


class ICDCSouthernClimatology(SnowBaseClass):
    """ Class for daily climatology fields from UHH ICDC """

    def __init__(self):
        super(ICDCSouthernClimatology, self).__init__()
        self._data = None
        self._current_date = [0, 0]
        self._requested_date = [-1, -1]
        self.error.caller_id = self.__class__.__name__

    def _get_along_track_snow(self, l2):

        # Extract requested day
        self._get_requested_date(l2)

        # Check if data for day is already loaded
        if self._requested_date != self._current_date:
            self._get_data(l2)
            self._current_date = self._requested_date
        else:
            self._msg = "ICDCSouthernClimatology: Daily grid already present"

        # Check if error with file I/O
        if self.error.status:
            return None, self._msg

        # Extract along track snow depth and density
        sd, sd_unc = self._get_snow_track(l2)

        # Apply along-track smoothing if required
        if self._options.smooth_snowdepth:
            filter_width = self._options.smooth_filter_width_m
            # Convert filter width to index
            filter_width /= l2.footprint_spacing
            # Round to odd number
            filter_width = np.floor(filter_width) // 2 * 2 + 1
            sd = idl_smooth(sd, filter_width)
            sd_unc = idl_smooth(sd_unc, filter_width)

        # Collect Parameters and return
        # (density and density uncertainty fixed from l2 settings)
        snow = SnowParameterContainer()
        snow.depth = sd
        snow.depth_uncertainty = sd_unc
        snow.density = np.full(sd.shape, self._options.snow_density)
        snow.density_uncertainty = np.full(
            sd.shape, self._options.snow_density_uncertainty)

        return snow, self._msg

    def _get_data(self, l2):
        """ Loads file from local repository only if needed """
        path = self._get_local_repository_filename(l2)
        # Validation
        if not os.path.isfile(path):
            self._msg = "ICDCSouthernClimatology: File not found: %s " % path
            self.error.add_error("auxdata_missing_snow", self._msg)
            return
        self._data = ReadNC(path)
        # This step is important for calculation of image coordinates
        # self._data.ice_conc = np.flipud(self._data.ice_conc)
        self._msg = "ICDCSouthernClimatology: Loaded SIC file: %s" % path

    def _get_requested_date(self, l2):
        """ Use first timestamp as reference, date changes are ignored """
        month = l2.track.timestamp[0].month
        day = l2.track.timestamp[0].day
        self._requested_date = [month, day]

    def _get_local_repository_filename(self, l2):
        """ Get the filename (no subfolders as climatology for each day)"""
        path = self._local_repository
        filename = self._filenaming.format(month=self.month, day=self.day)
        path = os.path.join(path, filename)
        return path

    def _get_snow_track(self, l2):

        # TODO: Move functionality to image coordinate class

        # Convert grid/track coordinates to grid projection coordinates
        kwargs = self._options[l2.hemisphere].projection
        p = Proj(**kwargs)
        x, y = p(self._data.lon, self._data.lat)
        l2x, l2y = p(l2.track.longitude, l2.track.latitude)

        # Convert track projection coordinates to image coordinates
        # x: 0 < n_lines; y: 0 < n_cols
        dim = self._options[l2.hemisphere].dimension
        x_min, y_min = np.nanmin(x), np.nanmin(y)
        ix, iy = (l2x-x_min)/dim.dx, (l2y-y_min)/dim.dy

        # Extract snow depth along track data from grid
        sd_parameter_name = self._options.snow_depth_nc_variable
        sdgrid = getattr(self._data, sd_parameter_name)[0, :, :]
        sdgrid = np.flipud(sdgrid)
        # sdgrid = np.roll(sdgrid, 16, axis=0)
        sd = ndimage.map_coordinates(sdgrid, [iy, ix], order=0)
        sd[sd < 0.0] = np.nan

        # Extract snow depth uncertainty
        unc_parameter_name = self._options.snow_depth_uncertainty_nc_variable
        uncgrid = getattr(self._data, unc_parameter_name)[0, :, :]
        uncgrid = np.flipud(uncgrid)
        # uncgrid = np.roll(uncgrid, 16, axis=0)
        unc = ndimage.map_coordinates(uncgrid, [iy, ix], order=0)
        unc[unc < 0.0] = np.nan

#        import matplotlib.pyplot as plt
#        from mpl_toolkits.basemap import Basemap

#        plt.figure("x")
#        plt.imshow(x, origin="lower")
#
#        plt.figure("y")
#        plt.imshow(y, origin="lower")
#
#
#        plt.figure("lat")
#        plt.imshow(self._data.lat, origin="lower")
#
#        plt.figure("lon")
#        plt.imshow(self._data.lon, origin="lower")
#        lat = self._data.lat
#        lat_diff = np.flipud(lat) - lat
#        plt.figure("lat flip test")
#        plt.imshow(lat_diff)
#        plt.colorbar(label="Latitude Diff (flipped - unflipped)")
#        plt.show()

#        width, height = dim.dy*dim.n_cols, dim.dx*dim.n_lines
#        print "width: ", width
#        print "height: ", height

#        plt.figure("map")
#        m = Basemap(projection='stere', lon_0=0, lat_0=-90, lat_ts=-70,
#                    resolution='i', width=width, height=height)
#        m.drawcoastlines(color='#4b4b4d', linewidth=2, zorder=20)
#        px, py = m(l2.longitude, l2.latitude)
#        m.imshow(sdgrid, vmin=0, vmax=0.5, interpolation="none", alpha=0.5)
#        im = m.imshow(np.flipud(self._data.lon), interpolation="none",
#                      alpha=0.4)
#        plt.colorbar(im, label="Longitude")
#        cs = m.contour(self._data.lon, self._data.lat, self._data.lat,
#                       latlon=True, levels=np.arange(-89, 30, 5),
#                       colors="white", linestyles="solid", linewidths=1)
#        plt.clabel(cs, fontsize=10)
#        m.plot(px, py, lw=2, color="red", zorder=100)
#        m.drawparallels(np.arange(-90., 120., 10.), latmax=88)
#        m.drawmeridians(np.arange(0., 420., 15.), latmax=88,
#                        labels=[True, True, True, True])
#
#        image_extent = (np.amin(x), np.amax(x),
#                        np.amin(y), np.amax(y))
#
#        plt.figure("map coordinates")
#        plt.imshow(sdgrid, vmin=0, vmax=0.5, interpolation="none",
#                   extent=image_extent, origin="lower")
#        plt.plot(l2x, l2y, lw=2, color="red", zorder=100)
#
#        image_extent = (0, dim.n_cols,
#                        0, dim.n_lines)
#
#        plt.figure("image coordinates")
#        plt.imshow(sdgrid, vmin=0, vmax=0.5, interpolation="none",
#                   extent=image_extent, origin="lower")
#        plt.plot(ix, iy, lw=2, color="red", zorder=100)
#        plt.scatter(ix[0], iy[0], color="red", zorder=101)
#
#        plt.figure("snow depth")
#        x = np.arange(l2.n_records)
#        plt.plot(x, sd)
#        plt.xlim(0, len(x))
#        plt.ylim(0, 0.6)
#        plt.show()

        return sd, unc

    @property
    def month(self):
        return "%02g" % self._requested_date[0]

    @property
    def day(self):
        return "%02g" % self._requested_date[1]


class SnowParameterContainer(object):

    def __init__(self):
        self.depth = None
        self.density = None
        self.depth_uncertainty = None
        self.density_uncertainty = None

    def set_invalid(self, indices):
        self.depth[indices] = np.nan
        self.density[indices] = np.nan
        self.depth_uncertainty[indices] = np.nan
        self.density_uncertainty[indices] = np.nan

    def set_dummy(self, n_records):
        self.depth = np.full((n_records), np.nan)
        self.density = np.full((n_records), np.nan)
        self.depth_uncertainty = np.full((n_records), np.nan)
        self.density_uncertainty = np.full((n_records), np.nan)


def get_l2_snow_handler(name):
    pyclass = globals().get(name, None)
    if pyclass is not None:
        return pyclass()
    else:
        return pyclass
#    try:
#        return globals()[name]()
#    except:
#        msg = "Unknown snow depth handler: %s" % name
#        error = ErrorStatus(caller_id="get_l2_snow_handler")
#        error.add_error("invalid-snow-handler", msg)
#        error.raise_on_error()
