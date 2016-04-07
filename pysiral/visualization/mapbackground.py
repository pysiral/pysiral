# -*- coding: utf-8 -*-
"""
Created on Thu Apr 07 16:18:02 2016

@author: shendric
"""

from pysiral.iotools import get_temp_png_filename

import os
import numpy as np
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, maskoceans
from PIL import Image


class BasemapWarpimageBackground(object):

    def __init__(self):
        self.basemap_projection_dict = {
            "projection": "cyl", "llcrnrlat": -90, "urcrnrlat": 90,
            "llcrnrlon": -180, "urcrnrlon": 180, "resolution": "i"}
        self.color = None
        self.shader = None
        self.folder = None
        self.filename = None
        self.log = None
        self.figsize = (12, 6)
        self.dpi = 1200
        self._temp_plain_background_map = None

    def export_to_png(self):
        # Switch interactive plotting off
        # plt.ioff()
        self._get_temp_file()
        self._create_plain_background_map()
        self._shade_plain_background_map()
        self._clean_up_temp_file()

    def _get_temp_file(self):
        self._temp_plain_background_map = get_temp_png_filename()
        self.log.info("Create temp file: %s" % self._temp_plain_background_map)

    def _create_plain_background_map(self):
        fig = plt.figure(figsize=self.figsize)
        plt.gca().set_position([0, 0, 1, 1])
        m = Basemap(**self.basemap_projection_dict)
        m.drawmapboundary(**self.color.mapboundary)
        m.fillcontinents(**self.color.continent)
        plt.savefig(self._temp_plain_background_map, dpi=self.dpi)
        plt.close(fig)

    def _shade_plain_background_map(self):
        # Get the background as rgb [0-1]
        background_im = Image.open(self._temp_plain_background_map)
        background_rgb = np.array(background_im).astype(float) / 256
        # Get the shading.
        shading_rgba = self.shader.get_shading_layer(background_im)
        # Overlay background rgb with shadind rgba
        shaded_background_rgb = np.copy(background_rgb)
        for i in np.arange(3):
            shaded_background_rgb[:, :, i] = shading_rgba[:, :, i] + \
                background_rgb[:, :, i] * (1.0 - shading_rgba[:, :, 3])
        # Write image to file
        filename = os.path.join(self.folder, self.filename)
        mpimg.imsave(filename, shaded_background_rgb, format="png")

    def _clean_up_temp_file(self):
        os.remove(self._temp_plain_background_map)


class BackgroundColorScheme(object):

    def __init__(self):
        self.mapboundary = {"color": "none", "fill_color": "#003e6e"}
        self.continent = {"color": "#4b4b4d", "lake_color": "#4b4b4d"}


class NSICDBackgroundColorScheme(BackgroundColorScheme):

    def __init__(self):
        super(NSICDBackgroundColorScheme, self).__init__()


class AWILightBackgroundColorScheme(BackgroundColorScheme):

    def __init__(self):
        super(AWILightBackgroundColorScheme, self).__init__()
        self.mapboundary["fill_color"] = "#ffffff"
        self.continent["color"] = "#bcbdbf"
        self.continent["lake_color"] = "#bcbdbf"


class BackgroundShader(object):

    def __init__(self):
        pass


class SpherePhongShading(BackgroundShader):

    def __init__(self, lat0, lon0):
        super(SpherePhongShading, self).__init__()
        # Nadir point of assumed illumination point
        self.lat0 = lat0 * -1.0  # Accounts for array orientation
        self.lon0 = lon0
        self.phong_wet = PhongSettings(0.0, 0.5, 0.5, 20)
        self.phong_dry = PhongSettings(0.1, 0.6, 0.4, 5)
        self.shader = ShaderSettings([0, 0, 0, 0])
        self.imsize = None
        self.landmask_resolution = "i"

    def get_shading_layer(self, background):
        self.imsize = background.size
        import timeit
        start_time = timeit.default_timer()
        shading_array = self._get_phong_rgba()
        elapsed = timeit.default_timer() - start_time
        print "Elapsed Time: ", elapsed
        return shading_array

    def _get_phong_rgba(self):
        # First get pixel coordinated (for angle & landmask)
        self._get_pixel_coordinates()
        # Calculate angle between surface point in assumed
        angle = self._get_surface_viewing_angle()
        # Get intensity from angle and simple phong model
        intensity = self._get_phong_intensity(angle, "phong_wet")
        intensity_dry = self._get_phong_intensity(angle, "phong_dry")
        # Merge intensity with land mask
        is_land = self._get_landmask()
        intensity[is_land] = intensity_dry[is_land]
        # Make sure intensity <= 1
        is_too_high = np.where(intensity > 1.)
        intensity[is_too_high] = 0.999
        is_too_low = np.where(intensity < 0.)
        intensity[is_too_low] = 0.
        # Create an image object [0-1] with intensity as alpha channel
        template = Image.new("RGBA", self.imsize, self.shader.start_rgba)
        rgba = np.array(template).astype(float) / 256.
        rgba[:, :, 3] = 1.0 - intensity
        return rgba

    def _get_pixel_coordinates(self):
        lon = np.deg2rad(np.linspace(-180, 180, self.imsize[0]))
        lat = np.deg2rad(np.linspace(-90, 90, self.imsize[1]))
        self.lons, self.lats = np.meshgrid(lon, lat)

    def _get_landmask(self):
        lons, lats = self.lons, self.lats
        data = np.ones(shape=lons.shape, dtype=bool)
        out = maskoceans(np.rad2deg(lons), np.rad2deg(lats),
                         data, resolution=self.landmask_resolution,
                         inlands=False, grid=1.25)
        is_land = np.logical_not(np.flipud(out.mask))
        return np.where(is_land)

    def _get_surface_viewing_angle(self):
        """ Central Angel Calculation in 3D polar coordinates """
        lons, lats = self.lons, self.lats
        # Convert illumination nadir latitude/longitude to radians
        lat0, lon0 = np.deg2rad(self.lat0), np.deg2rad(self.lon0)
        # Calculate central angle between surface point and illumination nadir
        angle = 2.0 * np.arcsin(np.sqrt(np.sin(0.5*(lats-lat0))**2.0 +
                                        np.cos(lats) * np.cos(lat0) *
                                        np.sin(0.5*(lons-lon0))**2.0))
        return angle

    def _get_phong_intensity(self, angle, phong_type):
        # Mask out the "dark" hemisphere
        dark = np.where(angle > np.pi/2.0)
        angle[dark] = np.pi/2.0
        # Get the phong type
        phong = getattr(self, phong_type)
        intensity = phong.ambient + phong.diffuse*np.cos(angle) + \
            phong.specular*np.cos(angle)**phong.roughness
        return intensity


class LatitudeShading(BackgroundShader):

    def __init__(self, hemisphere):
        super(LatitudeShading, self).__init__()


class PhongSettings(object):

    def __init__(self, k_a, k_d, k_s, r):
        self.ambient = k_a
        self.diffuse = k_d
        self.specular = k_s
        self.roughness = r


class ShaderSettings(object):

    def __init__(self, rgba):
        self.start_rgba = tuple(rgba)

    def set_color(self, rgba):
        self.start_rgba = tuple(rgba)