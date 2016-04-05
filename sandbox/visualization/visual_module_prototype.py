# -*- coding: utf-8 -*-
"""
Created on Mon Apr 04 15:43:38 2016

Sandbox for a more operational version of the visualization part of pysiral

@author: shendric
"""

from pysiral.path import folder_from_filename
from pysiral.logging import stdout_logger

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from mpl_toolkits.basemap import Basemap, maskoceans

import numpy as np
import uuid
import tempfile
import os

from PIL import Image


def test_background_image_generation():

    log = stdout_logger("test-background-image-generation")
    mapbackground_folder = get_map_background_folder()

#    # Create a dark themed background with phong shading
#    background = BasemapWarpimageBackground()
#    background.log = log
#    background.color = NSICDBackgroundColorScheme()
#    background.shader = SpherePhongShading(60, -45)
#    background.folder = mapbackground_folder
#    background.filename = "basemap-background-north-dark.png"
#    background.export_to_png()

    # Create a light themed background with
    background = BasemapWarpimageBackground()
    background.log = log
    background.color = AWILightBackgroundColorScheme()
    background.shader = SpherePhongShading(75, 0)
    background.shader.phong_wet = PhongSettings(1.0, 0.2, 0.1, 5)
    background.shader.phong_dry = PhongSettings(0.4, 0.4, 0.3, 5)
    background.folder = mapbackground_folder
    background.filename = "basemap-background-north-light.png"
    background.export_to_png()


def test_visualization_arctic():
    lat_0 = 75
    lon_0 = 0
    fontcolor = "#4b4b4b"
    grid_keyw = {"dashes": (None, None), "color": fontcolor,
                 "linewidth": 0.1, "latmax": 84, "zorder": 120}
#    wm = plt.get_current_fig_manager()
#    wm.
    figure = plt.figure("Grid Data Visualization",  figsize=(12, 12),
                        facecolor="#ffffff")
    #    m = Basemap(projection='nsper', lon_0=lon_0, lat_0=lat_0,
#                satellite_height=h*1000., resolution='l')
    m = Basemap(projection='ortho', lon_0=lon_0, lat_0=lat_0,
                resolution='l')
    m.warpimage(r"mapbackground\basemap-background-north-light.png",
                scale=1.0, zorder=100)
    m.drawmapboundary(color="#ffffff", fill_color="none", linewidth=2.0,
                      zorder=300)
#    x, y = m(grid.pcolor.longitude, grid.pcolor.latitude)
#    data = getattr(grid, "sea_ice_thickness")
#    dataMask = np.isnan(data)
#    data = np.ma.masked_where(dataMask, data)
#    m.pcolor(x, y, data, cmap=plt.get_cmap("plasma"), vmin=0, vmax=5,
#             zorder=110)

#    m.drawparallels(np.arange(-90., 120., 15.), **grid_keyw)
#    m.drawmeridians(np.arange(0., 420., 30.), **grid_keyw)

    plt.savefig("temp3.png", dpi=300, facecolor=figure.get_facecolor(),
                bbox_inches="tight")

    figure_buffer = Image.open("temp3.png")
    width_fract, height_fract = 0.60, 0.60
    width = int(width_fract*figure_buffer.size[0])
    height = int(height_fract*figure_buffer.size[1])
    xpad, ypad = 0.40*(1.-width_fract), 0.5*(1.-height_fract)
    ypad -= 0.4*(1.-height_fract)
    xoff = int(xpad*figure_buffer.size[0])
    yoff = int(ypad*figure_buffer.size[1])
    x1, x2 = xoff, xoff+width
    y1, y2 = yoff, yoff+height

    figure_buffer = np.array(figure_buffer)
    cropped_figure = figure_buffer[y1:y2, x1:x2, :]

    import matplotlib.patches as patches

    plt.figure(figsize=(12, 12), facecolor="#ffffff")
    plt.gca().set_position([0, 0, 1, 1])

    ax = plt.gca()
    im = ax.imshow(cropped_figure)
#    patch = patches.Circle((400, 400), radius=400, transform=ax.transData)

    pad = int(0.15*width)
    x, y = pad, pad
    dx, dy = width - 2*pad, height-2*pad

    patch = patches.FancyBboxPatch(
        [x, y], dx, dy,
        boxstyle=patches.BoxStyle("Round", pad=pad), transform=ax.transData)
    im.set_clip_path(patch)

    plt.annotate("CryoSat-2", (0.05, 0.92), xycoords="axes fraction",
                 color=fontcolor, fontsize=32)
    plt.annotate("March 2015", (0.05, 0.88), xycoords="axes fraction",
                 color=fontcolor, fontsize=24)
    plt.axis('off')

    sm = plt.cm.ScalarMappable(cmap=plt.get_cmap("plasma"),
                               norm=plt.Normalize(vmin=0, vmax=5))
    # fake up the array of the scalar mappable. Urgh...
    sm._A = []
    ax = plt.gca()
    cb_ax_kwargs = {'loc': 3,
                    'bbox_to_anchor': (0.05, 0.83, 1, 1),
                    'width': "30%",
                    'height': "2%",
                    'bbox_transform': ax.transAxes,
                    'borderpad': 0}
    ticks = MultipleLocator(1)
    axins = inset_axes(ax, **cb_ax_kwargs)
    cb = plt.colorbar(sm, cax=axins, ticks=ticks, orientation="horizontal")
    cl = plt.getp(cb.ax, 'xmajorticklabels')
    plt.setp(cl, fontsize=22, color=fontcolor)
    cb.set_label("Sea Ice Thickness (m)", fontsize=22, color=fontcolor)
    cb.outline.set_linewidth(0.2)
    cb.outline.set_alpha(0.0)
    for t in cb.ax.get_yticklines():
        t.set_color("1.0")
    cb.ax.tick_params('both', length=0.1, which='major', pad=10)
    plt.sca(ax)
    # Add the plane marker at the last point.
    from matplotlib.offsetbox import OffsetImage, AnnotationBbox
    logo = np.array(Image.open('AWI_Logo_Blau_RGB.png'))
    im = OffsetImage(logo, zoom=0.20, resample=True, alpha=0.75)
    ab = AnnotationBbox(im, (0.95, 0.89), xycoords='axes fraction',
                        frameon=False, box_alignment=(1, 0))
    # Get the axes object from the basemap and add the AnnotationBbox artist
    plt.gca().add_artist(ab)
    plt.show()


def test_visualization_antarctic():
    pass


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
        self._temp_plain_background_map = os.path.join(
            tempfile.tempdir, str(uuid.uuid4())+".png")
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
        self.continent["color"] = "#00ace5"
        self.continent["lake_color"] = "#00ace5"


class BackgroundShader(object):

    def __init__(self):
        pass


class SpherePhongShading(BackgroundShader):

    def __init__(self, lat0, lon0):
        super(SpherePhongShading, self).__init__()
        # Nadir point of assumed illumination point
        self.lat0 = lat0 * -1.0  # Accounts for array orientation
        self.lon0 = lon0
        self.phong_wet = PhongSettings(0.1, 0.5, 0.6, 35)
        self.phong_dry = PhongSettings(0.3, 0.6, 0.15, 5)
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


def get_map_background_folder():
    return os.path.join(os.path.abspath(
        folder_from_filename(__file__)), "mapbackground")


if __name__ == "__main__":
    mpl.rcParams['font.sans-serif'] = "arial"
    for target in ["xtick.color", "ytick.color", "axes.edgecolor",
                   "axes.labelcolor"]:
        mpl.rcParams[target] = "#4b4b4d"
    test_background_image_generation()
    test_visualization_arctic()
    # test_visualization_antarctic()
