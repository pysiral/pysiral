# -*- coding: utf-8 -*-
"""
Created on Mon Apr 04 15:43:38 2016

Sandbox for a more operational version of the visualization part of pysiral

@author: shendric
"""

from pysiral.path import folder_from_filename
from pysiral.logging import stdout_logger

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
    background.shader = SpherePhongShading(90, 0)
    background.shader.phong_wet = PhongSettings(1.0, 0.0, 0.0, 5)
    background.shader.phong_dry = PhongSettings(0.3, 0.7, 0.4, 10)
    background.folder = mapbackground_folder
    background.filename = "basemap-background-north-light.png"
    background.export_to_png()


def test_visualization_arctic():

    # projection needed for pcolor grid calculation
    # needs to be saved in NC file in the future
    projection = {
        "proj": "laea",
        "lat_0": -90.0,
        "lon_0": 0.0,
        "ellps": "WGS84",
        "datum": "WGS84",
        "units": "m"}
    # Quick and dirty read the netcdf file
    filename = os.path.join(
        r"testdata",
        r"cs2_smos_ice_thickness_20160208_20160214.nc")
    ncdata = NCMaskedGridData(filename)

    # Extract parameter and merge with metadata
    data = GridMapParameter()
    data.set_grid(ncdata.longitude, ncdata.latitude)
    data.set_parameter(ncdata.analysis_thickness, "sea_ice_thickness")
    data.set_projection(**projection)

    # Create the map
#    output = os.path.join(
#        r"testdata",
#        r"cs2smos_20160208_20160214_analysis_thickness_dark.png")
#    gridmap = ArcticGridPresentationMap()
#    gridmap.style = GridMapAWIDarkStyle()
#    gridmap.data = data
#    gridmap.label.title = "CryoSat-2 - SMOS (OI Analysis)"
#    gridmap.label.period = "Feb. 08 till Feb. 14, 2016"
#    gridmap.save2png(output)

    # Light style map
    output = os.path.join(
        r"testdata",
        r"cs2smos_20160208_20160214_analysis_thickness_light.png")
    gridmap = ArcticGridPresentationMap()
    gridmap.style = GridMapAWILightStyle()
    gridmap.data = data
    gridmap.label.title = "CryoSat-2 - SMOS (OI Analysis)"
    gridmap.label.period = "Feb. 08 till Feb. 14, 2016"
    gridmap.save2png(output)


def test_visualization_antarctic():
    pass


class NCMaskedGridData(object):

    def __init__(self, filename):
        self.filename = filename
        self.parse()

    def parse(self):
        from pysiral.iotools import ReadNC
        nc = ReadNC(self.filename)
        for parameter in nc.parameters:
            data = np.ma.array(getattr(nc, parameter))
            data.mask = np.isnan(data)
            setattr(self, parameter, data)


class GridMapParameter():
    """
    Contains data, pcolor grid calculation capability, colormap definition
    and standardized parameter naming
    """

    def __init__(self):
        from pysiral.config import get_parameter_definitions
        self._parameter_definitions = get_parameter_definitions()
        self._projection = None
        self.latitude = None
        self.longitude = None
        self.grid = None
        self.pgrid = None
        self.pardef = None

    def get_label(self):
        return self.pardef.label+" ("+self.pardef.unit+")"

    def set_grid(self, longitude, latitude):
        self.longitude = longitude
        self.latitude = latitude

    def set_parameter(self, grid, parameter_name):
        self.grid = grid
        self.pardef = self._parameter_definitions[parameter_name]

    def set_projection(self, **projection):
        from pysiral.maptools import GeoPcolorGrid
        self.pgrid = GeoPcolorGrid(self.longitude, self.latitude)
        self.pgrid.calc_from_proj(**projection)


class ArcticGridPresentationMap(object):

    def __init__(self):
        self.output = None
        self.temp_file = get_temp_png_filename()
        self.data = None
        self.style = GridMapAWIStyle()
        self.label = GridMapLabels()

    @property
    def projection(self):
        return {"projection": "ortho", "lon_0": 0, "lat_0": 75,
                "resolution": "l"}

    def save2png(self, output):
        self.output = output
        # 1. Create Orthographic map with data plot
        self._create_orthographic_map()
        # 2. Crop and clip Orthographic map, add labels, save
        self._crop_orthographic_map()
        # Remove tempory fils
        self._clean_up()

    def _create_orthographic_map(self):
        # switch off interactive plotting
        plt.ioff()
        figure = plt.figure(**self.style.figure.keyw)
        m = Basemap(**self.projection)
        # load the (shaded) background
        filename = self.style.background.get_filename("north")
        m.warpimage(filename, **self.style.background.keyw)
        # coastline
        if self.style.coastlines.is_active:
            coastlines = get_landcoastlines(m, **self.style.coastlines.keyw)
            plt.gca().add_collection(coastlines)
        # Plot the data as pcolor grid
        data = self.data
        x, y = m(data.pgrid.longitude, data.pgrid.latitude)
        cmap = plt.get_cmap(data.pardef.cmap.name)
        vmin, vmax = data.pardef.cmap.vmin, data.pardef.cmap.vmax
        m.pcolor(x, y, data.grid, cmap=cmap, vmin=vmin, vmax=vmax, zorder=110)
        # Draw the grid
        # XXX: Skip for noe
        plt.savefig(self.temp_file, dpi=600, facecolor=figure.get_facecolor(),
                    bbox_inches="tight")
        plt.close(figure)

    def _crop_orthographic_map(self):
        # Read the temporary files
        # (only this works, else the image size is capped by screen resolution
        # TODO: try with plt.ioff()
        image = Image.open(self.temp_file)
        imarr = np.array(image)
        # crop the full orthographic image, do not change projections
        x1, x2, y1, y2 = self.style.crop.get_crop_region(image.size)
        cropped_image = imarr[y1:y2, x1:x2, :]
        # Create a new figure
        figure = plt.figure(**self.style.figure.keyw)
        ax = plt.gca()
        # Display cropped image
        plt.axis('off')
        im = ax.imshow(cropped_image)
        ax.set_position([0, 0, 1, 1])
        # clip the image
        if self.style.clip.is_active:
            patch = self.style.clip.get_patch(x1, x2, y1, y2, ax)
            im.set_clip_path(patch)
        # Add labels
        plt.annotate(self.label.title, (0.04, 0.93), xycoords="axes fraction",
                     **self.style.font.title)
        plt.annotate(self.label.period, (0.04, 0.89), xycoords="axes fraction",
                     **self.style.font.period)
        # Add colorbar
        vmin, vmax = self.data.pardef.cmap.vmin, self.data.pardef.cmap.vmax
        sm = plt.cm.ScalarMappable(cmap=plt.get_cmap("plasma"),
                                   norm=plt.Normalize(vmin=vmin, vmax=vmax))
        sm._A = []
        cb_ax_kwargs = {
            'loc': 3, 'bbox_to_anchor': (0.04, 0.84, 1, 1),
            'width': "30%", 'height': "2%", 'bbox_transform': ax.transAxes,
            'borderpad': 0}
        ticks = MultipleLocator(1)
        axins = inset_axes(ax, **cb_ax_kwargs)
        cb = plt.colorbar(sm, cax=axins, ticks=ticks, orientation="horizontal")
        cl = plt.getp(cb.ax, 'xmajorticklabels')
        plt.setp(cl, **self.style.font.label)
        parameter_label = self.data.get_label()
        cb.set_label(parameter_label, **self.style.font.label)
        cb.outline.set_linewidth(0.2)
        cb.outline.set_alpha(0.0)
        for t in cb.ax.get_yticklines():
            t.set_color("1.0")
        cb.ax.tick_params('both', length=0.1, which='major', pad=10)
        plt.sca(ax)

        # Add the plane marker at the last point.
        from matplotlib.offsetbox import OffsetImage, AnnotationBbox
        logo_filename = self.style.logo.get_filename()
        logo = np.array(Image.open(logo_filename))
        im = OffsetImage(logo, **self.style.logo.keyw)
        ab = AnnotationBbox(im, (0.95, 0.89), xycoords='axes fraction',
                            frameon=False, box_alignment=(1, 0))
        ax.add_artist(ab)
        # Now save the map
        plt.savefig(self.output, dpi=self.style.dpi,
                    facecolor=figure.get_facecolor())
        plt.close(figure)

    def _clean_up(self):
        os.remove(self.temp_file)


class GridMapAWIStyle(object):

    def __init__(self):
        self.dpi = 300
        self.figure = BasemapStyleDef()
        self.figure.set_keyw(figsize=(12, 12), facecolor="#ffffff")
        self.mapboundary = BasemapStyleDef()
        self.coastlines = BasemapStyleDef()
        self.grid = BasemapGridDef()
        self.crop = BasemapOrthoCrop()
        self.clip = BasemapImageClip()
        self.background = GridMapBackground()
        self.logo = GridMapLogo()
        self.font = GridMapFontProp()


class GridMapAWIDarkStyle(GridMapAWIStyle):

    def __init__(self):
        super(GridMapAWIDarkStyle, self).__init__()
        self.figure.set_keyw(figsize=(12, 12), facecolor="#000000")
        self.coastlines.is_active = True
        self.coastlines.set_keyw(color="#bcbdbf", linewidth=0.1)


class GridMapAWILightStyle(GridMapAWIStyle):

    def __init__(self):
        super(GridMapAWILightStyle, self).__init__()
        self.figure.set_keyw(figsize=(12, 12), facecolor="#ffffff")
        self.font.set_color("#4b4b4b")
        self.logo.tag = "awi-blue"
        self.background.tag = "light"
        self.clip.is_active = True


class GridMapLabels(object):

    def __init__(self):
        self.title = ""
        self.period = ""
        self.copyright = ""


class GridMapFontProp(object):

    def __init__(self):
        self.color = "#bcbdbf"
        self.label = {
            "color": self.color,
            "fontproperties": self.get_custom_font(fontsize=22)}
        self.period = {
            "color": self.color,
            "fontproperties": self.get_custom_font(fontsize=22)}
        self.title = {
            "color": self.color,
            "fontproperties": self.get_custom_font(fontsize=32)}
        self.copyright = {
            "fontsize": 18, "color": self.color,
            "fontproperties": self.get_custom_font(fontsize=18)}

    def get_custom_font(self, fontsize=20):
        import matplotlib.font_manager as fm
        # See if AWI font does exist
        font_path = os.path.join("font", "NeoSansW1G-Regular.otf")
        # Fall back to OpenSans if AWI font not installed
        if not os.path.isfile(font_path):
            font_path = os.path.join("font", "OpenSans-Regular.ttf")
        prop = fm.FontProperties(fname=font_path, size=fontsize)
        return prop

    def set_color(self, color):
        self.color = color
        for target in ["label", "period", "title", "copyright"]:
            prop = getattr(self, target)
            prop["color"] = color
            setattr(self, target, prop)


class GridMapBackground(object):

    def __init__(self):
        self.hemispere = "north"
        self.tag = "dark"
        self.keyw = {"scale": 1.0, "zorder": 100}

    def get_filename(self, hemisphere):
        filename = os.path.join(
            "mapbackground",
            "basemap-background-{hemisphere}-{tag}.png".format(
                hemisphere=self.hemispere, tag=self.tag))
        return filename


class GridMapLogo(object):

    def __init__(self):
        self.is_active = True
        self.tag = "awi-white"
        self.keyw = {"zoom": 0.8, "resample": True, "alpha": 0.75}

    def get_filename(self):
        filename = os.path.join(
            "logo", "logo-{tag}.png".format(tag=self.tag))
        return filename


class BasemapStyleDef(object):

    def __init__(self):
        self.is_active = False
        self.keyw = {}

    def set_keyw(self, **keyw):
        self.keyw = keyw
        for key in keyw.keys():
            setattr(self, key, keyw[key])


class BasemapGridDef(object):

    def __init__(self):
        self.is_active = True
        self.parallels = []
        self.meridians = []
        self.keyw = []

    def add_grid(self, parallels, meridians, **keyw):
        self.parallels.append(parallels)
        self.meridians.append(meridians)
        self.keyw.append(keyw)

    def get_grids(self):
        return [self.parallels, self.meridians, self.keyw]


class BasemapOrthoCrop(object):

    def __init__(self):
        self.width_fract = 0.60
        self.height_fract = 0.60
        self.xpad_fact = 0.40
        self.ypad_fact = 0.50
        self.size = (0, 0)

    def get_crop_region(self, orig_size):
        self.size = orig_size
        self.width = int(self.width_fract*self.size[0])
        self.height = int(self.height_fract*self.size[1])
        xpad, ypad = 0.40*(1.-self.width_fract), 0.5*(1.-self.height_fract)
        ypad -= 0.4*(1.-self.height_fract)
        xoff = int(xpad*self.size[0])
        yoff = int(ypad*self.size[1])
        x1, x2 = xoff, xoff+self.width
        y1, y2 = yoff, yoff+self.height
        return x1, x2, y1, y2


class BasemapImageClip(object):

    def __init__(self):
        self.is_active = False

    def get_patch(self, x1, x2, y1, y2, ax):
        import matplotlib.patches as patches
        width, height = x2-x1, y2-y1
        patch = patches.Circle((width/2, height/2), radius=width/2,
                               transform=ax.transData)
        return patch


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
#        self.continent["color"] = "#00ace5"
#        self.continent["lake_color"] = "#00ace5"


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


def get_temp_png_filename():
    return os.path.join(tempfile.tempdir, str(uuid.uuid4())+".png")


def get_map_background_folder():
    return os.path.join(os.path.abspath(
        folder_from_filename(__file__)), "mapbackground")


def get_landcoastlines(basemap, color="0.0", linewidth=1):
    from matplotlib.collections import LineCollection
    landpolygons = np.where(np.array(basemap.coastpolygontypes) == 1)[0]
    landsegs = []
    for index in landpolygons:
        landsegs.append(basemap.coastsegs[index])
    landcoastlines = LineCollection(landsegs, antialiaseds=(1,), zorder=120)
    landcoastlines.set_color(color)
    landcoastlines.set_linewidth(linewidth)
    return landcoastlines

if __name__ == "__main__":
    test_background_image_generation()
    test_visualization_arctic()
    # test_visualization_antarctic()
