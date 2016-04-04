# -*- coding: utf-8 -*-
"""
Created on Mon Apr 04 15:43:38 2016

Sandbox for a more operational version of the visualization part of pysiral

@author: shendric
"""

from pysiral.path import folder_from_filename
from pysiral.logging import stdout_logger

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

import uuid
import tempfile
import os

from PIL import Image


def test_background_image_generation():

    log = stdout_logger("test-background-image-generation")
    mapbackground_folder = get_map_background_folder()

    # Create a dark themed background with phong shading
    background = BasemapWarpimageBackground()
    background.log = log
    background.color = NSICDBackgroundColorScheme()
    background.shader = SpherePhongShading(60, 60)
    background.folder = mapbackground_folder
    background.filename = "basemap-background-north-dark.png"
    background.export_to_png()

    # Create a light themed background with
    background = BasemapWarpimageBackground()
    background.log = log
    background.color = AWILightBackgroundColorScheme()
    background.shader = LatitudeShading("north")
    background.folder = mapbackground_folder
    background.filename = "basemap-background-north-light.png"
    background.export_to_png()


def test_visualization_arctic():
    pass


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
        background = Image.open(self._temp_plain_background_map)
#        shading = self.shader.get_shade_alpha_layer(background)
#        background.paste(shading, (0, 0), shading)
        filename = os.path.join(self.folder, self.filename)
        background.save(filename, 'PNG', facecolor="#000000")

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

    def __init__(self, lat_center, lon_center):
        super(SpherePhongShading, self).__init__()


class LatitudeShading(BackgroundShader):

    def __init__(self, hemisphere):
        super(LatitudeShading, self).__init__()


def get_map_background_folder():
    return os.path.join(os.path.abspath(
        folder_from_filename(__file__)), "mapbackground")

if __name__ == "__main__":
    test_background_image_generation()
    test_visualization_arctic()
    test_visualization_antarctic()
