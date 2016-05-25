# -*- coding: utf-8 -*-
"""
Created on Thu Apr 07 16:21:42 2016

@author: shendric
"""

from pysiral.path import get_module_folder

import os


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


class GridMapPaperStyle(GridMapAWIStyle):

    def __init__(self):
        super(GridMapPaperStyle, self).__init__()
        self.figure.set_keyw(figsize=(12, 12), facecolor="#ffffff")
        self.font.set_color("#4b4b4b")
        self.coastlines.is_active = True
        self.coastlines.set_keyw(color="#000000", linewidth=0.5)
        self.mapboundary.set_keyw(color="#000000", zorder=200,
                                  linewidth=0.25, fill_color="#ecedef")
        self.continents.set_keyw(color="#bcbdbf", lake_color="#bcbdbf")
        self.font.annotation["color"] = "#ecedef"
        self.font.annotation["fontproperties"] = get_custom_font(
            fontsize=32, awi_font=False)


class GridMapSICCILightStyle(GridMapAWIStyle):

    def __init__(self):
        super(GridMapSICCILightStyle, self).__init__()
        self.figure.set_keyw(figsize=(12, 12), facecolor="#ffffff")
        self.font.set_color("#4b4b4b")
        self.logo.tag = "sicci"
        self.background.tag = "light"
        self.clip.is_active = True


class GridMapFontProp(object):

    def __init__(self):
        self.color = "#bcbdbf"
        self.label = {
            "color": self.color,
            "fontproperties": get_custom_font(fontsize=22)}
        self.period = {
            "color": self.color,
            "fontproperties": get_custom_font(fontsize=22)}
        self.title = {
            "color": self.color,
            "fontproperties": get_custom_font(fontsize=32)}
        self.annotation = {
            "color": self.color,
            "fontproperties": get_custom_font(fontsize=22)}
        self.copyright = {
            "fontsize": 18, "color": self.color,
            "fontproperties": get_custom_font(fontsize=18)}

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
        folder = get_module_folder(__file__)
        filename = os.path.join(
            folder, "mapbackground",
            "basemap-background-{hemisphere}-{tag}.png".format(
                hemisphere=self.hemispere, tag=self.tag))
        return filename


class GridMapLogo(object):

    def __init__(self):
        self.is_active = True
        self.tag = "awi-white"
        self.keyw = {"zoom": 0.8, "resample": True, "alpha": 0.75}

    def get_filename(self):
        folder = get_module_folder(__file__)
        filename = os.path.join(
            folder, "logo", "logo-{tag}.png".format(tag=self.tag))
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


def get_custom_font(fontsize=20, awi_font=True):
    import matplotlib.font_manager as fm
    # See if AWI font does exist
    folder = get_module_folder(__file__)
    font_path = os.path.join(folder, "font", "NeoSansW1G-Regular.otf")
    # Fall back to OpenSans if AWI font not installed
    if not os.path.isfile(font_path) or not awi_font:
        font_path = os.path.join(folder, "font", "OpenSans-Regular.ttf")
    prop = fm.FontProperties(fname=font_path, size=fontsize)
    return prop
