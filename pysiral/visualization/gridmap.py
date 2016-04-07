# -*- coding: utf-8 -*-
"""
Created on Thu Apr 07 16:18:34 2016

@author: shendric
"""

from pysiral.iotools import get_temp_png_filename
from pysiral.maptools import get_landcoastlines
from pysiral.visualization.mapstyle import GridMapAWIStyle

import os
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.basemap import Basemap


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
        cmap = data.get_cmap()
        m.pcolor(x, y, data.grid, cmap=plt.get_cmap(cmap.name),
                 vmin=cmap.vmin, vmax=cmap.vmax, zorder=110)
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
        cmap = self.data.get_cmap()
        sm = plt.cm.ScalarMappable(cmap=plt.get_cmap(cmap.name),
                                   norm=plt.Normalize(vmin=cmap.vmin,
                                                      vmax=cmap.vmax))
        sm._A = []
        cb_ax_kwargs = {
            'loc': 3, 'bbox_to_anchor': (0.04, 0.84, 1, 1),
            'width': "30%", 'height': "2%", 'bbox_transform': ax.transAxes,
            'borderpad': 0}
        ticks = MultipleLocator(cmap.step)
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


class GridMapLabels(object):

    def __init__(self):
        self.title = ""
        self.period = ""
        self.copyright = ""
