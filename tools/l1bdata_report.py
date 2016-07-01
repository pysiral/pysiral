# -*- coding: utf-8 -*-
"""
Created on Fri Jul 01 14:31:02 2016

@author: shendric
"""

from pysiral.config import ConfigInfo
from pysiral.logging import DefaultLoggingClass
from pysiral.path import (validate_directory, folder_from_filename,
                          filename_from_path, file_basename)
from pysiral.l1bdata import L1bdataNCFile
from pysiral.iotools import get_temp_png_filename

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4, landscape, portrait
from reportlab.platypus import (SimpleDocTemplate, Table, TableStyle,
                                Paragraph, PageBreak, Image, Spacer)
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import mm

import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
import numpy as np

import os
import glob


def l1bdata_report():

    # get the pysiral configuration info
    config = ConfigInfo()

    # parse command line arguments
    jobdef = L1BDataReportSettings(config)
    jobdef.parse_command_line_arguments()
    jobdef.collect_job_parameter()

    # start the job
    job = L1BDataReport(jobdef)
    job.run()


class L1BDataReportSettings(DefaultLoggingClass):

    def __init__(self, config):
        super(L1BDataReportSettings, self).__init__(self.__class__.__name__)
        self.pysiral_config = config
        self.args = None
        self.output_folder = None
        self.l1b_files = []

    def parse_command_line_arguments(self):
        self.args = self.parser.parse_args()
        self._validate_args()

    def collect_job_parameter(self):
        # get list of l1b data files
        self.get_l1b_file_list()
        # Create output folder
        self.create_output_folder()

    def create_output_folder(self):
        validate_directory(self.output_folder)
        self.log.info("Export directory: %s" % self.output_folder)

    def get_l1b_file_list(self):
        self.l1b_files = glob.glob(self.args.l1bdata_search_pattern)

    def _validate_args(self):
        if self.args.output_directory is not None:
            path = self.args.output_directory
        else:
            path = folder_from_filename(self.args.l1bdata_search_pattern)
            path = os.path.join(path, "report")
        self.output_folder = path

    @property
    def n_l1b_files(self):
        return len(self.l1b_files)

    @property
    def parser(self):

        import argparse

        parser = argparse.ArgumentParser()

        parser.add_argument(
            '-o', action='store', dest='output_directory', default=None,
            help='Output directory for pdf report(s)')

        parser.add_argument(
            action='store', dest='l1bdata_search_pattern',
            help='link to l1bdata file(s)')

        return parser


class L1BDataReport(DefaultLoggingClass):

    def __init__(self, config):
        super(L1BDataReport, self).__init__(self.__class__.__name__)
        self.jobdef = config

    def run(self):

        for l1b_file in self.jobdef.l1b_files:

            # Read the l1b data file
            l1b = L1bdataNCFile(l1b_file)
            l1b.parse()

            # Get the output filename
            filename = self.get_report_filename(l1b_file)
            l1b_filename = filename_from_path(l1b_file)

            # Create the pdf document
            pdf = L1BDataReportPDF(filename, l1b)
            pdf.l1b_filename = l1b_filename
            pdf.open_document()
            pdf.create_overview_page()
            pdf.create_waveform_page()
            pdf.create_range_correction_pages()
            pdf.create_classifier_pages()
            pdf.close_document()

    def get_report_filename(self, l1b_file):
        outfile = filename_from_path(l1b_file)
        outfile = file_basename(outfile)
        return os.path.join(self.jobdef.output_folder, outfile+'.pdf')


class L1BDataReportPDF(object):

    def __init__(self, filename, l1b):
        self.doc = None
        self.filename = filename
        self.l1b_filename = None
        self.elements = []
        self.l1b = l1b
        self.style = getSampleStyleSheet()

    def open_document(self):
        self.doc = SimpleDocTemplate(self.filename, pagesize=A4)
        self.doc.pagesize = portrait(A4)

    def create_overview_page(self):

        # Table Title
        style_heading = self.style['Heading1']
        style_subtitle = self.style['Heading4']

        self.elements.append(Paragraph("l1bdata Report", style_heading))
        self.elements.append(Paragraph(self.l1b_filename, style_subtitle))
        self.elements.append(Spacer(width=0, height=5*mm))

        filename = get_temp_png_filename()
        l1bdata_orbit_plot(self.l1b, filename)
        self.elements.append(Image(filename, width=110*mm, height=110*mm))
        self.elements.append(Spacer(width=0, height=5*mm))
        # os.remove(filename)

        subtitle = Paragraph("Metadata", style_subtitle)
        self.elements.append(subtitle)

        # Add table with metadata info
        data = []
        header = ["Key", "Value"]
        data.append(header)

        # Add parameters
        for key in self.l1b.info.attribute_list:
            entry = [key, str(getattr(self.l1b.info, key))]
            data.append(entry)

        # Table Style
        table_style = TableStyle([
            ('FONT', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONT', (0, 1), (-1, -1), 'Helvetica'),
            ('FONTSIZE', (0, 0), (-1, -1), 6),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('INNERGRID', (0, 0), (-1, -1), 0.01, colors.black),
            ('BOX', (0, 0), (-1, -1), 0.5, colors.black)])

        col_widths = (30*mm, 125*mm)
        rowheights = [4*mm] * len(data)
        t = Table(data, hAlign='CENTER', colWidths=col_widths,
                  rowHeights=rowheights, repeatCols=1)
        t.setStyle(table_style)
        self.elements.append(t)
        self.elements.append(PageBreak())

    def create_waveform_page(self):
        filename = get_temp_png_filename()
        l1bdata_waveform_plot(self.l1b, filename)
        xsize = 150
        self.elements.append(Image(filename, width=xsize*mm,
                                   height=xsize/1.2*mm))
        self.elements.append(Spacer(width=0, height=10*mm))

    def create_range_correction_pages(self):
        pass

    def create_classifier_pages(self):
        pass

    def close_document(self):
        self.doc.build(self.elements)


def l1bdata_orbit_plot(l1b, filename):


    # AWI eisblau #00ace5
    # AWI tiefblau #003e6e
    # AWI grau 1 #4b4b4d
    # AWI grau 2 #bcbdbf

#    projection = {
#        "north": {
#            "projection": "ortho",
#            "lon_0": -45,
#            "lat_0": 80,
#            "resolution": "i"},
#        "south": {
#            "projection": "ortho",
#            "lon_0": 0,
#            "lat_0": -70,
#            "resolution": "i"}
#            }

    hemisphere = "north"
    if l1b.time_orbit.latitude[0] < 0:
        hemisphere = "south"

    xsize, ysize = 4000, 2000

    lats = np.linspace(-0.5*np.pi, 0.5*np.pi, num=ysize)

    if hemisphere == "south":
        shading = np.zeros(shape=(ysize, xsize))
        indices = np.where(lats > np.deg2rad(-90))[0]
        colormap_name = "Greys"
    elif hemisphere == "north":
        shading = np.ones(shape=(ysize, xsize)) * 2.0
        indices = np.where(lats < np.deg2rad(90))[0]
        colormap_name = "Greys_r"
    for i in indices:
        shading[i, :] = (np.sin(lats[i])+1.0)**3.


    temp_filename = get_temp_png_filename()
    fig = plt.figure(figsize=(10, 5))
    ax = plt.axes([0, 0, 1, 1])
    plt.axis('off')
    ax.imshow(np.flipud(shading), cmap=plt.get_cmap(colormap_name),
              vmin=0, vmax=1.1*np.amax(shading))
    plt.savefig(temp_filename, dpi=600)
    plt.close(fig)

    grid_keyw = {"dashes": (None, None),
                 "color": "#bcbdbf",
                 "linewidth": 0.5,
                 "latmax": 88, "zorder": 10}

    lat_0 = np.median(l1b.time_orbit.latitude)
    lon_0 = np.median(l1b.time_orbit.longitude)

#    if np.nanmean(l1b.time_orbit.latitude) < 0:
#        lat_0 *= -1.0

    plt.figure("Sentinel3A Orbit Quickview", figsize=(12, 12),
               facecolor="#4b4b4d")
    m = Basemap(projection='ortho', lon_0=lon_0, lat_0=lat_0, resolution='i')
    m.fillcontinents(color='#4b4b4d', lake_color='#4b4b4d', zorder=20)
    # draw parallels and meridians.
    m.drawparallels(np.arange(-90., 120., 10.), **grid_keyw)
    m.drawmeridians(np.arange(0., 420., 15.), **grid_keyw)

    x, y = m(l1b.time_orbit.longitude, l1b.time_orbit.latitude)
    m.scatter(x, y, c=l1b.surface_type.flag, s=10,
              cmap=plt.get_cmap("cool"), edgecolors="none", zorder=202)
    m.scatter(x[0], y[0], marker="D", s=80, edgecolors="#4b4b4d", zorder=201)
    # m.plot(x, y, color="#003e6e", linewidth=2.0)

    plt.title("")

    m.warpimage(temp_filename, zorder=200, alpha=0.80)
    os.remove(temp_filename)
    plt.tight_layout()
    plt.savefig(filename, dpi=300)


def l1bdata_waveform_plot(l1b, filename):

    plt.figure(facecolor="white", figsize=(12, 10))
    ax = plt.gca()

    image, elevation_range = align_echo_power(
        l1b.waveform.power, l1b.waveform.range, l1b.time_orbit.altitude)
    image_extent = (0, len(l1b.time_orbit.altitude),
                    np.amin(elevation_range), np.amax(elevation_range))

    ref_aspect = 0.4
    num_recs, num_range_bins = image.shape
    waveform_aspect = float(num_recs)/float(num_range_bins)
    aspect = waveform_aspect / ref_aspect

    im = ax.imshow(
        np.log10(image).transpose(), cmap=plt.get_cmap("magma"),
        interpolation='none', origin='lower', extent=image_extent,
        aspect=aspect)

    ax.yaxis.set_tick_params(direction='out')
    ax.yaxis.set_ticks_position('left')
    ax.set_ylabel("Elevation (m)")
    ax.xaxis.set_ticks([])
    ax.spines["left"].set_position(("data", -20))

    spines_to_remove = ["top", "right", "bottom"]
    for spine in spines_to_remove:
        ax.spines[spine].set_visible(False)

    cb = plt.colorbar(im, orientation="vertical", shrink=0.5)
    cb.set_label("Echo Power log(dB)")
    plt.tight_layout()
    plt.savefig(filename, dpi=300)


def align_echo_power(power, range, altitude):

    from scipy.interpolate import interp1d

    n_range_bins = len(power[0, :])
    n_records = len(power[:, 0])

    elevation = np.repeat(altitude, n_range_bins)
    elevation = elevation.reshape(n_records, n_range_bins)
    elevation -= range

    range_step = range[0, 1] - range[0, 0]

    elevation_range = np.arange(np.amin(elevation)-range_step*0.5,
                                np.amax(elevation)+range_step*0.5,
                                range_step)

    aligned_power = np.ndarray(shape=(n_records, len(elevation_range)))
    aligned_power *= np.nan

    for i in np.arange(n_records):
        f = interp1d(elevation[i, :].flatten(), power[i, :].flatten(),
                     bounds_error=False, fill_value=np.nan)
        aligned_power[i, :] = f(elevation_range)

    return aligned_power, elevation_range

if __name__ == "__main__":
    l1bdata_report()
