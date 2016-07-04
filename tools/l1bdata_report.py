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
                                Paragraph, PageBreak, Image, Spacer, Frame)
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
        self.temp_files = []
        self.figure_xsize = 150

    def open_document(self):
        margin = 10*mm
        self.doc = SimpleDocTemplate(
           self.filename, pagesize=A4, leftMargin=margin,
           rightMargin=margin, topMargin=margin, bottomMargin=margin)
        self.doc.pagesize = portrait(A4)

        self.style_heading = self.style['Heading1']
        self.style_subtitle = self.style['Heading3']

    def create_overview_page(self):

        # Table Title
        self.doc.l1b_filename = self.l1b_filename
        self.elements.append(Paragraph("l1bdata Report", self.style_heading))
        self.elements.append(Paragraph(self.l1b_filename, self.style_subtitle))
        self.elements.append(Spacer(width=0, height=5*mm))

        filename = get_temp_png_filename()
        l1bdata_orbit_plot(self.l1b, filename)
        self.elements.append(Image(filename, width=110*mm, height=110*mm))
        self.elements.append(Spacer(width=0, height=5*mm))
        self.temp_files.append(filename)

        subtitle = Paragraph("Metadata", self.style_subtitle)
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
        print colors.black

        col_widths = (30*mm, 125*mm)
        rowheights = [4*mm] * len(data)
        t = Table(data, hAlign='CENTER', colWidths=col_widths,
                  rowHeights=rowheights, repeatCols=1)
        t.setStyle(table_style)
        self.elements.append(t)
        self.elements.append(PageBreak())

    def create_waveform_page(self):

        # Settings
        xsize = self.figure_xsize

        # Add chapter title
        subtitle = Paragraph("Waveform Parameter", self.style_subtitle)
        self.elements.append(subtitle)

        # Add waveform plot
        filename_waveform = get_temp_png_filename()
        l1bdata_waveform_plot(self.l1b, filename_waveform)
        self.elements.append(Image(filename_waveform, width=xsize*mm,
                                   height=0.5*xsize*mm))

        # Add flag plot
        filename_waveform_flags = get_temp_png_filename()
        l1bdata_waveform_flag_plot(self.l1b, filename_waveform_flags)
        self.elements.append(Image(filename_waveform_flags, width=xsize*mm,
                                   height=xsize*mm))

        # new page
        self.elements.append(PageBreak())

        # Add temp files (figures) to clean up list
        self.temp_files.extend([filename_waveform, filename_waveform_flags])

    def create_range_correction_pages(self):

        # Settings
        xsize = self.figure_xsize

        # Add chapter title
        subtitle = Paragraph("Range Corrections", self.style_subtitle)
        self.elements.append(subtitle)

        # Add range corrections plot
        filename_range_corrections = get_temp_png_filename()
        l1bdata_range_corrections_plot(self.l1b, filename_range_corrections)
        self.elements.append(Image(filename_range_corrections, width=xsize*mm,
                                   height=1.75*xsize*mm))
        # new page
        self.elements.append(PageBreak())

        # Add temp files (figures) to clean up list
        self.temp_files.append(filename_range_corrections)

    def create_classifier_pages(self):

        # Settings
        xsize = self.figure_xsize

        # Add chapter title
        subtitle = Paragraph("L1b Classifier", self.style_subtitle)
        self.elements.append(subtitle)

        # Add range corrections plot
        filename_classifier = get_temp_png_filename()
        l1bdata_classifier_plot(self.l1b, filename_classifier)
        self.elements.append(Image(filename_classifier, width=xsize*mm,
                                   height=1.75*xsize*mm))
        # new page
        self.elements.append(PageBreak())

        # Add temp files (figures) to clean up list
        self.temp_files.append(filename_classifier)

    def close_document(self):
        # self.doc.multiBuild(self.elements, canvasmaker=FooterCanvas)
        self.doc.build(self.elements, onLaterPages=myLaterPages)
        # Clean up temp files
        for temp_file in self.temp_files:
            try:
                os.remove(temp_file)
            except:
                pass


def myLaterPages(canvas, doc):

    canvas.saveState()
    canvas.setFont('Helvetica', 9)
    canvas.drawString(256, 25, "Page %d %s" % (doc.page, doc.l1b_filename))
    canvas.restoreState()


def l1bdata_orbit_plot(l1b, filename):

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

    plt.figure(figsize=(12, 12), facecolor="#4b4b4d")
    m = Basemap(projection='ortho', lon_0=lon_0, lat_0=lat_0, resolution='i')
    m.fillcontinents(color='#4b4b4d', lake_color='#4b4b4d', zorder=20)
    # draw parallels and meridians.
    m.drawparallels(np.arange(-90., 120., 10.), **grid_keyw)
    m.drawmeridians(np.arange(0., 420., 15.), **grid_keyw)

    x, y = m(l1b.time_orbit.longitude, l1b.time_orbit.latitude)
    m.scatter(x, y, s=10, color="#FF1700", edgecolors="none", zorder=202)
    m.scatter(x[0], y[0], marker="D", s=80, edgecolors="#FF1700", zorder=201)
    # m.plot(x, y, color="#003e6e", linewidth=2.0)

    plt.title("")

    m.warpimage(temp_filename, zorder=200, alpha=0.80)
    os.remove(temp_filename)
    plt.tight_layout()
    plt.savefig(filename, dpi=300)


def l1bdata_waveform_plot(l1b, filename):

    plt.figure(facecolor="white", figsize=(12, 6))
    ax = plt.gca()
    ax.set_axis_bgcolor('0.98')
    ax.set_position([0.1, 0.05, 0.7, 0.9])

    image, elevation_range = align_echo_power(
        l1b.waveform.power, l1b.waveform.range, l1b.time_orbit.altitude,
        elevation_limit=[-75, 50])
    image_extent = (0, len(l1b.time_orbit.altitude),
                    np.amin(elevation_range), np.amax(elevation_range))

    im = ax.imshow(
        np.log10(image).transpose(), cmap=plt.get_cmap("magma"),
        interpolation='none', origin='lower', extent=image_extent,
        aspect='auto')

    ax.yaxis.set_tick_params(direction='out')
    ax.yaxis.set_ticks_position('left')
    ax.set_ylabel("Elevation (m)")
    ax.xaxis.set_ticks([])
    ax.spines["left"].set_position(("data", -20))

    spines_to_remove = ["top", "right", "bottom", "left"]
    for spine in spines_to_remove:
        ax.spines[spine].set_color("#4b4b4b")
        ax.spines[spine].set_linewidth(0.25)

    cax = plt.gcf().add_axes([0.85, 0.05, 0.02, 0.9])
    cb = plt.colorbar(im, cax=cax)
    cb.set_label("Echo Power log10(Watt)")
    plt.savefig(filename, dpi=300, facecolor=plt.gcf().get_facecolor())


def l1bdata_waveform_flag_plot(l1b, filename):

    from pysiral.visualization.cmap import (is_valid_cmap, surface_type_cmap,
                                            radar_mode_cmap)
    import matplotlib.gridspec as gridspec

    plt.figure(facecolor="white", figsize=(12, 12))

    axgrid = gridspec.GridSpec(3, 1)
    axgrid.update(left=0.1, right=0.8)

    flags = [l1b.waveform.is_valid,
             l1b.waveform.radar_mode,
             l1b.surface_type.flag]

    cmap_info = [is_valid_cmap(), radar_mode_cmap(), surface_type_cmap()]

    title = ["waveform is_valid flag", "waveform radar mode",
             "l1b surface type flag"]

    for i in np.arange(3):

        # Plot the flag as time line
        ax = plt.subplot(axgrid[i])

        flag = flags[i]
        cmap, labels = cmap_info[i]

        im = plt.imshow([flag], interpolation=None, aspect='auto',
                        vmin=-0.5, vmax=len(labels)-0.5, cmap=cmap,
                        extent=(0, 1, 4, 5))
        ax.set_title(title[i])
        ax.set_ylim([0, 5])
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        spines_to_remove = ["top", "right", "bottom", "left"]
        for spine in spines_to_remove:
            ax.spines[spine].set_visible(False)

        # Plot the colorbar with labels
        position = ax.get_position()
        height = position.y1 - position.y0
        cax = plt.gcf().add_axes([0.85, position.y0, 0.02, height])

        cbar = plt.colorbar(im, orientation="vertical", cax=cax)
        cbar_ticks = np.arange(len(labels))
        cbar.set_ticks(cbar_ticks)
        cbar.set_ticklabels(labels)
        cbar.solids.set_edgecolor("1.0")
        cbar.outline.set_alpha(0.0)
        cbar.ax.tick_params('both', length=0.1, which='major', pad=5)

        # Plot a pie diagram with fractions
        x0 = 0.5 * (position.x1 - position.x0 - 0.8*height) + position.x0
        pax = plt.gcf().add_axes([x0, position.y0, 0.8*height, 0.8*height])

        n = float(len(flag))
        fractions = []
        for i in np.arange(9):
            num = len(np.where(flag == i)[0])
            fractions.append(float(num)/n)

        # Pop fraction = 0
        pie_fractions = np.array(fractions)
        pie_colors = np.array(cmap.colors)

        non_zero = np.nonzero(pie_fractions)[0]
        pie_fractions = pie_fractions[non_zero]
        pie_colors = pie_colors[non_zero]

        wedges, texts, autotexts = pax.pie(
            pie_fractions, colors=pie_colors,
            autopct='%1.1f%%', startangle=90, labeldistance=1.1,
            pctdistance=0.75)

        for wedge in wedges:
            wedge.set_width(0.5)
            wedge.set_ec("none")

    plt.savefig(filename, dpi=300, facecolor=plt.gcf().get_facecolor())


def l1bdata_range_corrections_plot(l1b, filename):

    plt.figure(facecolor="white", figsize=(12, 21))

    from matplotlib.ticker import MultipleLocator

    n = len(l1b.correction.parameter_list)
    f, ax = plt.subplots(n, sharex=True, facecolor="white", figsize=(10, 16))
    for i in np.arange(n):
        ax[i].set_axis_bgcolor('0.98')
        correction, name = l1b.correction.get_parameter_by_index(i)
        ax[i].plot(correction, lw=2, color="#00ace5")
        ax[i].set_title(name)
        ax[i].yaxis.set_minor_locator(MultipleLocator(0.05))
        ax[i].yaxis.grid(True, which='minor')
        ax[i].yaxis.set_tick_params(direction='out')
        ax[i].yaxis.set_ticks_position('left')
        ax[i].xaxis.set_ticks([])

        spines_to_remove = ["top", "right", "bottom"]
        for spine in spines_to_remove:
            ax[i].spines[spine].set_visible(False)
    plt.tight_layout()
    plt.savefig(filename, dpi=300, facecolor=plt.gcf().get_facecolor())


def l1bdata_classifier_plot(l1b, filename):

    plt.figure(facecolor="white", figsize=(12, 21))
    x = np.arange(l1b.n_records)
    n = len(l1b.classifier.parameter_list)
    f, ax = plt.subplots(n, sharex=True, facecolor="white", figsize=(10, 16))
    for i, parameter_name in enumerate(l1b.classifier.parameter_list):
        ax[i].set_axis_bgcolor('0.98')
        classifier = l1b.classifier.get_parameter(parameter_name)
        ax[i].scatter(x, classifier, s=1, edgecolors="none", color="#00ace5")
        ax[i].set_title(parameter_name)
        ax[i].yaxis.set_tick_params(direction='out')
        ax[i].yaxis.set_ticks_position('left')
        ax[i].xaxis.set_ticks([])

        spines_to_remove = ["top", "right", "bottom"]
        for spine in spines_to_remove:
            ax[i].spines[spine].set_visible(False)
    plt.tight_layout()
    plt.savefig(filename, dpi=300, facecolor=plt.gcf().get_facecolor())


def align_echo_power(power, range, altitude, elevation_limit=None):

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

    if elevation_limit is not None:
        in_limit_flag = np.logical_and(elevation_range >= elevation_limit[0],
                                       elevation_range <= elevation_limit[1])
        in_limit = np.where(in_limit_flag)[0]
        elevation_range = elevation_range[in_limit]
        aligned_power = aligned_power[:, in_limit]

    return aligned_power, elevation_range

if __name__ == "__main__":
    l1bdata_report()
