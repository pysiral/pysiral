# -*- coding: utf-8 -*-
"""
Created on Fri Jul 01 14:31:02 2016

@author: shendric
"""

from pysiral.config import ConfigInfo
from pysiral.logging import DefaultLoggingClass
from pysiral.path import (validate_directory, folder_from_filename,
                          filename_from_path, file_basename)
from pysiral.l2data import L2iNCFileImport
from pysiral.iotools import get_temp_png_filename
from pysiral.visualization.l2ivis import (PlotL2IdataOrbit)

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4, portrait
from reportlab.platypus import (SimpleDocTemplate, Table, TableStyle,
                                Paragraph, PageBreak, Image, Spacer)
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import mm


import numpy as np

import os
import glob


def l2idata_report():

    # get the pysiral configuration info
    config = ConfigInfo()

    # parse command line arguments
    jobdef = L2IDataReportSettings(config)
    jobdef.parse_command_line_arguments()
    jobdef.collect_job_parameter()

    # start the job
    job = L2IDataReport(jobdef)
    job.run()


class L2IDataReportSettings(DefaultLoggingClass):

    def __init__(self, config):
        super(L2IDataReportSettings, self).__init__(self.__class__.__name__)
        self.pysiral_config = config
        self.args = None
        self.output_folder = None
        self.l2i_files = []

    def parse_command_line_arguments(self):
        self.args = self.parser.parse_args()
        self._validate_args()

    def collect_job_parameter(self):
        # get list of l1b data files
        self.get_l2i_file_list()
        # Create output folder
        self.create_output_folder()

    def create_output_folder(self):
        validate_directory(self.output_folder)
        self.log.info("Export directory: %s" % self.output_folder)

    def get_l2i_file_list(self):
        self.l2i_files = glob.glob(self.args.l2idata_search_pattern)

    def _validate_args(self):
        if self.args.output_directory is not None:
            path = self.args.output_directory
        else:
            path = folder_from_filename(self.args.l2idata_search_pattern)
            path = os.path.join(path, "report")
        self.output_folder = path

    @property
    def n_l2i_files(self):
        return len(self.l2i_files)

    @property
    def parser(self):

        import argparse

        parser = argparse.ArgumentParser()

        parser.add_argument(
            '-o', action='store', dest='output_directory', default=None,
            help='Output directory for pdf report(s)')

        parser.add_argument(
            action='store', dest='l2idata_search_pattern',
            help='link to l2i file(s)')

        return parser


class L2IDataReport(DefaultLoggingClass):

    def __init__(self, config):
        super(L2IDataReport, self).__init__(self.__class__.__name__)
        self.jobdef = config

    def run(self):

        for l2i_file in self.jobdef.l2i_files:

            # Read the l1b data file
            l2i = L2iNCFileImport(l2i_file)

            # Get the output filename
            filename = self.get_report_filename(l2i_file)
            l2i_filename = filename_from_path(l2i_file)

            # Create the pdf document
            pdf = L2IDataReportPDF(filename, l2i)
            pdf.l2i_filename = l2i_filename
            pdf.open_document()
            pdf.create_overview_page()
            pdf.create_auxilary_data_page()
            pdf.create_surface_type_page()
            pdf.create_elevation_page()
            pdf.create_freeboard_thickness_page()
            pdf.close_document()

    def get_report_filename(self, l2i_file):
        outfile = filename_from_path(l2i_file)
        outfile = file_basename(outfile)
        return os.path.join(self.jobdef.output_folder, outfile+'.pdf')


class L2IDataReportPDF(object):

    def __init__(self, filename, l2i):
        self.doc = None
        self.filename = filename
        self.l2i_filename = None
        self.elements = []
        self.l2i = l2i
        self.style = getSampleStyleSheet()
        self.temp_files = []
        self.figure_xsize = 175

    def open_document(self):
        top_bottom_margin = 20*mm
        left_right_margin = 10*mm
        self.doc = SimpleDocTemplate(
           self.filename, pagesize=A4,
           leftMargin=left_right_margin, rightMargin=left_right_margin,
           topMargin=top_bottom_margin, bottomMargin=top_bottom_margin)
        self.doc.pagesize = portrait(A4)

        self.style_heading = self.style['Heading1']
        self.style_subtitle = self.style['Heading3']

    def create_overview_page(self):

        xsize = 100

        # Table Title
        self.doc.l2i_filename = self.l2i_filename
        self.elements.append(Paragraph("l2i Orbitdata Report",
                                       self.style_heading))
        self.elements.append(Paragraph(self.l2i_filename, self.style_subtitle))
        self.elements.append(Spacer(width=0, height=5*mm))

        filename = get_temp_png_filename()
        plot = PlotL2IdataOrbit()
        plot.create_plot(self.l2i)
        plot.savefig(filename)

        width = xsize*mm
        height = plot.canvas_aspect*xsize*mm
        self.elements.append(Image(filename, width=width, height=height))
        self.temp_files.append(filename)

        subtitle = Paragraph("Metadata", self.style_subtitle)
        self.elements.append(subtitle)

        # Add table with metadata info
        data = []
        header = ["Key", "Value"]
        data.append(header)

        # Add parameters
        for key in self.l2i.attribute_list:
            entry = [key, str(getattr(self.l2i.info, key))]
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

    def create_auxilary_data_page(self):

        # Settings
        xsize = self.figure_xsize

        # Add chapter title
        subtitle = Paragraph("Auxiliary Parameter", self.style_subtitle)
        self.elements.append(subtitle)
        self.elements.append(Spacer(width=0, height=10*mm))

#        # Add waveform plot
#        filename = get_temp_png_filename()
#        plot = PlotL1bdataWaveform()
#        plot.create_plot(self.l1b)
#        plot.savefig(filename)
#
#        width = xsize*mm
#        height = plot.canvas_aspect*xsize*mm
#        self.elements.append(Image(filename, width=width, height=height))
#        self.temp_files.append(filename)
#
#        # Add flag plot
#        filename = get_temp_png_filename()
#        plot = PlotL1bdataWaveformFlags()
#        plot.create_plot(self.l1b)
#        plot.savefig(filename)
#
#        width = xsize*mm
#        height = plot.canvas_aspect*xsize*mm
#        self.elements.append(Image(filename, width=width, height=height))
#        self.temp_files.append(filename)

        # new page
        self.elements.append(PageBreak())

    def create_surface_type_page(self):

        # Settings
        xsize = self.figure_xsize

        # Add chapter title
        subtitle = Paragraph("Surface Type Classification",
                             self.style_subtitle)
        self.elements.append(subtitle)
        self.elements.append(Spacer(width=0, height=10*mm))

        # new page
        self.elements.append(PageBreak())

    def create_elevation_page(self):

        # Settings
        xsize = self.figure_xsize

        # Add chapter title
        subtitle = Paragraph("Surface Elevation",
                             self.style_subtitle)
        self.elements.append(subtitle)
        self.elements.append(Spacer(width=0, height=10*mm))

        # new page
        self.elements.append(PageBreak())

    def create_freeboard_thickness_page(self):

        # Settings
        xsize = self.figure_xsize

        # Add chapter title
        subtitle = Paragraph("Freeboard & Thickness",
                             self.style_subtitle)
        self.elements.append(subtitle)
        self.elements.append(Spacer(width=0, height=10*mm))

        # new page
        self.elements.append(PageBreak())

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

    page_width = A4[0]
    page_height = A4[1]
    canvas.saveState()
    canvas.setFont('Helvetica', 9)
    canvas.drawCentredString(page_width/2., 25, "Page %d" % (doc.page))
    canvas.drawCentredString(page_width/2.0, page_height-25, doc.l2i_filename)
    canvas.restoreState()


if __name__ == "__main__":
    l2idata_report()
