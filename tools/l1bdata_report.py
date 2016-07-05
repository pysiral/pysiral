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
from pysiral.visualization.l1bvis import (
    PlotL1bdataOrbit, PlotL1bdataWaveform, PlotL1bdataWaveformFlags,
    PlotL1bdataRangeCorrections, PlotL1bdataWaveformClassifier)

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4, portrait
from reportlab.platypus import (SimpleDocTemplate, Table, TableStyle,
                                Paragraph, PageBreak, Image, Spacer)
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import mm


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
            pdf.create_classifier_pages()
            pdf.create_range_correction_pages()
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
        self.doc.l1b_filename = self.l1b_filename
        self.elements.append(Paragraph("l1bdata Report", self.style_heading))
        self.elements.append(Paragraph(self.l1b_filename, self.style_subtitle))
        self.elements.append(Spacer(width=0, height=5*mm))

        filename = get_temp_png_filename()
        plot = PlotL1bdataOrbit()
        plot.filename = filename
        plot.create_plot(self.l1b)
        plot.savefig()
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

        # Settings
        xsize = self.figure_xsize

        # Add chapter title
        subtitle = Paragraph("Waveform Parameter", self.style_subtitle)
        self.elements.append(subtitle)
        self.elements.append(Spacer(width=0, height=10*mm))

        # Add waveform plot
        filename = get_temp_png_filename()
        plot = PlotL1bdataWaveform()
        plot.filename = filename
        plot.create_plot(self.l1b)
        plot.savefig()
        width = xsize*mm
        height = plot.canvas_aspect*xsize*mm
        self.elements.append(Image(filename, width=width, height=height))
        self.temp_files.append(filename)

        # Add flag plot
        filename = get_temp_png_filename()
        plot = PlotL1bdataWaveformFlags()
        plot.filename = filename
        plot.create_plot(self.l1b)
        plot.savefig()
        width = xsize*mm
        height = plot.canvas_aspect*xsize*mm
        self.elements.append(Image(filename, width=width, height=height))
        self.temp_files.append(filename)

        # new page
        self.elements.append(PageBreak())

    def create_range_correction_pages(self):

        # Settings
        xsize = self.figure_xsize

        # Add range corrections plot
        max_grc_plots_per_page = 5
        n_grc_parameter = len(self.l1b.correction.parameter_list)
        n_grc_plots = np.ceil(n_grc_parameter/max_grc_plots_per_page)

        for plotn in np.arange(n_grc_plots):

            # Add chapter title
            label = "Range Corrections"
            if n_grc_plots > 1:
                label = label + " (%g of %g)" % (plotn+1, n_grc_plots)
            subtitle = Paragraph(label, self.style_subtitle)
            self.elements.append(subtitle)

            # Decide which range corrections to plot
            if plotn < n_grc_plots-1:
                grc_indices = np.arange(
                    plotn*max_grc_plots_per_page,
                    (plotn+1)*max_grc_plots_per_page).astype(int)
            else:
                grc_indices = np.arange(
                    plotn*max_grc_plots_per_page,
                    n_grc_parameter).astype(int)

            # Create the plot and write to temp file
            filename = get_temp_png_filename()
            plot = PlotL1bdataRangeCorrections()
            plot.filename = filename
            plot.create_plot(self.l1b, grc_indices=grc_indices,
                             max_number_plots=max_grc_plots_per_page)
            plot.savefig()

            # Add the plot to the pdf
            width = xsize*mm
            height = plot.canvas_aspect*xsize*mm
            self.elements.append(Image(filename, width=width, height=height))

            # new page
            self.elements.append(PageBreak())

            # Add temp files (figures) to clean up list
            self.temp_files.append(filename)

    def create_classifier_pages(self):

        # Settings
        xsize = self.figure_xsize

        # Add range corrections plot
        max_clf_plots_per_page = 5
        n_clf_parameter = len(self.l1b.classifier.parameter_list)
        n_clf_plots = np.ceil(n_clf_parameter/max_clf_plots_per_page)

        for plotn in np.arange(n_clf_plots):

            # Add chapter title
            label = "Waveform Classifier"
            if n_clf_plots > 1:
                label = label + " (%g of %g)" % (plotn+1, n_clf_plots)
            subtitle = Paragraph(label, self.style_subtitle)
            self.elements.append(subtitle)

            # Decide which range corrections to plot
            if plotn < n_clf_plots-1:
                grc_indices = np.arange(
                    plotn*max_clf_plots_per_page,
                    (plotn+1)*max_clf_plots_per_page).astype(int)
            else:
                grc_indices = np.arange(
                    plotn*max_clf_plots_per_page,
                    n_clf_parameter).astype(int)

            # Create the plot and write to temp file
            filename = get_temp_png_filename()
            plot = PlotL1bdataWaveformClassifier()
            plot.filename = filename
            plot.create_plot(self.l1b, clf_indices=grc_indices,
                             max_number_plots=max_clf_plots_per_page)
            plot.savefig()

            # Add the plot to the pdf
            width = xsize*mm
            height = plot.canvas_aspect*xsize*mm
            self.elements.append(Image(filename, width=width, height=height))

            # new page
            self.elements.append(PageBreak())

            # Add temp files (figures) to clean up list
            self.temp_files.append(filename)

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
    canvas.drawCentredString(page_width/2.0, page_height-25, doc.l1b_filename)
    canvas.restoreState()


if __name__ == "__main__":
    l1bdata_report()
