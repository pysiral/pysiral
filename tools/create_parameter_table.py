# -*- coding: utf-8 -*-
"""
Created on Mon Jul 06 15:21:48 2015

@author: Stefan
"""

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4, landscape
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import mm

from pysiral.config import ConfigInfo

import os


def create_parameter_table():
    """ Reads the parameter definition file and creates overview table """

    info = ConfigInfo()

    filename = os.path.join(
        info.doc_path, "parameter", "pysiral_parameter_table.pdf")
    doc = SimpleDocTemplate(filename, pagesize=A4)
    doc.pagesize = landscape(A4)
    elements = []

    # Styles
    styles = getSampleStyleSheet()
    styleH = styles['Heading1']

    # Table Title
    elements.append(Paragraph("Parameter Definition in pysiral", styleH))

    data = []
    header = ["short_name",
              "long_name",
              "level",
              "dtype",
              "format",
              "unit",
              "docstr"]
    data.append(header)

    # Add parameters
    for key in sorted(info.parameter.definition.keys(branch_mode='only')):
        entry = [
            info.parameter.definition[key].short_name,
            key,
            info.parameter.definition[key].level,
            info.parameter.definition[key].ascii_format,
            info.parameter.definition[key].dtype,
            info.parameter.definition[key].unit,
            Paragraph(
                info.parameter.definition[key].docstr, styles['Normal'])
        ]
        data.append(entry)

    # Table Style
    style = TableStyle([
        ('FONT', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONT', (0, 1), (-1, -1), 'Helvetica'),
        ('FONTSIZE', (0, 0), (-1, -1), 10),
        ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
        ('VALIGN', (0, 0), (-1, -1), 'TOP'),
        ('INNERGRID', (0, 0), (-1, -1), 0.50, colors.black),
        ('BOX', (0, 0), (-1, -1), 0.25, colors.black)])

    col_widths = ('*', '*', '*', '*', '*', '*', 100*mm)
    t = Table(data, hAlign='CENTER', colWidths=col_widths, repeatCols=1)
    t.setStyle(style)

    # Send the data and build the file
    elements.append(t)
    doc.build(elements)

if __name__ == "__main__":
    create_parameter_table()
