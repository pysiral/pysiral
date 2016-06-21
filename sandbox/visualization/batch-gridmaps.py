# -*- coding: utf-8 -*-
"""
Created on Thu Apr 07 19:59:06 2016

@author: Stefan
"""

from pysiral.proj import EASE2North
from pysiral.iotools import NCMaskedGridData
from pysiral.logging import stdout_logger
from pysiral.path import filename_from_path, file_basename
from pysiral.visualization.gridmap import ArcticGridPresentationMap
from pysiral.visualization.mapstyle import GridMapAWILightStyle
from pysiral.visualization.parameter import (GridMapParameter,
                                             GridMapDiffParameter)

import os
import glob
import numpy as np
from datetime import datetime
from dateutil.relativedelta import relativedelta


MONTH = ["January", "February", "March", "April", "May", "June", "July",
         "August", "September", "October", "November", "December"]


def batch_gridmaps_cs2_monthly():
    """ Plot a series of parameters from CS-2 monthly grid files """

    # Parameter list dict (key: name in netCDF file, value: standardized name)
    parameter_list = {
        "sea_ice_thickness": "sea_ice_thickness"}

    # Projection info of the CryoSat-2 data
    projection_keyw = EASE2North().projection_keyw

    # Create the logging instance
    log = stdout_logger("batch-gridmaps-cs2-monthly")

    # Get list of input files
    folder = r"E:\altim\product\altimetry\cs2awi\cs2awi_north_version_1_2\grid"
    pattern = "*.nc"
    log.info("Searching folder: %s" % folder)
    ncfiles = sorted(glob.glob(os.path.join(folder, pattern)))
    log.info("Found %g netCDF files" % len(ncfiles))

    for i, ncfile in enumerate(ncfiles):

        # Read the file
        log.info("Read file (%g of %g) : %s" % (
            i+1, len(ncfiles), filename_from_path(ncfile)))
        ncdata = NCMaskedGridData(ncfile)

        # Loop over parameters
        for parameter_name in parameter_list.keys():

            standard_name = parameter_list[parameter_name]
            log.info("  Create map: %s" % standard_name)

            # Extract parameter and merge with metadata
            data = GridMapParameter()
            data.set_grid(ncdata.longitude, ncdata.latitude)
            parameter = ncdata.get_by_name(parameter_name)
            data.set_parameter(parameter, parameter_name)
            data.set_projection(**projection_keyw)
            # Light style map
            map_filename = file_basename(ncfile)+"_"+standard_name+".png"
            output = os.path.join(folder, "visuals", map_filename)

            gridmap = ArcticGridPresentationMap()
            gridmap.style = GridMapAWILightStyle()
            gridmap.data = data
            gridmap.label.title = "CryoSat-2"
            gridmap.label.period = get_cs2_period_from_filename(ncfile)
            gridmap.save2png(output)
            log.info("  Map figure written: %s" % output)


def batch_gridmaps_cs2_diff_previous_month():
    """ Plot a series of differences from CS-2 monthly grid files """

    # Parameter list dict (key: name in netCDF file, value: standardized name)
    parameter_list = {
        "sea_ice_thickness": "sea_ice_thickness"}

    # Projection info of the CryoSat-2 data
    projection_keyw = EASE2North().projection_keyw

    # Create the logging instance
    log = stdout_logger("batch-gridmaps-cs2-previous-month")

    # Get list of input files
    folder = r"E:\altim\product\altimetry\cs2awi\cs2awi_north_version_1_2\grid"

    for yyyy in np.arange(2015, 2016):
        for mm_inc in np.arange(11, 18):
            mb = datetime(yyyy, 1, 1) + relativedelta(months=mm_inc)
            ma = mb - relativedelta(months=1)
            # Read the file
            log.info("Read files for %4g%02g minus %4g%02g" % (
                mb.year, mb.month, ma.year, ma.month))
            # Read the later file
            ncfile_a = "cs2awi_nh_%4g%02g.nc" % (ma.year, ma.month)
            path_a = os.path.join(folder, ncfile_a)
            ncfile_b = "cs2awi_nh_%4g%02g.nc" % (mb.year, mb.month)
            path_b = os.path.join(folder, ncfile_b)

            # check if source files exist
            if not os.path.isfile(path_a) or not os.path.isfile(path_b):
                log.info("Source files do not exist, skipping")
                continue

            log.info("Read grid_a: %s" % ncfile_a)
            ncdata_a = NCMaskedGridData(path_a)
            # Read the earlier file
            log.info("Read grid_b: %s" % ncfile_b)
            ncdata_b = NCMaskedGridData(path_b)

            # Loop over parameters
            for parameter_name in parameter_list.keys():

                standard_name = parameter_list[parameter_name]
                log.info("  Create map: %s" % standard_name)

                # Extract parameter and merge with metadata
                data = GridMapDiffParameter()
                data.set_grid(ncdata_a.longitude, ncdata_a.latitude)
                parameter_a = ncdata_a.get_by_name(parameter_name)
                parameter_b = ncdata_b.get_by_name(parameter_name)
                data.set_parameter(parameter_a, parameter_b, parameter_name)
                data.set_projection(dx=25000, dy=25000, **projection_keyw)
                # Light style map
                diff_str = "cs2awi_nh_%4g%02g_diff_%4g%02g_%s" % (
                    mb.year, mb.month, ma.year, ma.month, standard_name)
                map_filename = diff_str+"_"+standard_name+".png"
                output = os.path.join(
                    folder, "visuals", "diff_previous_month", map_filename)

                gridmap = ArcticGridPresentationMap()
                gridmap.style = GridMapAWILightStyle()
                gridmap.data = data
                gridmap.label.title = "CryoSat-2"
                gridmap.label.period = get_cs2_diff_period(ma, mb)
                gridmap.save2png(output)
                log.info("  Map figure written: %s" % output)


def batch_gridmaps_cs2_diff_previous_year():
    """ Plot a series of differences from CS-2 monthly grid files """

    # Parameter list dict (key: name in netCDF file, value: standardized name)
    parameter_list = {
        "sea_ice_thickness": "sea_ice_thickness"}

    # Projection info of the CryoSat-2 data
    projection_keyw = EASE2North().projection_keyw

    # Create the logging instance
    log = stdout_logger("batch-gridmaps-cs2-previous-year")

    # Get list of input files
    folder = r"E:\altim\product\altimetry\cs2awi\cs2awi_north_version_1_2\grid"

    for yyyy in np.arange(2015, 2016):
        for mm_inc in np.arange(10, 18):
            mb = datetime(yyyy, 1, 1) + relativedelta(months=mm_inc)
            ma = mb - relativedelta(years=1)
            # Read the file
            log.info("Read files for %4g%02g minus %4g%02g" % (
                mb.year, mb.month, ma.year, ma.month))
            # Read the later file
            ncfile_a = "cs2awi_nh_%4g%02g.nc" % (ma.year, ma.month)
            path_a = os.path.join(folder, ncfile_a)
            ncfile_b = "cs2awi_nh_%4g%02g.nc" % (mb.year, mb.month)
            path_b = os.path.join(folder, ncfile_b)

            # check if source files exist
            if not os.path.isfile(path_a) or not os.path.isfile(path_b):
                log.info("Source files do not exist, skipping")
                continue

            log.info("Read grid_a: %s" % ncfile_a)
            ncdata_a = NCMaskedGridData(path_a)
            # Read the earlier file
            log.info("Read grid_b: %s" % ncfile_b)
            ncdata_b = NCMaskedGridData(path_b)

            # Loop over parameters
            for parameter_name in parameter_list.keys():

                standard_name = parameter_list[parameter_name]
                log.info("  Create map: %s" % standard_name)

                # Extract parameter and merge with metadata
                data = GridMapDiffParameter()
                data.set_grid(ncdata_a.longitude, ncdata_a.latitude)
                parameter_a = ncdata_a.get_by_name(parameter_name)
                parameter_b = ncdata_b.get_by_name(parameter_name)
                data.set_parameter(parameter_a, parameter_b, parameter_name)
                data.set_projection(dx=25000, dy=25000, **projection_keyw)
                # Light style map
                diff_str = "cs2awi_nh_%4g%02g_diff_%4g%02g_%s" % (
                    mb.year, mb.month, ma.year, ma.month,standard_name)
                map_filename = diff_str+"_"+standard_name+".png"
                output = os.path.join(
                    folder, "visuals", "diff_previous_year", map_filename)

                gridmap = ArcticGridPresentationMap()
                gridmap.style = GridMapAWILightStyle()
                gridmap.data = data
                gridmap.label.title = "CryoSat-2"
                gridmap.label.period = get_cs2_diff_period(ma, mb)
                gridmap.save2png(output)
                log.info("  Map figure written: %s" % output)


def batch_gridmaps_cs2_diff_winter_2011_2012():
    """ Plot a series of differences from CS-2 monthly grid files """

    # Parameter list dict (key: name in netCDF file, value: standardized name)
    parameter_list = {
        "sea_ice_thickness": "sea_ice_thickness"}

    # Projection info of the CryoSat-2 data
    projection_keyw = EASE2North().projection_keyw

    # Create the logging instance
    log = stdout_logger("batch-gridmaps-cs2-diff-winter-2011-2012")

    # Get list of input files
    folder = r"E:\altim\product\altimetry\cs2awi\cs2awi_north_version_1_2\grid"

    for yyyy in [2015]:
        for mm_inc in np.arange(10, 18):
            mb = datetime(yyyy, 1, 1) + relativedelta(months=mm_inc)
            ma = datetime(2011, 1, 1) + relativedelta(months=mm_inc)
            # Read the file
            log.info("Read files for %4g%02g minus %4g%02g" % (
                mb.year, mb.month, ma.year, ma.month))
            # Read the later file
            ncfile_a = "cs2awi_nh_%4g%02g.nc" % (ma.year, ma.month)
            path_a = os.path.join(folder, ncfile_a)
            ncfile_b = "cs2awi_nh_%4g%02g.nc" % (mb.year, mb.month)
            path_b = os.path.join(folder, ncfile_b)

            # check if source files exist
            if not os.path.isfile(path_a) or not os.path.isfile(path_b):
                log.info("Source files do not exist, skipping")
                continue

            log.info("Read grid_a: %s" % ncfile_a)
            ncdata_a = NCMaskedGridData(path_a)
            # Read the earlier file
            log.info("Read grid_b: %s" % ncfile_b)
            ncdata_b = NCMaskedGridData(path_b)

            # Loop over parameters
            for parameter_name in parameter_list.keys():

                standard_name = parameter_list[parameter_name]
                log.info("  Create map: %s" % standard_name)

                # Extract parameter and merge with metadata
                data = GridMapDiffParameter()
                data.set_grid(ncdata_a.longitude, ncdata_a.latitude)
                parameter_a = ncdata_a.get_by_name(parameter_name)
                parameter_b = ncdata_b.get_by_name(parameter_name)
                data.set_parameter(parameter_a, parameter_b, parameter_name)
                data.set_projection(dx=25000, dy=25000, **projection_keyw)
                # Light style map
                diff_str = "cs2awi_nh_%4g%02g_diff_%4g%02g_%s" % (
                    mb.year, mb.month, ma.year, ma.month,standard_name)
                map_filename = diff_str+"_"+standard_name+".png"
                output = os.path.join(
                    folder, "visuals", "diff_winter_2011_2012", map_filename)

                gridmap = ArcticGridPresentationMap()
                gridmap.style = GridMapAWILightStyle()
                gridmap.data = data
                gridmap.label.title = "CryoSat-2"
                gridmap.label.period = get_cs2_diff_period(ma, mb)
                gridmap.save2png(output)
                log.info("  Map figure written: %s" % output)


def batch_gridmaps_cs2_diff_all_winters():
    """ Plot a series of differences from CS-2 monthly grid files """

    # Parameter list dict (key: name in netCDF file, value: standardized name)
    parameter_list = {
        "sea_ice_thickness": "sea_ice_thickness"}

    # Projection info of the CryoSat-2 data
    projection_keyw = EASE2North().projection_keyw

    # Create the logging instance
    log = stdout_logger("batch-gridmaps-cs2-diff-winter-2011-2012")

    # Get list of input files
    folder = r"E:\altim\product\altimetry\cs2awi\cs2awi_north_version_1_2\grid"

    for yyyy in [2015]:
        for mm_inc in np.arange(10, 18):

            ma = datetime(2015, 1, 1) + relativedelta(months=mm_inc)
            # Read the file
            log.info("Read files for %4g%02g minus Winter Mean" % (
                ma.year, ma.month))
            # Read the later file
            ncfile_a = "cs2awi_nh_%4g%02g.nc" % (ma.year, ma.month)
            path_a = os.path.join(folder, ncfile_a)

            log.info("Read grid_a: %s" % ncfile_a)
            ncdata_a = NCMaskedGridData(path_a)

            # Loop over parameters
            for parameter_name in parameter_list.keys():

                standard_name = parameter_list[parameter_name]
                log.info("  Create map: %s" % standard_name)

                # Extract parameter and merge with metadata
                data = GridMapDiffParameter()
                data.set_grid(ncdata_a.longitude, ncdata_a.latitude)
                parameter_a = ncdata_a.get_by_name(parameter_name)

                # Read the earlier file
                parameter_b = None
                count = 0
                for yyyy in [2010, 2011, 2012, 2013, 2014]:
                    mb = datetime(yyyy, 1, 1) + relativedelta(months=mm_inc)
                    ncfile_b = "cs2awi_nh_%4g%02g.nc" % (mb.year, mb.month)
                    path_b = os.path.join(folder, ncfile_b)
                    log.info("Read grid_b: %s" % ncfile_b)
                    ncdata_b = NCMaskedGridData(path_b)
                    parameter_yyyy = ncdata_b.get_by_name(parameter_name)
                    if parameter_b is None:
                        parameter_b = parameter_yyyy
                    else:
                        parameter_b = parameter_b + parameter_yyyy
                    count = count + 1
                parameter_b = parameter_b/float(count)

                data.set_parameter(parameter_b, parameter_a, parameter_name)
                data.set_projection(dx=25000, dy=25000, **projection_keyw)
                # Light style map
                diff_str = "cs2awi_nh_%4g%02g_diff_2010_2015_mean_%s" % (
                    mb.year, mb.month, standard_name)
                map_filename = diff_str+"_"+standard_name+".png"
                output = os.path.join(
                    folder, "visuals", "diff_all_winters", map_filename)

                gridmap = ArcticGridPresentationMap()
                gridmap.style = GridMapAWILightStyle()
                gridmap.data = data
                gridmap.label.title = "CryoSat-2"
                gridmap.label.period = get_cs2_diff_mean(ma)
                gridmap.save2png(output)
                log.info("  Map figure written: %s" % output)


def get_cs2_period_from_filename(ncfile):
    basename = file_basename(ncfile)
    yyyymm = basename.split("_")[2]
    return MONTH[int(yyyymm[4:6])-1]+" "+yyyymm[0:4]


def get_cs2_diff_period(ma, mb):
    return "%s %4g - %s %4g" % (MONTH[mb.month-1], mb.year,
                                MONTH[ma.month-1], ma.year)


def get_cs2_diff_mean(ma):
    yyyy0 = 2010 + ma.year-2015
    yyyy1 = yyyy0 + 4
    return "%s %4g - Mean (%4g to %4g)" % (
        MONTH[ma.month-1], ma.year, yyyy0, yyyy1)

if __name__ == "__main__":
    batch_gridmaps_cs2_monthly()
    # batch_gridmaps_cs2_diff_previous_month()
    # batch_gridmaps_cs2_diff_previous_year()
    # batch_gridmaps_cs2_diff_winter_2011_2012()
    # batch_gridmaps_cs2_diff_all_winters()
