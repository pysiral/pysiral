# -*- coding: utf-8 -*-
"""
Created on Sun Apr 10 15:52:01 2016

@author: Stefan
"""

from pysiral.proj import EASE2North, EASE2South
from pysiral.iotools import NCMaskedGridData
from pysiral.visualization.parameter import (
    GridMapParameter, GridMapDiffParameter)

from pysiral.path import (file_basename, folder_from_filename,
                          validate_directory)
import glob
import os


mission_name_dict = {"cryosat2": "CryoSat-2", "envisat": "Envisat",
                     "ers2": "ERS-2", "sentinel3a": "Sentinel-3A"}

cs2awi_naming = {"sea_ice_freeboard": "freeboard",
                 "sea_ice_thickness": "sea_ice_thickness",
                 "sea_ice_concentration": "sea_ice_concentration",
                 "lead_fraction": "lead_fraction",
                 "sea_surface_height_anomaly": "sea_surface_anomaly"}

month_names = {
    "1": "January", "2": "February", "3": "March", "4": "April", "5": "May",
    "6": "June", "7": "July", "8": "August", "9": "September", "10": "October",
    "11": "November", "12": "December"}


PROJDICT = {"north": EASE2North(), "south": EASE2South()}


def l3s_map():

    parser = get_l3s_map_argparser()
    args = parser.parse_args()

    # Get files (can be file link or search pattern)
    ncfiles = sorted(glob.glob(args.l3s_filename))

    print "%g l3s files found" % len(ncfiles)

    # XXX: Temporary cs2awi fix
    xres, yres = None, None
    if args.cs2awi:
        xres, yres = 25000, 25000

    for ncfile in ncfiles:

        # TODO: batch processing
        ncdata = NCMaskedGridData(ncfile)

        try:
            hemisphere = ncdata.hemisphere
        except:
            hemisphere = "north"

        # Get sea ice concentration as background
        try:
            sic = GridMapParameter()
            sic.set_grid(ncdata.longitude, ncdata.latitude)
            sic.set_parameter(ncdata.get_by_name("sea_ice_concentration"),
                              "sea_ice_concentration")
        except:
            sic = None

        # Loop over parameters
        for parameter_name in args.parameter_list:

            # Extract parameter and merge with metadata
            if args.diff == "none":
                data = GridMapParameter()
                data.set_grid(ncdata.longitude, ncdata.latitude)
                parameter = ncdata.get_by_name(parameter_name)
                if args.cs2awi:
                    parameter_name = cs2awi_naming[parameter_name]
                data.set_parameter(parameter, parameter_name)

            else:
                ncdata_diff = NCMaskedGridData(args.diff)
                data = GridMapDiffParameter()
                data.set_grid(ncdata.longitude, ncdata.latitude)
                parameter_a = ncdata.get_by_name(parameter_name)
                parameter_b = ncdata_diff.get_by_name(parameter_name)
                if args.cs2awi:
                    parameter_name = cs2awi_naming[parameter_name]
                data.set_parameter(parameter_a, parameter_b, parameter_name)

            projection = PROJDICT[hemisphere]
            data.set_projection(xres=xres, yres=yres,
                                **projection.projection_keyw)

            if args.cs2awi:
                data.set_nan_mask(ncdata.sea_ice_freeboard)

            # Set output folder
            if args.destination is None:
                destination = folder_from_filename(args.l3s_filename)
                destination = os.path.join(
                    destination, "visuals", args.maptype)
            else:
                destination = args.destination
            validate_directory(destination)

            # Output filename
            annotation_str = get_annotation_str(args.annotation)
            map_filename = file_basename(ncfile) + "_" + data.short_name + \
                annotation_str + ".png"
            output = os.path.join(destination, map_filename)

            # Map Labels
            if args.cs2awi:
                mission_id = "cryosat2"
                # TODO: This is preliminary
                period = ncdata.description.split()[-1][1:-1]
                month, year = period.split("/")
                period_label = "%s %s" % (month_names[month], year)
            else:
                mission_id = ncdata.mission_ids
                period_label = ncdata.period_label

            # Map Title
            title = mission_name_dict[mission_id]
            try:  # Will fails for cs2awi files
                if (args.diff != "none") and (
                        ncdata_diff.mission_ids != ncdata_diff.mission_ids):
                    mission_id_b = ncdata_diff.mission_ids
                    title = mission_name_dict[mission_id_b] + " - " + \
                        mission_name_dict[mission_id]
            except:
                title = "CryoSat-2"

            # Light style map
            MapClass, StyleClass = get_map_classes(args, hemisphere)
            gridmap = MapClass()
            gridmap.style = StyleClass()
            gridmap.data = data
            gridmap.sic = sic
            gridmap.label.title = title
            gridmap.label.period = period_label
            gridmap.label.annotation = args.annotation
            gridmap.save2png(output)


def get_map_classes(args, hemisphere):

    from pysiral.visualization.gridmap import (
        ArcticGridPresentationMap, ArcticGridPaperMap,
        AntarcticGridPresentationMap, AntarcticGridPaperMap)
    from pysiral.visualization.mapstyle import (
        GridMapAWILightStyle, GridMapPaperStyle)
    import sys

    if args.maptype not in ["presentation", "paper"]:
        sys.exit("Invalid map type (%s), aborting ..." %
                 args.maptype)

    if hemisphere == "north":
        if args.maptype == "presentation":
            return ArcticGridPresentationMap, GridMapAWILightStyle
        elif args.maptype == "paper":
            return ArcticGridPaperMap, GridMapPaperStyle
    elif hemisphere == "south":
        if args.maptype == "presentation":
            return AntarcticGridPresentationMap, GridMapAWILightStyle
        elif args.maptype == "paper":
            return AntarcticGridPaperMap, GridMapPaperStyle
    else:
        sys.exit("Invalid hemisphere definition (%s), aborting ..." %
                 hemisphere)


def get_annotation_str(annotation):
    annotation_str = annotation.lower()
    annotation_str = annotation_str.replace(" ", "_")
    if annotation_str != "":
        annotation_str = "_" + annotation_str
    return annotation_str


def get_l3s_map_argparser():
    """ Handle command line arguments """
    import argparse

    parser = argparse.ArgumentParser()

    # Positional arguments
    parser.add_argument(
        action='store', dest='l3s_filename',
        help='l3s netcdf or search pattern if --batch is set')

    parser.add_argument(
        '-parameter', action='store', dest='parameter_list',
        nargs='+', required=True,
        help='list of parameter')

    parser.add_argument(
        '-o', '-output',
        action='store', dest='destination',
        help=r'output folder for maps (default: input\visuals)')

    # options
    parser.add_argument(
        '-type',
        choices=['presentation', 'paper'],
        action='store', dest='maptype', default='presentation',
        help='map type')

    # options
    parser.add_argument(
        '-diff',
        action='store', dest='diff', default='none',
        help='reference grid')

    parser.add_argument(
        '-annotation',
        action='store', dest='annotation', default='',
        help='additional label')

    parser.add_argument(
        '--cs2awi',
        action='store_true', dest='cs2awi',
        help='source file is from IDL cs2awi processor')

    # Add colobar
    cb_parser = parser.add_mutually_exclusive_group(required=False)
    cb_parser.add_argument('--cb', dest='colorbar', action='store_true')
    cb_parser.add_argument('--no-cb', dest='colorbar', action='store_false')
    parser.set_defaults(colorbar=True)

    # show preprocessor version
    parser.add_argument(
        '--version', action='version', version='%(prog)s 0.1a')

    return parser

if __name__ == "__main__":
    l3s_map()
