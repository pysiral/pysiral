# -*- coding: utf-8 -*-
"""
Created on Sun Apr 10 15:52:01 2016

@author: Stefan
"""

from pysiral.proj import EASE2North
from pysiral.iotools import NCMaskedGridData
from pysiral.visualization.parameter import GridMapParameter

from pysiral.path import (file_basename, folder_from_filename,
                          validate_directory)
import os


mission_name_dict = {"cryosat2": "CryoSat-2", "envisat": "Envisat",
                     "ers2": "ERS-2"}

cs2awi_naming = {"sea_ice_freeboard": "freeboard",
                 "sea_ice_thickness": "sea_ice_thickness",
                 "sea_ice_concentration": "sea_ice_concentration",
                 "lead_fraction": "lead_fraction",
                 "sea_surface_height_anomaly": "sea_surface_anomaly"}

month_names = {
    "1": "January", "2": "February", "3": "March", "4": "April", "5": "May",
    "6": "June", "7": "July", "8": "August", "9": "September", "10": "October",
    "11": "November", "12": "December"}


def l3s_map():

    parser = get_l3s_map_argparser()
    args = parser.parse_args()

    # TODO: batch processing
    ncfile = args.l3s_filename
    ncdata = NCMaskedGridData(ncfile)

    # Loop over parameters
    for parameter_name in args.parameter_list:

        # Extract parameter and merge with metadata
        data = GridMapParameter()
        data.set_grid(ncdata.longitude, ncdata.latitude)
        parameter = ncdata.get_by_name(parameter_name)
        if args.cs2awi:
            parameter_name = cs2awi_naming[parameter_name]
        data.set_parameter(parameter, parameter_name)
        data.set_projection(**EASE2North().projection_keyw)

        if args.cs2awi:
            data.set_nan_mask(ncdata.sea_ice_freeboard)

        # Set output folder
        if args.destination is None:
            destination = folder_from_filename(args.l3s_filename)
            destination = os.path.join(destination, "visuals", args.maptype)
        else:
            destination = args.destination
        validate_directory(destination)

        # Output filename
        map_filename = file_basename(ncfile)+"_"+data.short_name+".png"
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

        # Light style map
        MapClass, StyleClass = get_map_classes(args)
        gridmap = MapClass()
        gridmap.style = StyleClass()
        gridmap.data = data
        gridmap.label.title = mission_name_dict[mission_id]
        gridmap.label.period = period_label
        gridmap.label.annotation = args.annotation
        gridmap.save2png(output)


def get_map_classes(args):

    from pysiral.visualization.gridmap import (
        ArcticGridPresentationMap, ArcticGridPaperMap)
    from pysiral.visualization.mapstyle import (
        GridMapAWILightStyle, GridMapPaperStyle)
    import sys

    if args.maptype == "presentation":
        return ArcticGridPresentationMap, GridMapAWILightStyle
    elif args.maptype == "paper":
        return ArcticGridPaperMap, GridMapPaperStyle
    else:
        sys.exit("Invalid map type (%s), aborting ..." % args.maptype)


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

    parser.add_argument(
        '-annotation',
        action='store', dest='annotation', default='',
        help='additional label')

    # Batch Processing
    parser.add_argument(
        '--batch',
        action='store_true', dest='batch',
        help='positional argument is a search pattern')

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
