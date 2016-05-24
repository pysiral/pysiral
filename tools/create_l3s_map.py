# -*- coding: utf-8 -*-
"""
Created on Sun Apr 10 15:52:01 2016

@author: Stefan
"""

from pysiral.proj import EASE2North
from pysiral.iotools import NCMaskedGridData
from pysiral.visualization.gridmap import ArcticGridPresentationMap
from pysiral.visualization.mapstyle import (GridMapSICCILightStyle)
from pysiral.visualization.parameter import GridMapParameter

from pysiral.path import (file_basename, folder_from_filename,
                          validate_directory)
import os


mission_name_dict = {"cryosat2": "CryoSat-2", "envisat": "Envisat",
                     "ers2": "ERS-2"}

cs2awi_naming = {"sea_ice_freeboard": "freeboard",
                 "sea_ice_thickness": "sea_ice_thickness",
                 "sea_ice_concentration": "sea_ice_concentration"}

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
            period_label = ""
        else:
            mission_id = ncdata.mission_ids
            period_label = ncdata.period_label

        # Light style map
        gridmap = ArcticGridPresentationMap()
        gridmap.style = GridMapSICCILightStyle()
        gridmap.data = data
        gridmap.label.title = mission_name_dict[mission_id]
        gridmap.label.period = period_label
        gridmap.save2png(output)


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
