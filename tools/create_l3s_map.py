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


def l3s_map():

    parser = get_l3s_map_argparser()
    args = parser.parse_args()

    parameter_list = ["sea_ice_thickness", "freeboard"]

    mission_name_dict = {"cryosat2": "CryoSat-2", "envisat": "Envisat",
                         "ers2": "ERS-2"}

    ncfile = args.l3s_filename
    ncdata = NCMaskedGridData(ncfile)

    # Loop over parameters
    for parameter_name in parameter_list:

        # Extract parameter and merge with metadata
        data = GridMapParameter()
        data.set_grid(ncdata.longitude, ncdata.latitude)
        parameter = ncdata.get_by_name(parameter_name)
        data.set_parameter(parameter, parameter_name)
        data.set_projection(**EASE2North().projection_keyw)

        # Set output folder
        if args.destination is None:
            destination = folder_from_filename(args.l3s_filename)
            destination = os.path.join(destination, "visuals")
        else:
            destination = args.destination
        validate_directory(destination)

        # Light style map
        map_filename = file_basename(ncfile)+"_"+data.short_name+".png"
        output = os.path.join(destination, map_filename)
        gridmap = ArcticGridPresentationMap()
        gridmap.style = GridMapSICCILightStyle()
        gridmap.data = data
        gridmap.label.title = mission_name_dict[ncdata.mission_ids]
        gridmap.label.period = ncdata.period_label
        gridmap.save2png(output)


def get_l3s_map_argparser():
    """ Handle command line arguments """
    import argparse

    parser = argparse.ArgumentParser()
    # Mission id string: cryosat2, envisat, ...
    parser.add_argument(
        '-i', '-l3sfile', action='store', dest='l3s_filename',
        help='l3s netcdf file path', required=True)
    parser.add_argument(
        '-o', '-destination', action='store', dest='destination',
        help='destination folder for maps')
    # show preprocessor version
    parser.add_argument(
        '--version', action='version', version='%(prog)s 0.1a')

    return parser

if __name__ == "__main__":
    l3s_map()
