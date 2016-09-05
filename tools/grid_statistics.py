# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 16:40:54 2016

@author: shendric
"""

from pysiral.iotools import NCMaskedGridData
from pysiral.visualization.parameter import GridMapParameter


cs2awi_naming = {"sea_ice_freeboard": "freeboard",
                 "sea_ice_thickness": "sea_ice_thickness",
                 "sea_ice_concentration": "sea_ice_concentration",
                 "lead_fraction": "lead_fraction",
                 "sea_surface_height_anomaly": "sea_surface_anomaly"}

def grid_statistics():

    # get command line arguments
    parser = get_argparser()
    args = parser.parse_args()

    # Read the grod file
    ncdata = NCMaskedGridData(args.grid_file)

    # Get the parameter
    parameter_name = args.parameter_name
    data = GridMapParameter()
    data.set_grid(ncdata.longitude, ncdata.latitude)
    parameter = ncdata.get_by_name(parameter_name)
    if args.cs2awi:
        parameter_name = cs2awi_naming[parameter_name]
    data.set_parameter(parameter, parameter_name)

    import matplotlib.pyplot as plt

    cmap = data.get_cmap()

    plt.figure()
    plt.imshow(data.grid, interpolation="none", cmap=plt.get_cmap(cmap.name),
               vmin=cmap.vmin, vmax=cmap.vmax)
    plt.show()


def get_argparser():
    """ Handle command line arguments """

    import argparse

    parser = argparse.ArgumentParser()

    # Positional arguments
    parser.add_argument('-diff', action='store', default=None, dest='diff')
    parser.add_argument('-index-box', action='store', default=None,
                        dest='index_box', type=int, nargs=4)
    parser.add_argument('-lonlat-box', action='store', default=None,
                        dest='lonlat_box', type=int, nargs=4)
    parser.add_argument('--cs2awi', action='store_true', dest='cs2awi')
    parser.add_argument(action='store', dest='grid_file')
    parser.add_argument(action='store', dest='parameter_name')
    parser.add_argument(action='store', dest='output')

    return parser


if __name__ == "__main__":
    grid_statistics()
