# -*- coding: utf-8 -*-
"""
Created on Thu Feb 02 17:32:48 2017

@author: shendric
"""

from pysiral.l2data import L2iNCFileImport
from pysiral.path import (filename_from_path, file_basename,
                          validate_directory)
from pysiral.visualization.l2ivis import (PlotL2iDataOrbit,
                                          PlotL2iSeaSurfaceHeightFreeboard,
                                          PlotL2iFreeboardThickness)

import matplotlib.pyplot as plt

import numpy as np

import argparse
import glob
import os


def find_anomalous_frb_l2i():
    """ Look for anomalous l2i files that cause trackiness in l3s grids """

    args = get_command_line_arguments()

    # loop over files
    l2i_files = glob.glob(os.path.join(args.l2i_folder, "l2i*.nc"))

    # Create graph export folder
    export_directory = os.path.join(args.l2i_folder, "l2i_check")
    validate_directory(export_directory)

    ssa_var = np.full(len(l2i_files), np.nan)
    ssa_mean = np.full(len(l2i_files), np.nan)
    frb_mean = np.full(len(l2i_files), np.nan)
    n_anomalous_frb = np.full(len(l2i_files), np.nan)
    n_anomalous_frb_fract = np.full(len(l2i_files), np.nan)
    frb2sit = np.full(len(l2i_files), np.nan)

    for i, l2i_file in enumerate(l2i_files):
        l2i = L2iNCFileImport(l2i_file)
        ssa_var[i] = np.ptp(l2i.sea_surface_anomaly)
        ssa_mean[i] = np.nanmean(l2i.sea_surface_anomaly)
        frb_mean[i] = np.nanmean(l2i.freeboard)
        frb2sit[i] = np.nanmean(l2i.freeboard)/np.nanmean(l2i.sea_ice_thickness)

        n_anomalous_frb[i] = len(np.where(l2i.freeboard > 0.4)[0])
        n_anomalous_frb_fract[i] = n_anomalous_frb[i]/l2i.n_records

#        if np.abs(frb2sit[i]) > 0.5:
#            print filename_from_path(l2i_file)
#            plt.figure("freeboard")
#            plt.plot(l2i.freeboard)
#            plt.figure("sea_ice_thickness")
#            plt.plot(l2i.sea_ice_thickness)
#            plt.show()


#        # Make an orbit plot
#        orbitmap_filename = os.path.join(
#                export_directory, file_basename(l2i_file)+"_map.png")
#        orbitmap = PlotL2iDataOrbit()
#        orbitmap.create_plot(l2i)
#        orbitmap.savefig(orbitmap_filename)
#
#        # make a plot for freeboard and sea surface height anomaly
#        frbssafig_filename = os.path.join(
#                export_directory, file_basename(l2i_file)+"_frbssa.png")
#        frbssafig = PlotL2iSeaSurfaceHeightFreeboard()
#        frbssafig.create_plot(l2i)
#        frbssafig.savefig(frbssafig_filename)

        # make a plot for freeboard and thickness
        frbsitfig_filename = os.path.join(
                export_directory, file_basename(l2i_file)+"_frbsit.png")
        frbsitfig = PlotL2iFreeboardThickness()
        frbsitfig.create_plot(l2i)
        frbsitfig.savefig(frbsitfig_filename)



#        if n_anomalous_frb[i] > 500:
#            print filename_from_path(l2i_file)
#            plt.figure("sea_surface_anomaly")
#            plt.plot(l2i.sea_surface_anomaly)
#            plt.figure("freeboard")
#            plt.plot(l2i.freeboard)
#            plt.show()

#    plt.figure("ssa_var")
#    plt.plot(ssa_var)
#    plt.figure("ssa_mean")
#    plt.plot(ssa_mean)
#    plt.figure("frb_mean")
#    plt.plot(frb_mean)
#    plt.figure("anomalous frb detector")
#    plt.scatter(n_anomalous_frb, n_anomalous_frb_fract)
    plt.figure("frb2sit")
    plt.plot(frb2sit)
    plt.show()


def get_command_line_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            action='store',
            dest='l2i_folder',
            type=str,
            help='l2i folder to inspect')
    return parser.parse_args()

if __name__ == "__main__":
    find_anomalous_frb_l2i()
