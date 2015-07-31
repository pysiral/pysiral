# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 14:09:32 2015

@author: Stefan
"""

from pysiral.config import ConfigInfo
from pysiral.job import Level2Job
from pysiral.l2proc import Level2Processor

import os
import glob
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl


def l2_processing():

    # Get configuration
    config = ConfigInfo()

    # Get an L1B SAR file
    l1b_directory = config.local_machine.local_l1b_repository.cryosat2.sar
    l1b_directory = os.path.join(l1b_directory, "2015", "04")
    l1b_files = glob.glob(os.path.join(l1b_directory, "*.DBL"))

    # Simulate the jobconfig
    # This has to come from the job configuration file
    mission_settings = {
        "id": "cryosat2",
        "options": config.get_mission_defaults("cryosat2")}
    roi_settings = {
        "pyclass": "LowerLatLimit",
        "options": {
            "latitude_threshold": -50.0}}
    l2_settings = {
        "corrections": [
            "dry_troposphere",
            "wet_troposphere",
            "inverse_barometric",
            "dynamic_atmosphere",
            "ionospheric",
            "ocean_tide_elastic",
            "ocean_tide_long_period",
            "ocean_loading_tide",
            "solid_earth_tide",
            "geocentric_polar_tide"],
        "surface_type": {
            "pyclass": "RickerTC2014",
            "options": {
                "ocean": {
                    "peakiness_min": 0.0,
                    "peakiness_max": 10.0,
                    "stack_standard_deviation_min": 18.5,
                    "ice_concentration_min": 5.0,
                    "ocog_width_min": 38},
                "lead": {
                    "peakiness_l_min": 40.0,
                    "peakiness_r_min": 30.0,
                    "peakiness_min": 40.0,
                    "stack_kurtosis_min": 40.0,
                    "stack_standard_deviation_max": 4.0,
                    "ice_concentration_min": 70.0},
                "sea_ice": {
                    "peakiness_r_max": 15.0,
                    "peakiness_l_max": 20.0,
                    "peakiness_max": 30.0,
                    "stack_kurtosis_max": 8.0,
                    "ice_concentration_min": 70.0}}},
        "retracker": {
            "ocean": {
                "pyclass": "TFMRA",
                "options": {
                    "threshold": 0.5,
                    "offset": 0.0,
                    "wfm_oversampling_factor": 10,
                    "wfm_oversampling_method": "linear",
                    "wfm_smoothing_window_size": 11,
                    "first_maximum_normalized_threshold": 0.15,
                    "first_maximum_local_order": 1}},
            "lead": {
                "pyclass": "TFMRA",
                "options": {
                    "threshold": 0.5,
                    "offset": 0.0,
                    "wfm_oversampling_factor": 10,
                    "wfm_oversampling_method": "linear",
                    "wfm_smoothing_window_size": 11,
                    "first_maximum_normalized_threshold": 0.15,
                    "first_maximum_local_order": 1}},
            "sea_ice": {
                "pyclass": "TFMRA",
                "options": {
                    "threshold": 0.5,
                    "offset": 0.0,
                    "wfm_oversampling_factor": 10,
                    "wfm_oversampling_method": "linear",
                    "wfm_smoothing_window_size": 11,
                    "first_maximum_normalized_threshold": 0.15,
                    "first_maximum_local_order": 1}}},
        "mss": {},
        "filter": {},
        "post_processing": {},
        "output": {}}

    # Assemble the job order
    job = Level2Job()
    job.mission_settings(mission_settings)
    job.roi_settings(roi_settings)
    job.l2proc_settings(l2_settings)
    job.validate()

    # Start the processor
    l2proc = Level2Processor(job)
    l2proc.set_config(config)
    l2proc.error_handling(raise_on_error=True)
    l2proc.set_l1b_files(l1b_files[0:1])
    l2proc.run()

    # Test Plots
    create_surface_type_plot(l2proc.orbit[0])


def create_surface_type_plot(l2):

    from matplotlib.colors import ListedColormap

    # AWI eisblau #00ace5
    # AWI tiefblau #003e6e
    # AWI grau 1 #4b4b4d
    # AWI grau 2 #bcbdbf

    rgb_list = [
        "#bcbdbf",  # Unkown
        "#003e6e",  # Ocean
        "#00ace5",  # Lead
        "#DDA0DD",  # Polynya
        "#76FF7A",  # Sea Ice
        "#4b4b4d",  # Lakes
        "#FFCC00",  # Land Ice
        "#826644"   # Land
        ]
    cmap = ListedColormap(rgb_list, name='surface_type')

    # Calculate Fractions
    n = float(l2.surface_type.n_records)
    flag = l2.surface_type.flag
    fractions = []
    labels = [l2.surface_type.name(i) for i in np.arange(8)]
    for i in np.arange(8):
        num = len(np.where(flag == i)[0])
        fractions.append(float(num)/n)

    # Pop fraction = 0
    pie_fractions = np.array(fractions)
    pie_labels = np.array(labels)
    pie_colors = np.array(rgb_list)

    non_zero = np.nonzero(pie_fractions)[0]
    pie_fractions = pie_fractions[non_zero]
    pie_labels = pie_labels[non_zero]
    pie_colors = pie_colors[non_zero]

    plt.figure(facecolor="white")

    plt.subplot2grid((5, 15), (0, 0), colspan=14, rowspan=4)
    plt.gca().set_aspect(1.0)
    wedges, texts, autotexts = plt.pie(
        pie_fractions, labels=pie_labels, colors=pie_colors,
        autopct='%1.1f%%', startangle=90, labeldistance=1.1, pctdistance=0.75)
    for wedge in wedges:
        wedge.set_width(0.5)
        wedge.set_ec("none")

    plt.subplot2grid((5, 15), (4, 0), colspan=14)
    ax = plt.gca()
    im = plt.imshow([flag], aspect=100, interpolation=None,
                    vmin=-0.5, vmax=7.5, cmap=cmap)
    ax.xaxis.set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xlabel("Record #")
    ax.yaxis.set_ticks([])
    ax.spines["bottom"].set_position(("data", 0.75))

    spines_to_remove = ["top", "right", "left"]
    for spine in spines_to_remove:
        ax.spines[spine].set_visible(False)

    plt.subplot2grid((5, 15), (0, 14), rowspan=4)
    cbar = plt.colorbar(im, orientation="vertical", cax=plt.gca())
    cbar_ticks = np.arange(8)
    cbar.set_ticks(cbar_ticks)
    cbar.set_ticklabels(labels)
    cbar.solids.set_edgecolor("1.0")
    cbar.outline.set_linewidth(5)
    cbar.outline.set_alpha(0.0)
    cbar.ax.tick_params('both', length=0.1, which='major', pad=5)

    plt.tight_layout()
    plt.show()


if __name__ == '__main__':

    mpl.rcParams['font.sans-serif'] = "arial"
    for target in ["xtick.color", "ytick.color", "axes.edgecolor",
                   "axes.labelcolor"]:
        mpl.rcParams[target] = "#4b4b4d"

    l2_processing()
