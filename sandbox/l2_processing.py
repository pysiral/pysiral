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
    l1b_directory = config.local_machine.l1b_repository.cryosat2.l1bdata
    l1b_directory = os.path.join(l1b_directory, "north", "2015", "03")
    l1b_files = glob.glob(os.path.join(l1b_directory, "*.nc"))

    # Simulate the jobconfig
    # This has to come from the job configuration file
    mission_settings = {
        "id": "cryosat2",
        "options": config.get_mission_defaults("cryosat2")}
    roi_settings = {
        "pyclass": "LowerLatLimit",
        "hemisphere": "north",
        "options": {
            "latitude_threshold": 50.0}}
    l2_settings = {
        "auxdata": {
            "mss": {
                "name": "dtu15",
                "options": {}},
            "sic": {
                "name": "osisaf",
                "options": {}},
            "sitype": {
                "name": "osisaf",
                "options": {}},
            "snow": {
                "name": "warren99",
                "options": {
                    "fyi_correction_factor": 0.5,
                    "smooth_snow_depth": True,
                    "smooth_filter_width_m": 25000.0,
                    "valid_snow_depth_range": [0, 0.6],
                    "invalid_replace": "nan"}}},
        "corrections": [
            "dry_troposphere",
            "wet_troposphere",
            "inverse_barometric",
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
        "ssa": {
            "pyclass": "SSASmoothedLinear",
            "options": {
                "critical_number_of_tiepoints": 5,
                "use_ocean_wfm": False,
                "smooth_filter_width_m": 25000.0,
                "smooth_filter_width_footprint_size": 300.0}},
        "frb": {
            "pyclass": "SnowGeometricCorrection",
            "options": {
                "vacuum_light_speed_reduction": 0.22}},
        "filter": {
            "freeboard": {
                "frb_minmax": {
                    "pyclass": "FreeboardValidRange",
                    "options": {
                        "valid_minimum_point_value": -0.5,
                        "valid_maximum_point_value": 2.5}}}},
        "validator": {
            "surface_type": {
                "n_leads": {
                    "pyclass": "LeadFixedMinimumNumber",
                    "options": {"minimum_n_leads": 3}},
                }},
        "post_processing": {},
        "output": {}}

    # Assemble the job order
    job = Level2Job()
    job.local_machine_settings(config.local_machine)
    job.mission_settings(mission_settings)
    job.roi_settings(roi_settings)
    job.l2proc_settings(l2_settings)
    job.validate()

    # Start the processor
    l2proc = Level2Processor(job)
    l2proc.set_config(config)
    l2proc.error_handling(raise_on_error=True)
    l2proc.set_l1b_files(l1b_files[0:5])
    l2proc.run()

#    with open(r"D:\pysiral_tfmra_test.dat", "w") as fh:
#        for elevation in l2proc.orbit[0].elev:
#            fh.write("%14.5f\n" % elevation)

    # Test Plots
    create_orbit_map(l2proc.orbit[1], block=False)
    create_surface_type_plot(l2proc.orbit[1], block=False)
    create_surface_elevation_plot(l2proc.orbit[1])


def create_surface_type_plot(l2, block=True):

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
        "#826644",  # Land
        "#FF1700"   # invalid
        ]
    cmap = ListedColormap(rgb_list, name='surface_type')

    # Calculate Fractions
    n = float(l2.surface_type.n_records)
    flag = l2.surface_type.flag
    fractions = []
    labels = [l2.surface_type.name(i) for i in np.arange(9)]
    for i in np.arange(9):
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
    im = plt.imshow([flag], aspect=500, interpolation=None,
                    vmin=-0.5, vmax=8.5, cmap=cmap)
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
    cbar_ticks = np.arange(9)
    cbar.set_ticks(cbar_ticks)
    cbar.set_ticklabels(labels)
    cbar.solids.set_edgecolor("1.0")
    cbar.outline.set_linewidth(5)
    cbar.outline.set_alpha(0.0)
    cbar.ax.tick_params('both', length=0.1, which='major', pad=5)

    plt.tight_layout()
    plt.show(block=block)


def create_surface_elevation_plot(l2, block=True):

    plot_style = {
        "ocean": {
            "color": "#003e6e",
            "sym": "o",
            "label": "Ocean Waveforms"},
        "lead": {
            "color": "#00ace5",
            "sym": "s",
            "label": "Lead Waveforms"},
        "sea_ice": {
            "color": "#76FF7A",
            "sym": "h",
            "label": "Sea Ice Waveforms"}}

    plt.figure(facecolor="white")
    # plot elevations
    plt.subplot(211)
    ax = plt.gca()
    x = np.arange(l2.n_records)
    for surface_type_name in "ocean", "lead", "sea_ice":
        type_definition = l2.surface_type.get_by_name(surface_type_name)
        if len(type_definition.indices) == 0:
            continue
        plt.scatter(
            x[type_definition.indices], l2.elev[type_definition.indices],
            color=plot_style[surface_type_name]["color"],
            marker=plot_style[surface_type_name]["sym"])
    plt.plot(x, l2.mss, color="black", lw=2, label="Mean Sea Surface (DTU15)")
    plt.plot(x, l2.mss+l2.ssa, color="violet", lw=2,
             label="Sea Surface Anomaly")
    ssh_tiepoints = l2.surface_type.lead.indices
    for ssh_tiepoint in ssh_tiepoints:
        y = [l2.mss[ssh_tiepoint], l2.elev[ssh_tiepoint]]
        plt.vlines(x[ssh_tiepoint], min(y), max(y))

    spines_to_remove = ["top", "right"]
    for spine in spines_to_remove:
        ax.spines[spine].set_visible(False)

    ax.yaxis.set_tick_params(direction='out')
    ax.yaxis.set_ticks_position('left')
    ax.set_ylabel("Elevation (m)")
    ax.spines["left"].set_position(("axes", -0.01))

    ax.xaxis.set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines["bottom"].set_position(("axes", -0.01))
    ax.set_xlim([0, l2.n_records])
    plt.legend(frameon=False, loc=1, ncol=2, fontsize="medium")

    # plot freeboard
    plt.subplot(212)
    ax = plt.gca()
    plt.hlines(0.0, 0, l2.n_records)
    for surface_type_name in "ocean", "lead", "sea_ice":
        type_definition = l2.surface_type.get_by_name(surface_type_name)
        if len(type_definition.indices) == 0:
            continue
        plt.scatter(
            x[type_definition.indices], l2.rfrb[type_definition.indices],
            color=plot_style[surface_type_name]["color"],
            marker=plot_style[surface_type_name]["sym"],
            label=plot_style[surface_type_name]["label"])

    plt.legend(frameon=False, loc=1, ncol=2, fontsize="medium")

    spines_to_remove = ["top", "right"]
    for spine in spines_to_remove:
        ax.spines[spine].set_visible(False)

    ax.yaxis.set_tick_params(direction='out')
    ax.yaxis.set_ticks_position('left')
    ax.set_ylabel("Radar Freeboard (m)")
    ax.set_xlabel("Record #")
    ax.spines["left"].set_position(("axes", -0.01))

    ax.xaxis.set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines["bottom"].set_position(("axes", -0.01))
    # ax.set_ylim([-0.5, 2])
    ax.set_xlim([0, l2.n_records])
    plt.show(block=block)


def create_orbit_map(l2, block=True):

    from mpl_toolkits.basemap import Basemap
    import matplotlib.pyplot as plt

    grid_keyw = {"dashes": (None, None), "color": "#bcbdbf",
                 "linewidth": 0.5, "latmax": 88}

    gridb_keyw = {"dashes": (None, None), "color": "#003e6e",
                  "linewidth": 2, "latmax": 88}

    lon_0 = l2.info.lon_max
    lat_0 = l2.info.lat_max

    plt.figure("l2 Debug Map", facecolor="#ffffff")
    m = Basemap(projection='ortho', lon_0=lon_0, lat_0=lat_0, resolution='l')
    m.fillcontinents(color='#00ace5', lake_color='#00ace5')

    # draw parallels and meridians.
    m.drawparallels(np.arange(-90., 120., 10.), **grid_keyw)
    m.drawparallels(np.arange(-50., 51., 100.), **gridb_keyw)
    m.drawmeridians(np.arange(0., 420., 30.), **grid_keyw)

    x, y = m(l2.track.longitude, l2.track.latitude)
    m.plot(x, y, color="#003e6e", linewidth=2.0, zorder=100)
    m.scatter(x[0], y[0], s=50, color="#003e6e", zorder=100)

    plt.show(block=block)

if __name__ == '__main__':

    mpl.rcParams['font.sans-serif'] = "arial"
    for target in ["xtick.color", "ytick.color", "axes.edgecolor",
                   "axes.labelcolor"]:
        mpl.rcParams[target] = "#4b4b4d"

    l2_processing()
