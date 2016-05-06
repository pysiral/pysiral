# -*- coding: utf-8 -*-

from pysiral.iotools import get_l1bdata_files
import argparse
import os


def create_monthly_orbit_plot():

    # Get options from command line arguments
    parser = get_argparser()
    args = parser.parse_args()

    # Get the list of l2 files
    l1b_file_list = get_l1bdata_files(
        args.mission_id,
        args.hemisphere,
        args.date[0],
        args.date[1])

    # Creating the export filename
    map_filename = os.path.join(
        args.destination_folder,
        "orbit_coverage_{mission_id}_{hemisphere}_{year}{month}.png".format(
            mission_id=args.mission_id,
            hemisphere=args.hemisphere,
            year="%04g" % args.date[0],
            month="%02g" % args.date[1]))

    orbitmap = MonthlyOrbitMap()
    orbitmap.hemisphere = args.hemisphere
    orbitmap.mission_id = args.mission_id
    orbitmap.filename = map_filename
    orbitmap.l1b_file_list = l1b_file_list
    orbitmap.create()


def get_argparser():
    """ Handle command line arguments """

    parser = argparse.ArgumentParser()

    # Mission id string: cryosat2, envisat, ...
    parser.add_argument(
        '-m', required=True,
        action='store', dest='mission_id',
        help='[cryosat2|envisat|...]')

    parser.add_argument(
        '-hemisphere', required=True,
        action='store', dest='hemisphere',
        help='[north|south]')

    parser.add_argument(
        '-d', required=True,
        action='store', dest='destination_folder',
        help='')

    # Start month as list: [yyyy, mm]
    parser.add_argument(
        '-month', required=True, type=int,
        action='store', dest='date', nargs='+',
        help='month to plot given as year and month (-month yyyy mm)')

    return parser


class MonthlyOrbitMap(object):

    projection = {
        "north": {
            "projection": "ortho",
            "lon_0": -45,
            "lat_0": 80,
            "resolution": "i"},
        "south": {
            "projection": "ortho",
            "lon_0": 0,
            "lat_0": -70,
            "resolution": "i"}
            }
    color = {
        }

    def __init__(self):
        self.hemisphere = None
        self.filename = None
        self.l1b_file_list = []
        self._figure = None

    def create(self):

        from pysiral.l1bdata import L1bdataNCFile
        from pysiral.iotools import get_temp_png_filename

        from mpl_toolkits.basemap import Basemap
        import matplotlib.pyplot as plt
        import numpy as np

        # Calculate shading and save as temp file for warpimage
        lats = np.linspace(-0.5*np.pi, 0.5*np.pi, num=500)

        if self.hemisphere == "south":
            shading = np.zeros(shape=(500, 1000))
            indices = np.where(lats > np.deg2rad(-88))[0]
            colormap_name = "Blues"
        elif self.hemisphere == "north":
            shading = np.ones(shape=(500, 1000)) * 2.0
            indices = np.where(lats < np.deg2rad(88))[0]
            colormap_name = "Blues_r"
        for i in indices:
            shading[i, :] = np.sin(lats[i])+1.0

        temp_filename = get_temp_png_filename()
        fig = plt.figure(figsize=(10, 5))
        ax = plt.axes([0, 0, 1, 1])
        plt.axis('off')
        ax.imshow(np.flipud(shading), cmap=plt.get_cmap(colormap_name))
        plt.savefig(temp_filename)
        plt.close(fig)

        # Create the actual map
        grid_keyw = {"dashes": (None, None), "color": "#91e5ff",
                     "linewidth": 0.2, "latmax": 88, "zorder": 10}

        plt.figure(figsize=(12, 12), facecolor="#ffffff")
        m = Basemap(**self.projection[self.hemisphere])
        m.fillcontinents(color='#91e5ff', lake_color='#91e5ff',
                         zorder=90)
        m.drawmapboundary(color="#91e5ff", linewidth=0.2)

        # draw parallels and meridians.
        m.drawparallels(np.arange(-90., 120., 10.), **grid_keyw)
        m.drawmeridians(np.arange(0., 420., 30.), **grid_keyw)

        for i, l1b_file in enumerate(self.l1b_file_list):
            l1b = L1bdataNCFile(l1b_file)
            l1b.parse()
            x, y = m(l1b.time_orbit.longitude, l1b.time_orbit.latitude)
            m.plot(x, y, color="#FF1700", linewidth=0.5,
                   zorder=300, alpha=0.75)

        m.warpimage(temp_filename, zorder=200, alpha=0.5)

        plt.savefig(self.filename, bbox_inches="tight")

        os.remove(temp_filename)

if __name__ == "__main__":
    create_monthly_orbit_plot()
