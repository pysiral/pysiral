# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 13:57:56 2016

@author: Stefan

Important Note:

    All sitype data handlers must be subclasses of pysiral.auxdata.AuxdataBaseClass in order to work
    for the Level-2 Processor. If the auxiliary class is based on a static dataset, this should be parsed
    in `__init__`.

    Please review the variables and properties in the parent class, as well as the correspodning config and
    support classes for grid track interpolation in the pysiral.auxdata module for additional guidance.

    The only other hard requirements is the presence of on specific method in order to be a valid subclass of
    AuxdataBaseClass:


        get_l2_track_vars(l2)

            This method will be called during the Level-2 processor. The argument is the Level-2 data object and
            the purpose of the method is to compute the auxilary variable(s) and associated uncertainty. These
            variable need to be registered using the `register_auxvar(id, name, value, uncertainty)` method of
            the base class. All sitype subclasses need to register at minimum the following variable:

                sea ice type (fraction of multi year ice):
                    id: sitype
                    name: sea_ice_type

            e.g., this code line is mandatory for `get_l2_track_vars` (uncertainty can be None):

                # Register Variables
                self.register_auxvar("sitype", "sea_ice_type", value, uncertainty)

"""

from pysiral.auxdata import AuxdataBaseClass, GridTrackInterpol
from pysiral.iotools import ReadNC
from pysiral.sla import SLABaseFunctionality

import scipy.ndimage as ndimage
from scipy import interpolate
from pyproj import Proj
import numpy as np
from pathlib import Path


class OsiSafSIType(AuxdataBaseClass):
    """ This is a class for the OSI-403 product with variables ice_type and confidence_level """

    def __init__(self, *args, **kwargs):
        super(OsiSafSIType, self).__init__(*args, **kwargs)
        self._data = None
        self.start_time = None
        self.hemisphere_code = None

    def get_l2_track_vars(self, l2):
        """ Default grid auxiliary data set"""

        # These properties are needed to construct the product path
        self.start_time = l2.info.start_time
        self.hemisphere_code = l2.hemisphere_code

        # Set the requested data
        self.set_requested_date_from_l2(l2)

        # Update the external data
        self.update_external_data()

        # Check if data is available
        if self.error.status:
            self.error.raise_on_error()

        # Get the data
        sitype, uncertainty = self._get_sitype_track(l2)

        # (Optional) Fill gaps in the sea ice type product
        # where the sea-ice concentration data indicates the
        # presence of sea ice. This usually happens close to
        # coastlines since the sea-ice type product has a larger
        # land exclusion zone than SIC (which leads to loss of
        # data, e.g. in the Canadian Archipelago)
        fill_gaps = self.cfg.options.get("fill_valid_sic_gaps", False)
        if fill_gaps:
            sitype, uncertainty = fill_sitype_gaps(sitype, uncertainty, l2.sic[:])

        # Register the data
        self.register_auxvar("sitype", "sea_ice_type", sitype, uncertainty)

    def load_requested_auxdata(self):
        """ Required subclass method: Load the data file necessary to satisfy condition for requested date"""

        # Retrieve the file path for the requested date from a property of the auxdata parent class
        path = Path(self.requested_filepath)

        # Validation
        if not path.is_file():
            msg = "OsiSafSIType: File not found: %s " % path
            self.add_handler_message(msg)
            self.error.add_error("auxdata_missing_sitype", msg)
            return

        # --- Read the data ---
        self._data = ReadNC(path)

        # Report
        self.add_handler_message("OsiSafSIType: Loaded SIType file: %s" % path)

    def _get_sitype_track(self, l2):

        # Extract from grid
        griddef = self.cfg.options[l2.hemisphere]
        grid_lons, grid_lats = self._data.lon, self._data.lat
        grid2track = GridTrackInterpol(l2.track.longitude, l2.track.latitude, grid_lons, grid_lats, griddef)
        sitype = grid2track.get_from_grid_variable(self._data.ice_type[0, :, :], flipud=True)
        confidence_level = grid2track.get_from_grid_variable(self._data.confidence_level[0, :, :], flipud=True)

        # set fill values to flag 0 -> nan
        fillvalues = np.where(sitype == -1)[0]
        sitype[fillvalues] = 0

        # --- Translate sea-ice type codes into myi fraction ---
        # flag_meanings: -1: fill value, 1: open_water, 2: first_year_ice, 3: multi_year_ice, 4: ambiguous
        translator = np.array([np.nan, np.nan, 0.0, 1.0, 0.5])
        sitype = np.array([translator[value] for value in sitype])

        # Translate confidence level into myi fraction uncertainty
        # flag_meaning: 0: unprocessed, 1: erroneous, 2: unreliable, 3: acceptable, 4: good, 5: excellent
        translator = np.array([np.nan, 1., 0.5, 0.2, 0.1, 0.0])
        sitype_uncertainty = np.array([translator[value] for value in confidence_level])

        return sitype, sitype_uncertainty

    @property
    def requested_filepath(self):
        """ Note: this overwrites the property in the super class due to some
        peculiarities with the filenaming (hemisphere code) """
        path = Path(self.cfg.local_repository)
        for subfolder_tag in self.cfg.subfolders:
            subfolder = getattr(self, subfolder_tag)
            path = path / subfolder
        filename = self.cfg.filenaming.format(
            year=self.year, month=self.month, day=self.day,
            hemisphere_code=self.hemisphere_code)
        path = path / filename
        return path


class OsiSafSITypeCDR(AuxdataBaseClass):
    """ Class for reprocessed OSISAF sea ice type products (e.g. for C3S).
    Needs to be merged into single OsiSafSitype class at some point """

    def __init__(self, *args, **kwargs):
        super(OsiSafSITypeCDR, self).__init__(*args, **kwargs)
        self._data = None
        self.start_time = None
        self.hemisphere_code = None

    def get_l2_track_vars(self, l2):
        """
        Mandadory method of AuxdataBaseClass subclass
        :param l2: Level-2 Data object
        :return:
        """

        # These properties are needed to construct the product path
        self.start_time = l2.info.start_time
        self.hemisphere_code = l2.hemisphere_code

        # Set the requested data
        self.set_requested_date_from_l2(l2)

        # Update the external data
        self.update_external_data()

        # Check if error with file I/O
        if self.error.status or self._data is None:
            sitype = self.get_empty_array(l2)
            uncertainty = self.get_empty_array(l2)
        else:
            # Get and return the track
            sitype, uncertainty = self._get_sitype_track(l2)

        # (Optional) Fill gaps in the sea ice type product
        # where the sea-ice concentration data indicates the
        # presence of sea ice. This usually happens close to
        # coastlines since the sea-ice type product has a larger
        # land exclusion zone than SIC (which leads to loss of
        # data, e.g. in the Canadian Archipelago)
        fill_gaps = self.cfg.options.get("fill_valid_sic_gaps", False)
        if fill_gaps:
            sitype, uncertainty = fill_sitype_gaps(sitype, uncertainty, l2.sic[:])

        # Register the data
        self.register_auxvar("sitype", "sea_ice_type", sitype, uncertainty)

    def load_requested_auxdata(self):
        """ Loads file from local repository only if needed """

        # Retrieve the file path for the requested date from a property of the auxdata parent class
        path = Path(self.requested_filepath)

        # Validation
        if not path.is_file():
            msg = "%s: File not found: %s " % (self.__class__.__name__, path)
            self.add_handler_message(msg)
            self.error.add_error("auxdata_missing_sitype", msg)
            return

        # Read and prepare input data
        self._data = ReadNC(path)

    def _get_sitype_track(self, l2):
        """ Extract ice type and ice type uncertainty along the track """

        # Extract from grid
        griddef = self.cfg.options[l2.hemisphere]
        grid_lons, grid_lats = self._data.lon, self._data.lat
        grid2track = GridTrackInterpol(l2.track.longitude, l2.track.latitude, grid_lons, grid_lats, griddef)
        sitype = grid2track.get_from_grid_variable(self._data.ice_type[0, :, :], flipud=True)
        uncertainty = grid2track.get_from_grid_variable(self._data.uncertainty[0, :, :], flipud=True)

        # Convert flags to myi fraction
        translator = np.array([np.nan, np.nan, 0.0, 1.0, 0.5, np.nan])
        fillvalues = np.where(sitype == -1)[0]
        sitype[fillvalues] = 5
        sitype = np.array([translator[value] for value in sitype])

        # Uncertainty in product is in %
        sitype_uncertainty = uncertainty / 100.

        return sitype, sitype_uncertainty

    @property
    def requested_filepath(self):
        """ Note: this overwrites the property in the super class due to some
        peculiarities with the filenaming (auto product changes etc) """

        # Unique to this class is the possibility to auto merge
        # products. The current implementation supports only two products
        path = Path(self.cfg.local_repository)

        # The path needs to be completed if two products shall be used
        if "auto_product_change" in self.cfg.options:
            opt = self.cfg.options.auto_product_change
            product_index = int(self.start_time > opt.date_product_change)
            product_def = opt.osisaf_product_def[product_index]
            path = path / product_def["subfolder"]
            self.cfg.filenaming = product_def["filenaming"]
            self.cfg.long_name = product_def["long_name"]

        for subfolder_tag in self.cfg.subfolders:
            subfolder = getattr(self, subfolder_tag)
            path = path / subfolder

        filename = self.cfg.filenaming.format(
            year=self.year, month=self.month, day=self.day,
            hemisphere_code=self.hemisphere_code)
        path = path / filename
        return path


class ICDCNasaTeam(AuxdataBaseClass):
    """ MYI Fraction from NASA Team Algorithm (from ICDC UHH) """

    def __init__(self, *args, **kwargs):
        super(ICDCNasaTeam, self).__init__(*args, **kwargs)
        self._data = None

    def get_l2_track_vars(self, l2):
        self._get_requested_date(l2)
        self._get_data(l2)
        if self.error.status:
            sitype, sitype_uncertainty = self.get_empty_array(l2), self.get_empty_array(l2)
        else:
            sitype, sitype_uncertainty = self._get_sitype_track(l2)
        self.register_auxvar("sitype", "sea_ice_type", sitype, None)

    def _get_requested_date(self, l2):
        """ Use first timestamp as reference, date changes are ignored """
        year = l2.track.timestamp[0].year
        month = l2.track.timestamp[0].month
        day = l2.track.timestamp[0].day
        self._requested_date = [year, month, day]

    def _get_data(self, l2):
        """ Loads file from local repository only if needed """

        opt = self.cfg.options

        # Check if file is already loaded
        if self._requested_date == self._current_date:
            # Data already loaded, nothing to do
            self.add_handler_message("ICDCNasaTeam: Daily grid already present")
            return

        # construct filename
        path = Path(self._get_local_repository_filename(l2))

        # Check if the file exists, add an error if not
        # (error is not raised at this point)
        if not path.isfile():
            msg = "ICDCNasaTeam: File not found: %s " % path
            self.add_handler_message(msg)
            self.error.add_error("auxdata_missing_sitype", msg)
            return

        # Bulk read the netcdf file
        self._data = ReadNC(path)

        # There are multiple myi concentrations fields in the product
        # The one used here is defined in the auxdata definition file
        # in the pysiral config folder (`auxdata_def.yaml`)
        # -> root.sitype.icdc_nasateam.options.variable_name
        myi_fraction = getattr(self._data, opt.variable_name)
        self._data.ice_type = myi_fraction[0, :, :]

        # Same for the uncertainty variable
        # (see description directly above for how to access variable name
        #  definition)
        myi_fraction_unc = getattr(self._data, opt.uncertainty_variable_name)
        self._data.ice_type_uncertainty = myi_fraction_unc[0, :, :]

        # Report and save current data period
        self.add_handler_message("ICDCNasaTeam: Loaded SIType file: %s" % path)
        self._current_date = self._requested_date

    def _get_local_repository_filename(self, l2):
        path = Path(self.cfg.local_repository)
        for subfolder_tag in self.cfg.subfolders:
            subfolder = getattr(self, subfolder_tag)
            path = path / subfolder
        filename = self.cfg.filenaming.format(
            year=self.year, month=self.month, day=self.day,
            hemisphere_code=l2.hemisphere_code)
        path = path / filename
        return path

    def _get_sitype_track(self, l2):

        # Convert grid/track coordinates to grid projection coordinates
        kwargs = self.cfg.option[l2.hemisphere].projection
        p = Proj(**kwargs)
        x, y = p(self._data.longitude, self._data.latitude)
        l2x, l2y = p(l2.track.longitude, l2.track.latitude)

        # Convert track projection coordinates to image coordinates
        # x: 0 < n_lines; y: 0 < n_cols
        dim = self.cfg.option[l2.hemisphere].dimension

        x_min = x[0, 0]-(0.5*dim.dx)
        y_min = y[0, 0]-(0.5*dim.dy)
        ix, iy = (l2x-x_min)/dim.dx, (l2y-y_min)/dim.dy

        # Extract along track data from grid
        myi_concentration_percent = ndimage.map_coordinates(self._data.ice_type, [iy, ix], order=0)

        myi_concentration_uncertainty = ndimage.map_coordinates(self._data.ice_type_uncertainty, [iy, ix], order=0)

        # Convert percent [0-100] into fraction [0-1]
        sitype = myi_concentration_percent/100.
        sitype_uncertainty = myi_concentration_uncertainty/100.

        # Remove invalid (negative) values
        sitype[np.where(sitype < 0)] = 0.0
        sitype_uncertainty[np.where(sitype_uncertainty < 0)] = 0.0

        return sitype, sitype_uncertainty


class MYIDefault(AuxdataBaseClass):
    """ Returns myi for all ice covered regions """

    def __init__(self, *args, **kwargs):
        super(MYIDefault, self).__init__(*args, **kwargs)

    def get_l2_track_vars(self, l2):
        """ Every ice is myi (sitype = 1) """
        sitype = np.zeros(shape=l2.sic.shape, dtype=np.float32)
        is_ice = np.where(l2.sic > 0)[0]
        sitype[is_ice] = 1.0
        uncertainty = np.full(sitype.shape, self.uncertainty_default)
        # Register the data
        self.register_auxvar("sitype", "sea_ice_type", sitype, uncertainty)

    @property
    def uncertainty_default(self):
        if "uncertainty_default" in self.cfg.option:
            return self.cfg.option.uncertainty_default
        else:
            return 0.0


class FYIDefault(AuxdataBaseClass):
    """ Returns myi for all ice covered regions """

    def __init__(self, *args, **kwargs):
        super(FYIDefault, self).__init__(*args, **kwargs)

    def get_l2_track_vars(self, l2):
        """ Every ice is fyi (sitype = 0) """
        sitype = np.zeros(shape=l2.sic.shape, dtype=np.float32)
        uncertainty = np.full(sitype.shape, self.uncertainty_default)
        self.register_auxvar("sitype", "sea_ice_type", sitype, uncertainty)

    @property
    def uncertainty_default(self):
        try:
            if "uncertainty_default" in self.cfg.options:
                return self.cfg.options.uncertainty_default
        except TypeError:
            pass
        return 0.0


def fill_sitype_gaps(sitype, sitype_uncertainty, sic, sic_threshold=70., gap_filled_uncertainty=0.5,
                     max_valid_nn_dist=150, ambiguos_sitype_value=0.5):
    """
    Fill gaps in sea ice type that are indicated as sea ice in the sea ice concentration product.
    Gaps close to a valid sea-ice type value will be filled using a nearest neighbour approach,
    further gaps will be labeled as ambiguous. Any gap-filled value will be associated with a
    large uncertainty.
    :param sitype: sea-ice type array
    :param sitype_uncertainty: sea-ice type uncertainty array
    :param sic: sea-ice concentration array
    :param sic_threshold: threshold in percent as to when a waveform can be intepreted as sea ice
    :param gap_filled_uncertainty: the sea-ice type uncertainty value associated with gap filled values
    :param max_valid_nn_dist: The maximum distance (in points) for valid nearest neighbour interpolation
    :param ambiguos_sitype_value: The value associated to ambiguous gaps
    :return: Gap-filled sitype, sitype_uncertainty
    """

    sitype_0 = np.copy(sitype)

    # Step 1: Find gaps (if any)
    is_seaice = sic >= sic_threshold
    is_sitype_gap = np.isnan(sitype)
    is_valid_sitype = np.logical_not(is_sitype_gap)
    fillable_gap = np.logical_and(is_seaice, is_sitype_gap)
    all_gap_indices = np.where(fillable_gap)[0]

    # Check if there is anything to do
    if len(all_gap_indices) == 0:
        return sitype, sitype_uncertainty

    # Step 2: Compute the distance of a sea ice type gap to the next valid value
    # Must have at at least some valid sea ice type values
    if not is_sitype_gap.all():
        gap_dist = SLABaseFunctionality().get_tiepoint_distance(is_valid_sitype)

        # Step 3: Fill closer gaps with nearest neighbour approach
        fillable_gap_nn = np.logical_and(fillable_gap, gap_dist <= max_valid_nn_dist)
        nn_gap_indices = np.where(fillable_gap_nn)[0]
        if len(nn_gap_indices) > 0:
            x = np.arange(len(sitype))
            nn_interp = interpolate.interp1d(x[is_valid_sitype], sitype[is_valid_sitype],
                                             kind="nearest", fill_value="extrapolate",
                                             bounds_error=False)
            sitype[nn_gap_indices] = nn_interp(x[nn_gap_indices])

    # If no valid sea-ice type information is available than all
    # gap values exceed the maximum nearest neighbour distance
    else:
        gap_dist = np.full(sitype.shape, max_valid_nn_dist + 1)

    # Step 4: Fill rest of the gaps with ambiguous
    ambiguos_gaps = np.logical_and(fillable_gap, gap_dist > max_valid_nn_dist)
    ambiguos_gap_indices = np.where(ambiguos_gaps)[0]
    sitype[ambiguos_gap_indices] = ambiguos_sitype_value

    # Step 5: Set uncertainty of all gaps
    sitype_uncertainty[all_gap_indices] = gap_filled_uncertainty

    # import matplotlib.pyplot as plt
    # plt.figure()
    # plt.plot(sitype_0, lw=2, label="sea ice type (orig)")
    # plt.plot(sitype, lw=1, label="sea ice type (gap filled)")
    # plt.plot(sic/100, label="sea ice concentration")
    # plt.legend()
    #
    # plt.figure()
    # plt.plot(gap_dist)
    #
    # plt.show()
    # breakpoint()

    # Return
    return sitype, sitype_uncertainty
