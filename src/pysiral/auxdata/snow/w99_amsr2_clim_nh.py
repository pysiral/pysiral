# -*- coding: utf-8 -*-

import re
from datetime import datetime
from pathlib import Path

import numpy as np
from loguru import logger

from xarray import open_dataset

from pysiral.auxdata import AuxdataBaseClass, GridTrackInterpol
from pysiral.auxdata.snow._data import SnowParameterContainer
from pysiral.core.legacy_classes import ErrorStatus


class Warren99AMSR2Clim(AuxdataBaseClass):
    """
    Class for monthly snow depth & density climatology based on merged Warren99 climatology and
    monthly AMSR2 snow depth composite (source: IUP). The source data is organized as
    netCDF files that contain a:

        1. monthly climatological snow depth and density
        2. the mask of the W99 region of influence
        3. uncertainties of geophysical parameters

    This class reads the netCDF data and applies a correction for first-year sea ice areas in regions
    solely influenced by the Warren99 climatalogy.

    REQUIREMENTS:

        - Level-2 object has attribute `sitype`
          (sea ice type classication: 0 [fyi] <) sitype <= 1 [myi])
        - options has attribute `fyi_correction_factor`
          (factor 0-1 for reduction of snow depth in FYI W99 areas)

    OPTIONAL:

        - options has boolean attribute `daily_scaling`
          (if True snow depth and density will be scaled between successive month)

    NOTES:

        - This class is mainly designed to be calles from the Level-2 processor for trajectory based data sets.
          An more generalized version is in planning.

    UPDATES:

        - [July 2020] With AWI CryoSat-2 v2.3 a change was introduced to allow a daily change of the snow depth fields
          rather than monthly (see github issue: https://github.com/shendric/pysiral/issues/40)
    """

    def __init__(self, *args, **kwargs):
        """
        Init the class for getting snow depth on sea ice with density and respective uncertainties.
        :param args: Arguments for AuxdataBaseClass
        :param kwargs: Keyword arguments for AuxdataBaseClass
        """
        super(Warren99AMSR2Clim, self).__init__(*args, **kwargs)
        self._data = Warren99AMSR2ClimDataContainer(self.cfg, self.use_daily_scaling)

    def get_l2_track_vars(self, l2):
        """
        This is the method that will be evoked by the Level-2 processor. This method will extract the
        geophysical parameters

        :param l2: The Level-2 data container
        :return: None
        """

        # Set the requested date
        self.set_requested_date_from_l2(l2)

        # Update the external data
        # NOTE: This will only be done once as the climatology has the same period as the Level-2 processor
        self.update_external_data()

        # Check if data has not been loaded
        if not self.has_data_loaded:
            logger.error(f"- Source data cannot be loaded {self.requested_filepath}")

        # Check if requested date is within validity period of the
        # climatology (winter month of October -> April only)
        valid_month = True
        if int(self.month) not in [10, 11, 12, 1, 2, 3, 4]:
            valid_month = False
            logger.warning("- Target month is outside snow climatology coverage (October - April)")

        # Ensure that snow depth and uncertainty variables are added (even if just NaN)
        if not valid_month or not self.has_data_loaded:
            logger.warning("- Previous errors and/or warnings, snow depth -> NaN")
            snow = SnowParameterContainer()
            snow.set_dummy(l2.n_records)
        else:
            # Extract along track snow depth and density
            snow = self._get_snow_track(l2)

        # Register Variables
        self.register_auxvar("sd", "snow_depth", snow.depth, snow.depth_uncertainty)
        self.register_auxvar("sdens", "snow_density", snow.density, snow.density_uncertainty)

    def update_external_data(self):
        """
        This method overwrites the default behaviour of loading a dedicated file per
        l2 track object. The functionality is delegated to a specific data class for this
        auxiliary data set.
        :return:
        """
        if not self._data.has_data_loaded:
            self._data.load()
            if self._data.has_data_loaded:
                for filepath in self._data.filepaths:
                    self.add_handler_message(self.__class__.__name__ + ": Load {}".format(filepath))
            else:
                msg = ": Loading data has failed failed"
                self.add_handler_message(self.__class__.__name__ + msg)
        else:
            if self._data.has_data_loaded:
                self.add_handler_message(self.__class__.__name__+": Data already present")
            else:
                msg = ": No Data: Loading failed in an earlier attempt"
                self.add_handler_message(self.__class__.__name__ + msg)

    def _get_snow_track(self, l2):
        """
        Extract the snow depth and density track along the l2 track
        :param l2:
        :return:
        """

        # Extract track data from grid
        griddef = self.cfg.options[l2.hemisphere]
        grid_lons, grid_lats = self._data.get_lonlat()
        grid2track = GridTrackInterpol(l2.track.longitude, l2.track.latitude, grid_lons, grid_lats, griddef)

        # Extract data (Map the extracted tracks directly on the snow parameter container)
        var_map = self.cfg.options.variable_map
        snow = SnowParameterContainer()
        for var_name in var_map.keys():
            source_name = var_map[var_name]
            sdgrid = self._data.get_var(source_name, self._requested_date)
            setattr(snow, var_name, grid2track.get_from_grid_variable(sdgrid))

        # Extract the W99 weight for the specific track
        w99_weight = grid2track.get_from_grid_variable(self._data.w99_weight)

        # Apply the same modification as the Warren climatology
        # Apply ice_type but this time modified by the regional weight of the Warren climatology.
        # The weight ranges from 0 to 1 and ensures that snow reduction on first-year sea ice is
        # only applied to the W99 contribution of the snow depth climatology.
        # NOTE: sea ice type here is defined as MYI fraction: 0 [fyi] <= sitype <= 1 [myi]
        scale_factor = (1.0 - l2.sitype) * self.cfg.options.fyi_correction_factor * w99_weight

        # The scaling factor affects the snow depth ...
        snow.depth = snow.depth - scale_factor * snow.depth

        # ... and the uncertainty. Here it is assumed that the uncertainty
        # is similar affected by the scaling factor.
        snow.depth_uncertainty = snow.depth_uncertainty - scale_factor * snow.depth_uncertainty

        # The uncertainty of the MYI fraction is acknowledged by adding
        # a term that depends on snow depth, the magnitude of scaling
        # and the sea ice type uncertainty
        scaling_uncertainty = snow.depth * scale_factor * l2.sitype.uncertainty * w99_weight
        snow.depth_uncertainty = snow.depth_uncertainty + scaling_uncertainty

        return snow

    @property
    def use_daily_scaling(self):
        """
        Return a flag that indicates whether to use daily scaling (absence of flag in options will be treated as no)
        :return:
        """
        if "daily_scaling" not in self.cfg.options:
            return False
        else:
            return self.cfg.options.daily_scaling


class Warren99AMSR2ClimDataContainer(object):
    """
    A dedicated data container for the merged W99/AMSR2 snow climatology. This class has been introduced
    with the use of daily scaling that requires data to loaded also from month adjacent to the month
    of the current Level-2 data object
    """

    def __init__(self, cfg, use_daily_scaling):
        """
        Init the class
        :param cfg: A copy of the auxdata class configuration
        :param use_daily_scaling:
        """

        # Properties
        self.cfg = cfg
        self.use_daily_scaling = use_daily_scaling
        self.data = None
        self.filepaths = []
        self.error = ErrorStatus()

    def load(self):
        """
        Load the required data. This will load the data for all winter month into memory and the return
        either a weighted fiels (if `use_daily_scaling` is True) or just the field from the corresponding month
        :return:
        """

        # Check if data is already loaded
        if self.has_data_loaded:
            return

        # Load the data of all month
        self.data = []
        for month_num in self.month_nums:

            # Get the target file path
            filepath = self.get_filepath(month_num)

            # Read the data set (and raise hard error if input is missing)
            try:
                nc = open_dataset(filepath)
                self.data.append(nc)
                self.filepaths.append(filepath)
            except FileNotFoundError:
                msg = "Could not locate file: {}".format(filepath)
                self.error.add_error("invalid-filepath", msg)
                self.error.raise_on_error()

    def get_lonlat(self):
        """
        Return longitude and latitude variables
        :return:
        """

        # The grid is the same for all month, therefore we can just retrieve the fields
        # from the first data sets
        dset = self.data[0]
        return dset.longitude.values, dset.latitude.values

    def get_var(self, parameter_name, date_tuple):
        """
        Get the a geophysical variable from the netCDF(s). If daily scaling is activated, the date information
        given by date tuple will be used to create output fields that are interpolated between adjacent month.
        :param parameter_name:
        :param date_tuple:
        :return:
        """

        # There are three cases that requires a different handling:
        #
        # 1. daily scaling is off
        #    -> return the single field of the single data set for the corresponding month
        if not self.use_daily_scaling:
            return self.get_monthly_field(date_tuple[1], parameter_name)

        # 2. daily scaling is on and requested date is a reference date
        #    -> return the field of the single data set for the reference date
        is_reference_date = date_tuple[1:] in self.reference_dates
        if self.use_daily_scaling and is_reference_date:
            return self.get_monthly_field(date_tuple[1], parameter_name)

        # 3. daily scaling is on and requested date is between reference dates
        #   -> return a linear interpolated field based on the distance to the two enclosing
        #      reference dates
        if self.use_daily_scaling and not is_reference_date:
            return self.get_weighted_variable(date_tuple, parameter_name)

    def get_filepath(self, month_num):
        """
        Return the file path for a given month
        :param month_num: Number of month (1-12)
        :return:
        """

        # Create a dictionary for automatic filepath completion
        date_dict = dict(month="{:02g}".format(month_num))

        # Main directory
        path = Path(self.cfg.local_repository)

        # Add the subfolders
        for subfolder_tag in self.cfg.subfolders:
            subfolder = date_dict[subfolder_tag]
            path = path / subfolder

        # Get the period dict (will be constructed from filenaming)
        period_dict = {}
        attrs = re.findall("{.*?}", self.cfg.filenaming)
        for attr_def in attrs:
            attr_name = attr_def[1:-1]
            period_dict[attr_name] = date_dict[attr_name]
        filename = self.cfg.filenaming.format(**period_dict)
        path = path / filename
        return path

    def get_monthly_field(self, month_num, parameter_name):
        """
        Return the monthly field for given parameter name
        :param month_num:
        :param parameter_name:
        :return:
        """
        index = self.month_nums.index(month_num)
        variable = getattr(self.data[index], parameter_name, None)
        if variable is None:
            msg = "Dataset has no variable: {}".format(parameter_name)
            self.error.add_error("invalid-variable", msg)
            self.error.raise_on_error()
        return variable.values

    def get_reference_month_nums(self, date_tuple):
        """
        Return the two month required for the interpolation.
        :param date_tuple: [year, month, day] as integer
        :return: month_left, month_right, weight_factor
        """

        # Compute the difference in days between requested days
        requested_date_dt = datetime(*date_tuple)

        # Sanity check for period coverage
        if requested_date_dt.month not in self.month_nums:
            logger.error("Target month is outside data coverage, snow depth will be NaN")
            return 10, 11, np.nan

        ref_datetimes = self.get_reference_datetimes(date_tuple)
        ref_date_offset = [(requested_date_dt-dt).days for dt in ref_datetimes]

        # Find the index of the first months, where the difference in day is 0 or less (right boundary)
        month_right_index = int(np.argmax(np.array(ref_date_offset) <= 0))
        month_left_index = month_right_index - 1
        month_left, month_right = self.month_nums[month_left_index], self.month_nums[month_right_index]

        # Compute the weighting factor
        period_n_days = (ref_datetimes[month_right_index] - ref_datetimes[month_left_index]).days
        weight_factor = float(ref_date_offset[month_left_index])/float(period_n_days)

        # All done
        return month_left, month_right, weight_factor

    def get_reference_datetimes(self, date_tuple):
        """
        Creates datetimes objects for the reference dates for the actual winter season
        :param date_tuple:
        :return:
        """

        # Get the winter id (year of October for October - April winter)
        winter_id = date_tuple[0] - int(date_tuple[1] < 10)
        year_vals = [winter_id]*3 + [winter_id+1]*4
        ref_dts = [datetime(yyyy, mm, dd) for yyyy, (mm, dd) in zip(year_vals, self.reference_dates)]
        return ref_dts

    def get_weighted_variable(self, date_tuple, parameter_name):
        """
        Compute the weighted variable between two reference dates
        :param date_tuple:
        :param parameter_name:
        :return:
        """

        # Get the fields of both reference month
        month_num_left, month_num_right, weight_factor = self.get_reference_month_nums(date_tuple)
        var_left = self.get_monthly_field(month_num_left, parameter_name)
        var_right = self.get_monthly_field(month_num_right, parameter_name)

        # Get the relative distance (0: var_left, 1: var_right)
        var = var_left + weight_factor*(var_right-var_left)

        # Done
        return var

    @property
    def w99_weight(self):
        """
        Return the static regional mask for the merged climatology
        :return:
        """
        return self.data[0].w99_weight.values

    @property
    def has_data_loaded(self):
        """
        Status flag if data is present for the current data period
        :return:
        """
        return self.data is not None

    @property
    def month_nums(self):
        return [10, 11, 12, 1, 2, 3, 4]

    @property
    def reference_dates(self):
        """
        Return the reference dates for the
        :return:
        """
        return [[10, 1],   # October 1st (to get full coverage of October)
                [11, 15],
                [12, 15],
                [1, 15],
                [2, 15],
                [3, 15],
                [4, 30]]   # April 30th (to get full coverage of April)