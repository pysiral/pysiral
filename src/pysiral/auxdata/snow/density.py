

from datetime import datetime

import numpy as np

from pysiral.auxdata import AuxdataBaseClass, GridTrackInterpol


class SeasonalArcticSnowDensityMallett2020(AuxdataBaseClass):
    """
    Returns a seasonal varying but regionally static snow density from Mallett et al. 2020. The density
    is computed as a linear function of time (elapsed month since the October for the Arctic winter season):

        rho_snow = 6.50 * t + 274.51

    with t: month since October.

    Reference:

        Mallett, R. D. C., Lawrence, I. R., Stroeve, J. C., Landy, J. C., and Tsamados, M.:
        Brief communication: Conventional assumptions involving the speed of radar waves in snow
        introduce systematic underestimates to sea ice thickness and seasonal growth rate estimates,
        The Cryosphere, 14, 251–260, https://doi.org/10.5194/tc-14-251-2020, 2020.

    Notes:

        - This parametrization is valid stricly only for the northern hemisphere

        - Here we interpret the time t not as month since October but as fractional month
          (computed via days since Oct 15) to avoid artifical jumps of geophysical paramters
          between the last day of a month and the first of a next one.

    Example entry in l2 proc config files:

        - snow:
            name: snow_density_seasonal_mallett
            options:
                snow_density_uncertainty: 50.

    """

    def __init__(self, *args, **kwargs):
        super(SeasonalArcticSnowDensityMallett2020, self).__init__(*args, **kwargs)

    def subclass_init(self):
        """
        Not requires
        :return:
        """
        pass

    def get_l2_track_vars(self, l2):
        """
        Mandatory method for the Level-2 processor. Registers snow density (and only density)
        as the auxiliary parameter.
        :param l2:
        :return:
        """

        # Set the requested date
        self.set_requested_date_from_l2(l2)

        # --- Compute the factor t for the date of the l2 object ---
        date_tuple = self._requested_date
        winter_id = date_tuple[0] - int(date_tuple[1] < 10)   # The year of the winter in October

        # Number of days since October 15.
        epoch, l2date = datetime(winter_id, 10, 15), datetime(*date_tuple)
        n_days = (l2date-epoch).days

        # Compute the daily t-factor based on that t should the value of 6 on April 15.
        # of the winter seasons. There are 182 days between Oct. 15 and the following
        # 15th of April.
        #
        # Notes:
        #   - In case of a leap year the factor t will be slightly above 6 on April 15.
        #     This is intended.
        #   - The t factor will be negative before Oct. 15. This is intended
        t = float(n_days)/182. * 6.

        # Compute the fixed snow depth
        static_sdens = 6.5 * t + 274.51

        # Get the snow density uncertainty
        if "snow_density_uncertainty" in self.cfg.options:
            static_sdens_unc = self.cfg.options.snow_density_uncertainty
        else:
            msg = ": Warning - snow_density_uncertainty missing in processor definition (set to 0.0)"
            self.add_handler_message(self.__class__.__name__ + msg)
            static_sdens_unc = 0.0

        # Create variables
        sdens = np.full(l2.n_records, static_sdens)
        sdens_unc = np.full(l2.n_records, static_sdens_unc)

        # Only use values for records, that are sea ice
        no_sea_ice = np.isnan(l2.sic[:])
        for var in [sdens, sdens_unc]:
            var[no_sea_ice] = np.nan

        # Register the auxiliary variable
        self.register_auxvar("sdens", "snow_density", sdens, sdens_unc)


class FixedSnowDepthDensity(AuxdataBaseClass):
    """
    Returns constant depth & density (values from l2 processor definition).

    Example entry in l2 proc config files:

        - snow:
            name: constant
            options:
                snow_density_uncertainty: 50

    TODO: Add uncertainties
    """

    def __init__(self, *args, **kwargs):
        super(FixedSnowDepthDensity, self).__init__(*args, **kwargs)

    def subclass_init(self):
        pass

    def get_l2_track_vars(self, l2):
        snow_depth = np.full(l2.n_records, self.cfg.options.fixed_snow_depth)
        snow_density = np.full(l2.n_records, self.cfg.options.fixed_snow_density)
        snow = SnowParameterContainer()
        snow.depth = snow_depth
        snow.density = snow_density
        # Register Variables
        self.register_auxvar("sd", "snow_depth", snow.depth, snow.depth_uncertainty)
        self.register_auxvar("sdens", "snow_density", snow.density, snow.density_uncertainty)