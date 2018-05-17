import numpy as np

from pysiral.catalog import L3CProductCatalog
from pysiral.errorhandler import ErrorStatus
from pysiral.filter import smooth_2darray
from pysiral.logging import DefaultLoggingClass
from pysiral.iotools import ReadNC



class L3ParameterCollection(DefaultLoggingClass):

    def __init__(self, variable_name, repo_dir=None, ctlg=None, squeeze_empty_dims=True):

        super(L3ParameterCollection, self).__init__(self.__class__.__name__)

        # Name of the parameter from the netCDF files
        self.variable_name = variable_name
        self.squeeze_empty_dims = squeeze_empty_dims
        self._product = {}
        self._mask = {}
        self.error = ErrorStatus()

        # Simple consistency check
        if repo_dir is None and ctlg is None:
            msg = "Either repo_dir or ctlg must be specified"
            self.error.add_error("invalid-args", msg)
            self.error.raise_on_error()

        if repo_dir is not None and ctlg is not None:
            msg = "Both repo_dir or ctlg are specified, using ctlg"
            self.error.add_error("invalid-args", msg)
        
        # Construct L3 product catalog
        if repo_dir is not None:
            self.ctlg = L3CProductCatalog(repo_dir)

        # use existing catalog
        if ctlg is not None:
            self.ctlg = ctlg

        # parse the files
        for product in self.ctlg.product_list:
            nc = ReadNC(product.path)
            var = getattr(nc, self.variable_name)
            if self.squeeze_empty_dims:
                var = np.squeeze(var)
            self._product[product.id] = L3Parameter(self.variable_name, var, product)
            self.log.info("Add product: %s"% product.id)

    def get_monthly_mean(self, month_num, exclude_years=None, **kwargs):
        # Get list of product ids for given month
        product_ids = self.ctlg.get_month_products(month_num, exclude_years=exclude_years)
        return self._get_product_ids_mean(product_ids, **kwargs)

    def get_month(self, year_num, month_num, **kwargs):
        
        # Search for all products in give month
        product_ids = self.ctlg.get_time_range_ids([year_num, month_num], [year_num, month_num])
        
        # Create product stack
        n_products = len(product_ids)

        if n_products == 0:
            return None
        else:
            return self._get_product_ids_mean(product_ids, **kwargs)

    def get_monthly_anomaly(self, year_num, month_num, exclude_ref_month=True, filter_width=None):
        
        # Get mean conditions
        # TODO: Add exclude_ref_month
        if exclude_ref_month:
            exclude_year = [year_num]
        else:
            exclude_year = []
        mean_grid = self.get_monthly_mean(month_num, exclude_years=exclude_year)

        # Get target month
        month_grid = self.get_month(year_num, month_num)
        
        if month_grid is not None:
            anomaly = month_grid-mean_grid
        else:
            anomaly = np.full(mean_grid.shape, np.nan)

        if filter_width is not None:
            anomaly = smooth_2darray(anomaly, filter_width=filter_width)
            
        return anomaly

    def get_winter_season_mean(self, start_year, anomaly=False, **kwargs):
        """ Returns the mean of a winter season (Oct. start_year till Apr. start_year+1)
        for the parameter or the anomaly (anomaly=True) """
        
        # TODO: Simple solution for the moment
        winter_season_month = [10, 11, 12, 1, 2, 3, 4]
        n_month = len(winter_season_month)

        # Init the array
        dims = [n_month]
        dim = dims.extend(self.grid_dims)
        grid_array = np.full(dims, np.nan)


        year_num = start_year
        for i, month_num in enumerate(winter_season_month):

            # Increase year number in January
            if month_num == 1:
                year_num += 1

            if anomaly:
                grid_array[i, :, :] = self.get_monthly_anomaly(year_num, month_num, **kwargs)
            else:
                grid_array[i, :, :] = self.get_month(year_num, month_num, **kwargs)

        return np.nanmean(grid_array, axis=0)


    def get_by_period_id(self, period_id, raise_if_multiple=True):
        """ Returns product(s) for a given period id (see catalog.ProductMetadata) """
        prd_list = [prd for prd in self.products if prd.ctlg.period_id == period_id]
        # Post processing
        if len(prd_list) == 0:
            value = None
        elif len(prd_list) == 1:
            value = prd_list[0]
        else:
            value = prd_list
            if raise_if_multiple:
                self.error.add_error("l3collect-multiple-products-per-period-id", "Multiple entries for period_id: %s" % period_id)
                self.error.raise_on_error()
        return value 

    def set_mask(self, mask, mask_name):
        """ Adds a boolean array that needs to be of the same dimension than
        the product (True=Masked) """
        for prd in self.products:
            prd.set_mask(mask, mask_name)

    def _get_product_ids_mean(self, product_ids, filter_width=None):
        """ Return the mean of a list of product ids """
        dims = self.grid_dims
        shape = (len(product_ids), dims[0], dims[1])
        varstack = np.full(shape, np.nan)
        for i, product_id in enumerate(product_ids):
            varstack[i, :, :] = self._product[product_id].variable
        meanvar = np.nanmean(varstack, axis=0)
        if filter_width is not None:
            meanvar = smooth_2darray(meanvar, filter_width)
        return meanvar

    @property
    def grid_dims(self):
        product_id = self.ctlg.product_ids
        return self._product[product_id[0]].grid_dims

    @property
    def years(self):
        years = sorted([prd.ctlg.tcs.year for prd in self.products])
        return np.unique(years)

    @property
    def products(self):
        return [self._product[prd] for prd in sorted(self._product.keys())]

    @property
    def platforms(self):
        platforms = sorted([prd.ctlg.platform for prd in self.products])
        return np.unique(platforms)


class L3Parameter(DefaultLoggingClass):

    def __init__(self, variable_name, variable, ctlg):
        self.variable_name = variable_name
        self._variable = variable
        self._masks = {}
        self.ctlg = ctlg

    def set_mask(self, mask, mask_name):
        self._masks[mask_name] = mask

    def _get_masked_variable(self):
        masked_variable = np.copy(self._variable)
        for mask_name in self.mask_names:
            mask = self._masks[mask_name]
            masked_variable[np.where(mask)] = np.nan          
        return masked_variable

    @property
    def variable(self):
        return self._get_masked_variable()

    @property
    def grid_dims(self):
        return self._variable.shape

    @property
    def has_mask(self):
        return len(self.mask_names) > 0

    @property
    def mask_names(self):
        return sorted(self._masks.keys())


class L3ParameterPair(DefaultLoggingClass):

    def __init__(self, param_a, param_b):
        self.param_a = param_a
        self.param_b = param_b

    def get_common_grid_points(self):
        """ Return a vector of all grid values with coverage in both pairs """
        a_has_data = np.isfinite(self.param_a.variable)
        b_has_data = np.isfinite(self.param_b.variable)
        common_grid_indices = np.where(np.logical_and(a_has_data, b_has_data))
        return np.ravel(self.param_a.variable[common_grid_indices]), np.ravel(self.param_b.variable[common_grid_indices])

    @property
    def grid_dims(self):
        return self.variable.shape


class L3ParamPairCollection(DefaultLoggingClass):

    def __init__(self, collect_a, collect_b):

        super(L3ParamPairCollection, self).__init__(self.__class__.__name__)
        
        # Init parameters
        self._l3_pairs = {}

        # Store parameters
        self._collect_a = collect_a
        self._collect_b = collect_b

        # Search for overlap between both collections
        self._get_pairs()

    def get_all_pairs(self):
        return self._get_point_list(self.pairs)

    def get_month_pairs(self, month_num):
        pairs = [pair for pair in self.pairs if pair.param_a.ctlg.tcs.month == month_num]
        return self._get_point_list(pairs)

    def _get_point_list(self, pair_list):
        a = np.array([])
        b = np.array([])
        for pair in pair_list:
            prd_a, prd_b = pair.get_common_grid_points()
            a, b = np.append(a, prd_a), np.append(b, prd_b)
        return a, b

    def _get_pairs(self):
        # Loop over all product in collection A and see if counterpart exists in collection B
        for prd_a in self._collect_a.products:
            prd_b = self._collect_b.get_by_period_id(prd_a.ctlg.period_id)
            if prd_b is not None:
                self._l3_pairs[prd_a.ctlg.period_id] = L3ParameterPair(prd_a, prd_b)

    @property
    def pairs(self):
        return [self._l3_pairs[period_id] for period_id in self.period_ids]

    @property
    def period_ids(self):
        return sorted(self._l3_pairs.keys())