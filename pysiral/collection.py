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
            self._product[product.id] = L3Parameter(self.variable_name, var)
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

    def get_monthly_anomaly(self, year_num, month_num, exclude_ref_month=False, filter_width=None):
        
        # Get mean conditions
        # TODO: Add exclude_ref_month
        if exclude_ref_month:
            exclude_year = [year_num]
        else:
            exclude_year = False
        mean_grid = self.get_monthly_mean(month_num, exclude_years=exclude_year)

        # Get target month
        month_grid = self.get_month(year_num, month_num)
        
        anomaly = month_grid-mean_grid

        if filter_width is not None:
            anomaly = smooth_2darray(anomaly, filter_width=filter_width)
            
        return anomaly

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


class L3Parameter(DefaultLoggingClass):

    def __init__(self, variable_name, variable):
        self.variable_name = variable_name
        self.variable = variable

    @property
    def grid_dims(self):
        return self.variable.shape