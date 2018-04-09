# -*- coding: utf-8 -*-
""" Catalog module for Level-2 and Level-2 product repositories"""

import os
import re
import sys
import fnmatch
import uuid
import time
import datetime
import dateutil
import numpy as np

from pysiral.iotools import ReadNC
from pysiral.logging import DefaultLoggingClass
from pysiral.config import TimeRangeRequest
from pysiral.errorhandler import ErrorStatus


class SIRALProductCatalog(DefaultLoggingClass):
    """ Parent catalog class for product catalogs for different
    data processing levels """

    def __init__(self, repo_path, auto_id=True, repo_id=None):
        super(SIRALProductCatalog, self).__init__(self.__class__.__name__)
        self.error = ErrorStatus(caller_id=self.__class__.__name__)
        self.repo_path = repo_path
        self.auto_id = auto_id
        self._repo_id = repo_id
        self._catalog = {}

    def run_checks(self, check_list, raise_on_failure=True):
        """Runs a list of built-in queries. Valid checks: `is_single_hemisphere`, `is_single_version`
        
        Arguments:
            check_list {str list} -- list of checks to run
        
        Keyword Arguments:
            raise_on_failure {bool} -- Raise an SystemExi exception if a check is negative (default: {True})
        
        Returns:
            [bool list] -- flag list with check results (True: check passed) Only returned when `raise` keyword is set to `False`)
        """

        check_passed = np.zeros(np.shape(check_list), dtype=bool)
        for index, check in enumerate(check_list):
            check_passed[index] = getattr(self, check)
            try: 
                check_passed[index] = getattr(self, check)
            except AttributeError:
                self.log.error("invalid check: %s" % str(check))
                check_passed[index] = False
            finally:
                if raise_on_failure and not check_passed[index]:
                    self.log.error("failed check: %s" % str(check))
                    sys.exit()

        return check_passed

    def append(self, ctlg, duplication=False):
        """Add the information from another catalog with rules for handling of duplicates 
        (same period, same mission)
        
        Arguments:
            ctlg {object} -- of type SIRALProductCatalog
        
        Keyword Arguments:
            duplication {bool} -- Allows (True) / Prevents (False) duplicate entries (same platform, same period) in the merged catalog (default: {False})
        
        Returns:
            
        """

        # Perform Consistency checks
        # 1. argument need to be of type SIRALProductCatalog
        if not isinstance(ctlg, (L2PProductCatalog, L3CProductCatalog)):
            msg = "Invalid catalog (%s), must be from pysiral.catalog module"
            msg = msg % (str(ctlg.__class__.__name__))
            self.error.add_error("ctlg-invld-ctlg", msg)
            self.error.raise_on_error()

        # 2. Catalogs need to be of the same processing level (l2p, l3c, ....)
        if self.processing_level != ctlg.processing_level:
            msg = "Invalid processing level (%s) of new catalog, %s required for appending"
            msg = msg % (str(ctlg.processing_level), str(self.processing_level))
            self.error.add_error("ctlg-proclevel-mismatch", msg)
            self.error.raise_on_error()

        # Merge the catalogs
        for new_product in ctlg.product_list:
            
            # Check if new product is duplication in current catalog
            is_duplication = self._is_duplication(new_product)

            # if duplication ok & is duplication -> add
            if duplication and is_duplication:
                self._catalog[new_product.id] = new_product

            # if not duplication -> add
            elif not is_duplication:
                self._catalog[new_product.id] = new_product

            # don't add (is duplication and duplication not ok)
            else:
                continue

    def query_datetime(self, dt, return_value="bool"):
        """ Searches the repository for products for a given datetime
        
        Arguments:
            dt {datetime} -- datetime definition for the query
        
        Keyword Arguments:
            return_value {str} -- Defines the type of output: `bool` for True/False flag and `products` for
                                  product path (list) (default: {"bool"})
        
        Returns:
            [bool or str] -- see keyword `return`
        """

        if not isinstance(dt, datetime.datetime):
            raise ValueError("Argument dt needs to be datetime (was: %s)" % (type(dt)))

        product_files = [prd.path for prd in self.product_list if prd.has_coverage(dt)]

        if return_value == "products":
            result = product_files
        else:
            result = len(product_files) > 0
        return result

    def query_overlap(self, tcs, tce):
        """ Searches the repository for products that have overlapping coverage
        with a given time range
        
        Arguments:
            tcs {datetime} -- time coverage start of search period
            tce {datetime} -- time coverage end of search period
        
        Keyword Arguments:
            return_value {str} -- Defines the type of output: `bool` for True/False flag and `products` for
                                  product path (list) (default: {"bool"})
        
        Returns:
            [bool or str] -- see keyword `return`
        """

        if not isinstance(tcs, datetime.datetime):
            raise ValueError("Argument tcs needs to be datetime (was: %s)" % (type(tce)))

        if not isinstance(tce, datetime.datetime):
            raise ValueError("Argument tce needs to be datetime (was: %s)" % (type(tce)))

        # Search files
        product_files = [prd.path for prd in self.product_list if prd.has_overlap(tcs, tce)]
        return product_files

    def get_northern_winter_netcdfs(self, start_year):
        """Returns a list for northern winter data for the period october through april
        
        Arguments:
            start_year {int} -- start year for the winter (yyyy-oct till yyyy+1-apr)

        Returns: 
            product_files {str list} -- list of files for specified winter season
        """

        # Construct time range 
        winter_start_tuple = [start_year, 10]
        winter_end_tuple = [start_year+1, 4]
        time_range = TimeRangeRequest(winter_start_tuple, winter_end_tuple, period="custom")

        # Query time range
        product_files = self.query_overlap(time_range.start_dt, time_range.stop_dt)

        # Reporting
        msg = "Found %g %s files for winter season %g/%g"
        msg = msg % (len(product_files), self.processing_level, start_year, start_year+1)
        self.log.info(msg)

        return product_files

    def get_month_netcdfs(self, month):
        """Returns a list all products that have coverage for a given month
        
        Arguments:
            month {int} -- month number (1-Jan, ..., 12:Dec)

        Returns: 
            product_files {tuple list} -- ((year, month), [list of files month])
        """

        # Query time range
        years = self.years
        product_files = []

        n_files = 0
        for year in self.years:
            time_range = TimeRangeRequest((year, month), (year, month))
            tcs, tce = time_range.start_dt, time_range.stop_dt
            files = [prd.path for prd in self.product_list if prd.has_overlap(tcs, tce)]
            n_files += len(files)
            product_files.append(((year, month), files))

        # Reporting
        msg = "Found %g %s files for %s"
        month_name = datetime.datetime(2000, month, 1).strftime("%B")
        msg = msg % (n_files, self.processing_level, month_name)
        self.log.info(msg)

        return product_files

    def _catalogize(self):
        """Create the product catalog of the repository"""

        # Get the list of netcdf product files
        nc_files = self.nc_files

        # Create the catalog entries as dictionary with product id as keys
        self.log.info("Catalogizing %s repository: %s (%g files)" % (self.processing_level, str(self.repo_path), len(nc_files)))

        if self.auto_id:
            subfolders = self.repo_path.split(os.sep)
            try:
                search  = [bool(re.search(self.processing_level, subfolder)) for subfolder in subfolders]
                index = search.index(True)
                repo_id = subfolders[index-1]
            except:
                self.log.warning("Auto id failed, did not find `%s` in list of subfolders" % (self.processing_level))
                repo_id = None
            self._repo_id = repo_id
            self.log.info("%s repository ID: %s" % (self.processing_level, str(self.repo_id)))
        
        nc_access_times = []
        t0 = time.clock()
        for nc_file in self.nc_files:
            product_info = ProductMetadata(nc_file, target_processing_level=self.processing_level)
            nc_access_times.append(product_info.nc_access_time)
            self._catalog[product_info.id] = product_info
        t1 = time.clock()
        self.log.info("... done in %.1f seconds" % (t1-t0))
        self.log.info("... average netCDF access time: %.4f sec" % np.mean(nc_access_times))

    def _is_duplication(self, product_info):
        """ Tests if specified product is a duplicate to the current catalog """
        
        # Condition 1: period exists in current catalog
        period_exists = product_info.period_id

    @property
    def nc_files(self):
        """Lists all netcdf files (*.nc) in the repository.
        
        Returns:
            [str] -- list of netcdf files
        """
        nc_files = []
        for root, dirnames, filenames in os.walk(self.repo_path):
            for filename in fnmatch.filter(filenames, "*.nc"):
                nc_files.append(os.path.join(root, filename))
        return sorted(nc_files)

    @property
    def repo_id(self):
        return str(self._repo_id)

    @property
    def n_product_files(self):
        return len(self._catalog.keys())

    @property
    def product_ids(self):
        return sorted(self._catalog.keys())

    @property
    def period_ids(self):
        return [self._catalog[idstr].period_id for idstr in self.product_ids]

    @property
    def product_list(self):
        return [self._catalog[idstr] for idstr in self.product_ids]

    @property
    def versions(self):
        return [prd.product_version for prd in self.product_list]

    @property
    def hemispheres(self):
        return [prd.source_hemisphere for prd in self.product_list]

    @property
    def hemisphere_list(self):
        return np.unique(self.hemispheres)

    @property
    def is_single_version(self):
        return len(np.unique(self.versions)) == 1

    @property
    def is_single_hemisphere(self):
        return len(self.hemisphere_list) == 1

    @property
    def is_north(self):
        hemisphere_list = self.hemisphere_list
        return self.is_single_hemisphere and hemisphere_list[0] == "north" 

    @property
    def is_south(self):
        hemisphere_list = self.hemisphere_list
        return self.is_single_hemisphere and hemisphere_list[0] == "south" 

    @property
    def years(self):
        years = sorted([prd.time_coverage_start.year for prd in self.product_list])
        return np.unique(years)

    @property
    def time_coverage_start(self):
        tcs = [prd.time_coverage_start for prd in self.product_list]
        return np.min(tcs)

    @property
    def tcs(self):
        """ Abbrevivation for self.coverage_start """
        return self.time_coverage_start

    @property
    def time_coverage_end(self):
        tce = [prd.time_coverage_end for prd in self.product_list]
        return np.max(tce)

    @property
    def tce(self):
        """ Abbrevivation for self.coverage_end """
        return self.time_coverage_end


class L2PProductCatalog(SIRALProductCatalog):
    """Catalog class for L2P product repositories

    Arguments:
            repo_path {str} -- path to repository"""
    
    processing_level = "l2p"
    
    def __init__(self, *args, **kwargs):
        super(L2PProductCatalog, self).__init__(*args, **kwargs)
        self._catalogize()


class L3CProductCatalog(SIRALProductCatalog):
    """Catalog class for L3C product repositories

    Arguments:
            repo_path {str} -- path to repository"""
    
    processing_level = "l3c"
    period_id_level = "daily"
    
    def __init__(self, *args, **kwargs):
        super(L3CProductCatalog, self).__init__(*args, **kwargs)
        self._catalogize()


class ProductMetadata(DefaultLoggingClass):
    """Metadata data container for pysiral product files."""

    VALID_PROCESSING_LEVELS = ["l2i", "l2p", "l3c"]
    VALID_CDM_DATA_LEVEL = ["Trajectory", "Grid"]
    NC_PRODUCT_ATTRIBUTES = [
        "processing_level", "product_version", "cdm_data_level", 
        "time_coverage_start", "time_coverage_end", "product_timeliness",
        "time_coverage_duration", "source_mission_id", "source_hemisphere"]

    def __init__(self, path, target_processing_level=None):
        """
        Arguments:
            local_path {str} -- local path to the product netcdf
        
        Keyword Arguments:
            target_processing_level {str} -- Target processing level for the product netcdf. Settings 
                                             a specific processing level will cause an exception in 
                                             the case of a mismatch (default: {None})
        """
        super(ProductMetadata, self).__init__(self.__class__.__name__)

        self.path = path
        self.unique_str = str(uuid.uuid4())[0:8]

        if target_processing_level in self.VALID_PROCESSING_LEVELS or target_processing_level is None:
            self._targ_proc_lvl = target_processing_level
        else:
            raise ValueError("Invalid target processing level: %s" % str(target_processing_level))

        # Fetch attributes (if possible) from netcdf
        t0 = time.clock()
        nc = ReadNC(self.path, global_attrs_only=True)
        t1 = time.clock()
        self.nc_access_time = t1-t0

        for attribute in self.NC_PRODUCT_ATTRIBUTES:
            
            # Extract value from netcdf global attributes
            try:
                value = getattr(nc, attribute)
            except AttributeError:
                value = None

            # Now follow a few special rules from some attributes
            if attribute == "processing_level":
                value = self._validate_proc_lvl(value)

            if attribute in ["time_coverage_start", "time_coverage_end"]:
                value = self._parse_datetime_definition(value)

            setattr(self, attribute, value)

    def has_coverage(self, dt):
        """Test if datetime object is covered by product time coverage
        
        Arguments:
            dt {datetime} -- A datetime object that will be tested for coverage

        Returns:
            [bool] -- A True/False flag
        """
        flag = dt >= self.time_coverage_start and dt <= self.time_coverage_end
        return flag

    def has_overlap(self, tcs, tce):
        """Test if product has overlap with period defined by start & end datetime
        
        Arguments:
            tcs {datetime} -- A datetime object that defines start of time coverage test
            tce {datetime} -- A datetime object that defines end of time coverage test

        Returns:
            [bool] -- A True/False flag
        """

        # Validity check
        if tce <= tcs:
            raise ValueError("Invalid overlap test: tce <= tcs")

        no_coverage = tce <= self.time_coverage_start or tcs >= self.time_coverage_end
        flag = not no_coverage
        return flag

    def _parse_datetime_definition(self, value):
        """Converts the string representation of a date & time into a 
        datetime instance
        
        Arguments:
            value {str} -- [description]
        
        Returns:
            [datetime] -- [description]
        """
        value = dateutil.parser.parse(value)
        return value

    def _validate_proc_lvl(self, value):
        """Validates the processing level str from the netcdf file: a) only save the id and 
        not any potential description
        
        Arguments:
            value {str} -- The attribute value from the product netcdf
        
        Raises:
            ValueError -- if mismatch between processing level in the product and the target level 
        
        Returns:
            [str] -- The validated processing level id string
        """

        # Make sure the value for processing level is only the id
        # search for the occurance of all valid processing levels in the processing_level attribute
        pl_match = [pl in str(value).lower() for pl in self.VALID_PROCESSING_LEVELS]
        try: 
            index = pl_match.index(True)

        # NOTE: if the processing_level attribute does not exist, there is no choice but to trust the repo
        except ValueError:
            return self._targ_proc_lvl

        # Check if processing level in file matches target processing level
        value = self.VALID_PROCESSING_LEVELS[index]
        if value != self._targ_proc_lvl and self._targ_proc_lvl is not None:
            msg = "Invalid product processing level: %s (target: %s)"
            raise ValueError(msg % (value, self._targ_proc_lvl))

        return value

    def _get_datetime_label(self, dt):
        if self.processing_level in ["l2p", "l3c", "l3s"]:
            return dt.strftime("%Y%m%d")
        else:
            return dt.strftime("%Y%m%d%H%M%S")

    @property
    def id(self):
        """Generates an id string for the product.
        
        Returns:
            [str] -- id str of the product
        """
        
        identifier = (
            self.processing_level, str(self.source_mission_id),
            self.period_id, self.unique_str)
        idstr = "%s-%s-%s-%s" % identifier
        return idstr

    @property
    def tcs_label(self):
        return self._get_datetime_label(self.time_coverage_start)

    @property
    def tce_label(self):
        return self._get_datetime_label(self.time_coverage_end)

    @property
    def period_id(self):
        """ Generates a period id """
        identifier = (self.tcs_label, self.tce_label, self.time_coverage_duration)
        period_id = "%sT%s-%s" % identifier
        return period_id