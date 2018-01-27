""" Catalog module for Level-2 and Level-2 product repositories"""

import os
import sys
import fnmatch
import uuid
import time

import numpy as np

from pysiral.iotools import ReadNC
from pysiral.logging import DefaultLoggingClass

class L2PProductCatalog(DefaultLoggingClass):
    """Catalog class for L2P product repositories

    Arguments:
            repo_path {str} -- path to repository"""
    
    def __init__(self, repo_path):
        super(L2PProductCatalog, self).__init__(self.__class__.__name__)
        self.repo_path = repo_path
        self._catalog = {}
        self._catalogize()

    def _catalogize(self):
        """Create the product catalog of the repository"""

        # Get the list of netcdf product files
        nc_files = self.nc_files

        # Create the catalog entries as dictionary with product id as keys
        self.log.info("Catalogizing l2p repository: %s (%g files)" % (str(self.repo_path), len(nc_files)))
        t0 = time.clock()
        for nc_file in self.nc_files:
            product_info = ProductMetadata(nc_file, target_processing_level="l2p")
            self._catalog[product_info.id] = product_info
        t1 = time.clock()
        self.log.info("... done in %.1f seconds" % (t1-t0))

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
    def n_product_files(self):
        return len(self._catalog.keys())

    @property
    def product_ids(self):
        return sorted(self._catalog.keys())

    @property
    def product_list(self):
        return [self._catalog[idstr] for idstr in self.product_ids]

    @property
    def versions(self):
        return [prd.product_version for prd in self.product_list]

    @property
    def has_unique_version(self):
        return len(np.unique(self.versions)) == 1


class ProductMetadata(DefaultLoggingClass):
    """Metadata data container for pysiral product files."""

    VALID_PROCESSING_LEVELS = ["l2i", "l2p", "l3c"]
    VALID_CDM_DATA_LEVEL = ["Trajectory", "Grid"]
    NC_PRODUCT_ATTRIBUTES = [
        "processing_level", "product_version", "cdm_data_level", 
        "time_coverage_start", "time_coverage_end", "product_timeliness",
        "time_coverage_duration", "source_mission_id", "source_hemishere"]

    def __init__(self, local_path, target_processing_level=None):
        """
        Arguments:
            local_path {str} -- local path to the product netcdf
        
        Keyword Arguments:
            target_processing_level {str} -- Target processing level for the product netcdf. Settings 
                                             a specific processing level will cause an exception in 
                                             the case of a mismatch (default: {None})
        """
        super(ProductMetadata, self).__init__(self.__class__.__name__)

        self.local_path = local_path

        if target_processing_level in self.VALID_PROCESSING_LEVELS or target_processing_level is None:
            self._targ_proc_lvl = target_processing_level
        else:
            raise ValueError("Invalid target processing level: %s" % str(target_processing_level))

        # Fetch attributes (if possible) from netcdf
        nc = ReadNC(self.local_path)
        for attribute in self.NC_PRODUCT_ATTRIBUTES:
            
            # Extract value from netcdf global attributes
            try:
                value = getattr(nc, attribute)
            except AttributeError:
                value = None

            # Now follow a few special rules from some attributes
            if attribute == "processing_level":
                value = self._validate_proc_lvl(value)

            # Save the attributes
            setattr(self, attribute, value)

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
        # (search for the occurance of all valid processing levels in the processing_level attribute)
        pl_match = [pl in str(value).lower() for pl in self.VALID_PROCESSING_LEVELS]
        try: 
            index = pl_match.index(True)
        except ValueError:
            self.log.error("Invalid processing level str (%s) in product: %s" % (str(value), self.local_path))
            sys.exit()
        value = self.VALID_PROCESSING_LEVELS[index]

        if value != self._targ_proc_lvl and self._targ_proc_lvl is not None:
            msg = "Invalid product processing level: %s (target: %s)"
            raise ValueError(msg % (value, self._targ_proc_lvl))

        return value

    @property
    def id(self):
        """Generates an id string for the product.
        
        Returns:
            [str] -- id str of the product
        """
        random_str = str(uuid.uuid4())
        identifier = (
            self.processing_level, str(self.source_mission_id),
            self.time_coverage_start, self.time_coverage_end,
            random_str[0:8])
        idstr = "%s-%s-%s-%s-%s" % identifier
        return idstr