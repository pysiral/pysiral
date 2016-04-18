# -*- coding: utf-8 -*-


from pysiral.path import file_basename

from netCDF4 import Dataset, date2num
from datetime import datetime
import numpy as np
import os


class NCDateNumDef(object):
    """
    Holds definition for datetime conversion to numbers and vice versa
    for netCDF operations
    """

    def __init__(self):
        self.units = "seconds since 1970-01-01"
        self.calendar = "standard"


class L1bDataNC(object):
    """
    Class to export a L1bdata object into a netcdf file
    """

    def __init__(self):
        self.datagroups = ["waveform", "surface_type", "time_orbit",
                           "classifier", "correction"]
        self.output_folder = None
        self.l1b = None
        self.filename = None
        self.zlib = False
        self.time_def = NCDateNumDef()
        self._rootgrp = None

    def export(self):
        self._validate()
        self._create_root_group()
        self._populate_data_groups()
        self._write_to_file()

    def _validate(self):
        if self.filename is None:
            self._create_filename()
        self.path = os.path.join(self.output_folder, self.filename)

    def _create_filename(self):
        self.filename = file_basename(self.l1b.filename)+"nc"

    def _create_root_group(self):
        """
        Create the root group and add l1b metadata as global attributes
        """
        self._rootgrp = Dataset(self.path, "w")
        # Save the l1b info data group as global attributes
        attdict = self.l1b.info.attdict
        self._convert_datetime_attributes(attdict)
        self._convert_bool_attributes(attdict)
        self._convert_nonetype_attributes(attdict)
        self._set_global_attributes(attdict)

    def _populate_data_groups(self):
        for datagroup in self.datagroups:
            print datagroup.upper()
            # Create the datagroup
            dgroup = self._rootgrp.createGroup(datagroup)
            content = getattr(self.l1b, datagroup)
            # Create the dimensions
            # (must be available as OrderedDict in Datagroup Container
            dims = content.dimdict.keys()
            for key in dims:
                dgroup.createDimension(key, content.dimdict[key])
            # Now add variables for each parameter in datagroup
            for parameter in content.parameter_list:
                data = getattr(content, parameter)
                # Convert datetime objects to number
                if type(data[0]) is datetime:
                    data = date2num(data, self.time_def.units,
                                    self.time_def.calendar)
                # Convert bool objects to integer
                if data.dtype.str == "|b1":
                    data = np.int8(data)
                dimensions = tuple(dims[0:len(data.shape)])
                print " "+parameter, dimensions, data.dtype.str
                var = dgroup.createVariable(
                    parameter, data.dtype.str, dimensions, zlib=self.zlib)
                var[:] = data

    def _write_to_file(self):
        self._rootgrp.close()

    def _convert_datetime_attributes(self, attdict):
        """
        Replace l1b info parameters of type datetime.datetime by a double
        representation to match requirements for netCDF attribute data type
        rules
        """
        for key in attdict.keys():
            content = attdict[key]
            if type(content) is datetime:
                attdict[key] = date2num(
                    content, self.time_def.units, self.time_def.calendar)

    def _convert_bool_attributes(self, attdict):
        """
        Replace l1b info parameters of type bool ['b1'] by a integer
        representation to match requirements for netCDF attribute data type
        rules
        """
        for key in attdict.keys():
            content = attdict[key]
            if type(content) is bool:
                attdict[key] = int(content)

    def _convert_nonetype_attributes(self, attdict):
        """
        Replace l1b info parameters of type bool ['b1'] by a integer
        representation to match requirements for netCDF attribute data type
        rules
        """
        for key in attdict.keys():
            content = attdict[key]
            if content is None:
                attdict[key] = ""

    def _set_global_attributes(self, attdict):
        """ Save l1b.info dictionary as global attributes """
        for key in attdict.keys():
            print key, type(attdict[key])
            self._rootgrp.setncattr(key, attdict[key])


def l1bnc_filenaming(l1b, config):
    """
    Returns the standard export folder and filename of a level-1b netCDF
    data file

    Export Folder
    -------------

    Based on the ``l1bdata`` file for the specific mission in
    ``local_machine_def.yaml`` in the pysiral main directory
    with sub-folders for hemisphere, year, month

    Filename
    --------

    Default l1bdata filename in pysiral::

        l1bdat_v$VERS_$REG_$MISSION_$ORBIT_$YYYYMMDDHHMISS_$YYYYMMDDHHMISS.nc

    :$VERS: l1bdata version [00 (beta)]
    :$REG: region [north | south]
    :$MISSION: mission short name [cryosat2 | envisat | ers1 | ...]
    :$ORBIT: orbit/cylce number [ 00026000 ]
    :$YYYYMMDDHHMISS: start and end time

    """
    # export folder: $mission_l1bdata_folder/YYYY/MM (start time)
    export_folder = config.local_machine.l1b_repository[l1b.mission].l1bdata
    yyyy = "%04g" % l1b.info.start_time.year
    mm = "%02g" % l1b.info.start_time.month
    hemisphere = l1b.info.hemisphere
    export_folder = os.path.join(export_folder, hemisphere, yyyy, mm)
    # construct filename from
    export_filename = "l1bdata_v{version:02g}_{region}_{mission}_" + \
                      "{orbit:06g}_{startdt:%Y%m%dT%H%M%S}_" + \
                      "{startdt:%Y%m%dT%H%M%S}.nc"
    export_filename = export_filename.format(
        version=0, region=hemisphere, mission=l1b.info.mission,
        orbit=l1b.info.orbit, startdt=l1b.info.start_time,
        stopdt=l1b.info.stop_time)
    return export_folder, export_filename
