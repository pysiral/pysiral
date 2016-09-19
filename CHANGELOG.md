# History of changes

## Version 0.2.0 (TBD)

**New Features**
* [l1b pre-processing] time range for can be specified for individual days
* [l1b pre-processing] pre-processor can be set to specific hemispheres
* [l1b pre-processing] obsolete l1bdata files can be automatically removed (`--remove-old` option to pysiral-l1bpreproc)
* [l1b pre-processing] new classifier parameter for CryoSat-2, Envisat and ERS
* [l2 processing] time range for can be specified for individual days
* [l2 processing] obsolete l2 output files can be automatically removed (`--no-overwrite-protection` & `--remove-old` options to pysiral-l2proc)
* [l2 processing] l2i output files log the settings of the level2 processor in the global attributes
* [l2 processing] l2i output variables with proper netCDF attributes
* [retracker] cythonized implementation of the TFRMA (3 times faster than pure python)
* [tools] `l1bdata_report.py`: creates a pdf summary of l1bdata files
* [tools] `zipl2month.py` creates monthly archives of l2i orbit files

**Bugfixes**

* [l1b pre-processing] TAI to UTC timestamp conversion for CryoSat-2 incorrect
* [l1b pre-processing] Wrong values for open ocean percentage when subsetting and merging l1b segments
* [l1b pre-processing] pysiral version tag was not properly passed to l1bdata netcdf files

**Changes**

* [l1b pre-processing] moved hard coded settings from python parser to mission config file (Envisat, ERS, Sentinel-3)
* [l1b pre-processing] All geophysical corrections are now piped into l1bdata files instead of selection in the mission configuration. Requires changes in the l2 settings files.
* [l2 processing] Internal handling of processor settings


