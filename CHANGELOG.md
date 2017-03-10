# History of changes

## Version 0.3.0 (10. March 2017)

**New Features**
* [l2 processing] Added SICCI-2 retracker class for Envisat (retracker.SICC2TfmraEnvisat)
* [l2 processing] Added SICCI-2 surface type class for Envisat and CryoSat-2 (surface_type.SICCI2Envisat, surface_type.SICCI2Cryosat2)
* [l2 processing] Added support for NASA-Team based MYI concentrations produced by Integrated Climate Data Center (ICDC) (-> sitype.ICDCNasaTeam)
* [l2 processing] Added support for Antarctic daily snow climatology by Integrated Climate Data Center (ICDC) (-> snow.ICDCSouthernClimatology)
* [l2 processing] Added error catching and reporting for algorithm classes. Run statistics, software environment parameters and a breakdown of encountered error is logged and written to `pysiral-l2proc-summary.txt` in the l2i export directory
* [l2 processing] Added capability to pass a list of input files using the `-l1b-files` argument (mutually exclusive with `-start` and `-stop`)
* [l2 processing] surface type classification algorithm can now access radar mode flag
* [l2 processing] retracker algorithm can now retrieve parameters from l1b object
* [l3 processing] Added surface type statistics (level-3) parameter (`n_total_waveforms`, `n_valid_waveforms`, `valid_fraction`, `lead_fraction`, `ice_fraction`)
* [l3 processing] Added sea ice concentration masking (grid parameter to nan if `sic < 5` or `sic == nan`)
* [l3 processing] Added netCDF variable attributes to L3s output
* [tools] Added customization options for the colorbar generator
* [scripts] Added bash script to download, extract and sort CryoSat-2 NRT files (requires username & password for ESA ftp)

**Bugfixes**
* [config] Bug #280: month exclusion ignored for the month given by -start option in the l1bpreproc, l2proc and l3proc calls leading to incorrect start and stop dates of the processors

**Changes**
* [l2 processing] reorganization of l2 processing scheme to allow computation of uncertainties
* [l2 processing] Improved error handling and reporting
* [l3 processing] Improved logging and documentation
* [l3 processing] Additional field `sea_ice_concentration_mask_targets` in level-3 output definition required
* [visualization] Minor changes to colormaps
* [visualization] Added sea ice concentration to map background
* [visualization] Added capability to create Antarctic maps

**Settings**
* [level-2] Added SICCI-2 prototype 002 settings for Arctic and Antarctic
* [level-2] Added SICCI-2 prototype 003 settings for Arctic and Antarctic


## Version 0.2.0 (19. September 2016)

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


