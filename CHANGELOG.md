# History of changes

## Version 0.5.0 (9 February 2018) 

**Bugfixes**
* [auxdata] Longnames of sitype in `auxdata_def.yaml` were not recognized properly (incorrect key)
* [datahandler] Catch exception caused by feeding empty directory to datahandler
* [l1bpreproc] catch a rare exception when the l1b input file is just 1 data point
* [l2data] radar_mode information was lost during the Level-2 processing
* [l2proc] bugfix: Improper handling of pysiral-l2proc call for specific list of l1b input files
* [l2preproc] Catch exception if error during parsing a l2i file
* [l2preproc] Catch exception if l2i file(s) exist for a single day but none contains valid freeboard data
* [l2preproc] Fixed a bug that caused all uncertainties in l2p products to be zeros
* [l3proc] Catch exception if error during parsing a l2i file
* [l3proc] Fixed a bug in the generation of lat/lon grid coordinates that caused the grid to be incorrect by 50% of the grid size at the outer grid limits
* [sentinel3] Fixed several issues in the Sentinel-3A preprocessor
* [sitype] sitype api method did not return value, uncertainty, msg
* [general] Fixed a few of incorrect variable names

**New Features**
* [catalog] Added catalog module for cataloging/querying pysiral product files.(Currently only working for l2p files, that will be extend to all data levels)
* [config] Added default_week (Monday - Sunday) as period in TimeRangeRequest
* [config] Added duration string in ISO format to TimeRangeRequest
* [config] Added monthly, weekly and daily iterators to TimeRangeRequest
* [config] Added isodate labels to TimeRangeRequest
* [datahandler] Added capability to find all l2i files for an arbitrary time_range
* [general] Implement timeliness attribute throughout the processing chain
* [general] pysiral-wide disabling of python warnings
* [general] keep internal record of official mission names (platform)
* [grid] new grid definition property: pyresample area definition
* [icesat] Added icesat module for converting GLAH13 data into a pysiral compliant data format
* [icesat] Added ICESat pre-processor
* [l1dbata] Added method to fill gaps in Level-1 data 
* [l1dbata] Added method to split l1bdata objects at discontinuities 
* [l2proc] Added class for snow freeboard class for converting radar freeboard to freeboard (pysiral.frb.SnowFreeboardAssumption)
* [l2proc] Added surface type classificator for ICESat based on (Khvorostovsky, K. and Rampal, P.: On retrieving sea ice freeboard from ...)
* [l2data] Added method to init l2data object from l2i netcdf file
* [l3proc] Added option to export grid description parameters as variable in netcdf output
* [l3proc] Added time bound property for Level-3 netcdf output
* [l3proc] Added default time and grid parameters to Level-3 netcdf output
* [l3proc] Allow flipping of Level-3 variable in netcdf output (to be consistent with osisaf cdr)
* [mask] new module handling masks
* [mask] Add mechanic to compute fixed masks for each l3 grid definition 
* [mask] Add l3 mask for Arctic area in which Warren99 is considered valid
* [mask] Add land/sea mask with is_land, is_sea & is_mixed flag
* [sentinel3] Add parser for EUMETSAT L2_WAT product
* [settings] New attribute template {uuid} to assign universally unique identifier (uuid) tracking id for l2/l3 products
* [settings] Attribute template {source_hemisphere} has now a select option to select text based on the hemisphere. (Usage: {source_mission_id:select;id1:Text id1;id2:Text id2;...}`)
* [settings] Attribute template {source_mission_id} has now a select option to select text based on the mission id. (Usage: {source_hemisphere:select;Text north;Text south}`)
* [sitype] Added sitype class for osi-saf sitype cdr products

**Changes**
* [l2data] Renamed `timestamp` variable to `time` due to variable naming conventions (backward compability with older files exist)
* [l2data] Added version tag of the algorithm (required in all l2 settings files)
* [l3proc] Expanded the list of possible Level-3 parameters in preparation to status/quality flags
* [l3proc] Expanded the options for default Level-3 processing with external masks and post-processing
* [l3proc] Allowed to load external masks in Level-3 processor
* [l3proc] Allowed computation of quality indicator flags
* [l3proc] Allowed computation of status flags
* [l3proc] Updated Level-3 settings with improved variable initializuation (init value, dtype)
* [settings] Global attributes in l2/l3 products now follow the order prescribed in the output definition files (was alphabetic)
* [settings] Attribute template options are now separated by semikolon (was comma)
* [settings] Time coverage start/end in l2p files are now always defined by requested data period (full day) and not actually available data points
* [settings] Filenaming in settings files is now depends on data period (e.g. different filenaming for weekly or monthly l3c products)
* [general] Set the creation time for products once at the initialization of the python class (was multiple calls of datetime.now() during output generation)
* [griddef] Added grid name in grid definitions
* [general] Countless code documentation, style and logging improvements

**Settings**
* [l2proc] Added Sentinel-3 l2proc setting files
* [l2proc] Added Level-2 algorithm settings recipe for C3S cdr v1
* [l3proc] Added default Level-3 settings for southern hemisphere without Warren99 is_valid check. Needs to be selected manually with the `-l3-settings` option. 
* [output] Added L3C output definition for C3S
* [output] Added l2p/l3c definition for AWI product in version 2.0
* [output] ccicdr output definitions (l2p/l3c) now follow CCI data standard requirements

## Version 0.4.7 (17 August 2017) [bugfix release]

**Bugfixes**
* [settings] L3C filename was incorrectly written (request period start time instead of iteration period start time)


## Version 0.4.6 (16 August 2017)

**New Features**
* [l2proc] Added input validation for l1bdata files (time coverage and versions)

**Settings**
* [l3proc] Changed SICCI-2 L3C global attributes towards pysiral and Data Discovery standards


## Version 0.4.5 (24 July 2017) [bugfix release]

**Bugfixes**
* [settings] Default uncertainty value for sitype south was missing


## Version 0.4.4 (24 July 2017) [bugfix release]

**Bugfixes**
* [settings] Wrong SICCI-2 L3C version
* [l2proc] ICDC sitype class crashed on missing file


## Version 0.4.3 (17 July 2017) [bugfix release]

**Bugfixes**
* [settings] Invalid default l2i settings solved


## Version 0.4.2 (7 July 2017)

**New Features**
* [sitype] MYI concentration uncertainty support Added
* [l2 processing] Added l1b prefilter to the l2proc workflow
* [l2 processing] OsiSafSIC now allows auto switch between products (e.g. osisaf-409 until 2015-04-15 and osisaf-430 after)

**Bugfixes**
* [l3 processing] #396: sea ice type uncertainty not filtered by freeboard mask in default l3 settings
* [output] #397: valid_fraction missing in SICCI-2 L3C output
* [output] #399: attribute tag `utcnow` not working

**Changes**
* [sit] computation of sea ice density uncertainty updated
* [surface_type] Merged SICCI2Envisat and SICCI2CryoSat2 surface type classes

**Settings**
* [l2 processing] Added Envisat backscatter drift correction as l1b prefilter
* [l2 processing] osisafcdr (merged osisaf-409 & osisaf-430) sic for all ccicdr settings


## Version 0.4.1 (3. July 2017)

**New Features**
* [grid] New `grid` module for grid related tasks (work in progress)
* [l2 data] Added `sea_surface_height` property (ssh = mss + ssa)
* [auxdata] Added a `NoneHandler` (dummy handler that only returns NaN's)
* [l3 processing] Added masking of l2i data before gridding

**Bugfixes**
* [l3 processing] Level-3 processor crash when no l2i data was available for given iteration

**Changes**
* [l2 processing] Updated uncertainty description
* [l2 processing] Can now chose output handler (option `-l2-output` for `pysiral-l3proc.py`)
* [l3 processing] Updated `pysiral-l3proc` command line arguments (use `python pysiral-l3proc.py -help` for description)
* [l3 processing] Level-3 processor settings are now splitted in output definition and grid settings
* [l3 processing] Numerous changes to Level3Processing workflow
* [config] Update API of TimeRangeRequest

**Settings**
* [auxdata] Added definition for OSI-SAF 430 SIC data
* [l2 settings] Added definition (l2 settings & outputdef) for SPICES NRT ssh product
* [outputdef] Minor updated in SICCI-2 output definition


## Version 0.4.0 (9. June 2017)

**New Features**
* [visualization] Added mean filter for smooting in l3 maps (`create_l3s_map.py`: new option `-avfilt int:window_size`)
* [auxdata] Auxiliary data handler Added (default auxiliary handler needs type and id of data handler and gets all information from config file)
* [output] Level-2 output handler Added (reads new output definition files and automatically generates filenames and output folder for various output definitions)

**Bugfixes**
* [maptools] incorrect computation of pcolor grid corner coordinates

**Changes**
* [l2 processing] Option `-run-tag` has been made optional. If omitted, the run tag will be set to L2 settings id
* [l2 processing] Larger internal overhaul of the l2 processing approach (modularization of auxdata, l2 processing settings and output definitions)
* [l2 processing] new output definition files for generating output from the level-2 processor. Can be stacked. 
* [retracker] Added polyplane fit option of SICCI TFMRA retracker
* [uncertainties] computation of uncertainties improved
* [config] various changes to TimeRangeRequest
* [ssh] Added additional filters to SSASmoothedLinear


**Settings**
* [ccicdr] Finalized version 1 of ccicdr settings files (Envisat & CryoSat, North & South)
* [general] separated l2 settings and output definitions
* [general] settings files can now be put in subfolders for better organization
* [general] settings in subfolders named `obsolete` will not be automatically recognized by the setting lookup mechanismen


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


