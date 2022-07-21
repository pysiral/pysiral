# History of changes

## [0.10] - TBD

### Added
- high resolution land/ocean masks in the Level-1 pre-processor
- Support for 2021 version of the NSIDC sea ice region mask
- Support for SAMOSA+ retracker (external package) has been added (@djbwork)
- mypy type hint check have been added
- Support for ERS-1/ERS-2/Envisat ML based TFMRA thresholds (@sp-awi)
- southern hemisphere grids (12.5 & 25 km)
- filter class to detect and flag wave-based freeboard biases in the marginal ice zone
- release candidate 1 of the CryoSat-2 CCI v3.0 CDR
- distance to ice edge computed from ice concentration data
- support for Cryo-TEMPO phase 2.1 sea ice and polar ocean products (@djbwork)
- `ml` sub-module to `auxdata` for accessing trained models in the Level-2 processor
- type hints throughout the project 
- citation file (CITATION.cff)

### Changed
- Level-1 pre-processor workflow has been changed to allow processor items at different stages
- Data type for waveforms has been changed from float32 to float64 (@djbwork)
- Level-3 processor now only iterates over grid cells with data to increase execution speed
- Allow flexible epoch/units for Level-2 output
- Numerous changes to configuration files 
- TFMRA retracker may now exclude certain range bins for leading edge detection

### Fixed
- Fixed a SEGV when running cythonized TFMRA under macOS (@djbwork)
- Fixed an issue with a hard-coded nominal tracking bin for Envisat data
- Fixed incorrect date in L2P output filename when L2 data object contained a few data points of the previous day
- Fixed an issue of noise/FFT artefacts confusing the computation of leading edge width for pulse-limited waveforms
- 

### Deprecated
- class `pysiral.auxdata.GridTrackInterpol`. To be replaced by `pysiral.grid.GridImageExtract`
- class `pysiral.core.DefaultLoggingClass` 

### Removed
- TAI2UTC converter class has been removed in favor of the astropy.tyme module
- module `pysiral.logging` (replaced by use of loguru)
- Deprecated processor and output definition files

## [0.9.6] - 2021-10-12

### Added
- Added L1 prep-rocessor class for CryoSat-2 baseline-E L1b files 
- Including support for CryoSat-2 baseline D and E LRM data (requires waveform padding)
- Added modified L1 pre-processor class that patches the data with FES2014b after reading (@djbwork).
- ESA CryoTEMPO A001 processor and output definitions for both sea ice and polar ocean themes (@djbwork)
- AWI CryoSat-2 v2.4 processor and output definitions
- C3S SIThick CDR/ICDR processor and output definitions
- ESA CCI ERS TDS processor settings
- Auxiliary class for mean dynamic topography
- Parameter smoother class for the Level-2 processor pipeline
- Unfiltered sea level anomaly estimation  class for the Level-2 processor pipeline (@djbwork)
- Flag `instrument_mode` defined as `radar_mode` flag value + 1 (@djbwork)
- Support for relative and absolute orbit numbers in Level-2 output (@djbwork)
- Global package variable `__software_version__` intendend to contain the git hash value of the last commit as a unique identifier of software version (currently only meaningful when code is executed from a git repository)
- Added support for multi-dimensional data sets in the Level-2 output
- Level-2 processor item that generates multi-threshold freeboards
- Expanded set out output variables in CCI l2i output
- Dedicated TMFRA settings for Envisat and ERS-1/2
- Computation of leading edge quality indicator in the Level-1 pre-processor defined as fraction of cummulative power raise on leading edge to peak power (leq = 1 -> no artefacts on leading edge)
- Orbit flag (ascending/descending) to l1p and l2i output files
- Support for handling multiple l1p versions. (No action needed if only one version exists, for multiple the new `-l1p-version` option must be set when starting the Level-2 processor)
- Surface type classification class `pysiral.surface.ClassifierThresholdSurfaceType` for the Level-2 processor that allows to fully parametrize the list of variables and corresponding thresholds from the Level-2 processor configuration file. 
- Support for DTU21 mean sea surface
- CodeQL github action


### Changed
- C3S auxiliary data has been simplified to a single data handler. The version of the CDR/ICDR is to be specified in the options dictionary (See https://github.com/shendric/pysiral/discussions/57).
- Ensure that output flag attributes are always the correct data type (@djbwork)
- Allow different epochs to be used for the `time` variable in pysiral Level-2 output (@djbwork)
- Used the term `sea_ice_freeboard` instead of `freeboard` throughout output and code.
- Radar freeboard to sea ice freeboard conversion class now also stores the geometric correction as a Level-2 variable (`sgcorr:snow_geometric_correction`)
- Geophysical range correction class now outputs the total range correction to the Level-2 data container
- Convention for Level-3 status flag is now that value 0 means nominal acquision
- Convention for Level-3 quality flag is now that value 0 means nominal acquision
- Refactoring of the cTFMRA code (requires new set of options)
- Started to use python 3 type hints when adding/changing code.
- Computation of sigma0 takes now options for different radar mode types (and is no longer hard-coded to CryoSat-2 sensor specifications)
- OCOG power computation in Level-1 pre-processor is now based on physical units and no longer scaled waveform counts. 
- Changed behaviour of `pysiral.auxdata.sitype.OsiSafSIType` class to allow for a change between OSI-403-c and OSI-403-d (`confidence` variable replaced by `uncertainty`)
- L3 variable attributes of C3S CDR/ICDR v2.0 and AWI v2.4 now contain `coverage_content_type` attribute (ACDD convention)
- AWI l2i output is now using the CF compliant `platform` and `sensor` global attributes. Changes to l2p and l3 processor needed to be made as they rely on the previously used `source_mission_id` attribute. 
- Added `data_record_type` global attribute to AWI l2i output as replacement for the deprecated `source_timeliness` attribute. 
- various code style, documentation and logging improvements.

### Fixed
- An issue with newer python versions that used datetime of the `cftime` package (instead of `datetime`) was fixed
- A method to measure execution time (`time.clock()`) is deprecated and had to be replaced with `time.process_time()` in pysiral scripts. 
- Level-3 processor now checks data availability near the pole hole and does not rely on the grid cell center coordinate being beyond orbit coverage
- Computation of waveform noise level by the TFMRA is now longer done at a fixed set of range bins at the start of the waveform but has to be configured for each radar altimeter platform (@sp-awi)
- FFT artefacts in waveforms for older radar altimeter platforms were mistaken as the start of the leading edge by the TFMRA. The options now include a mandatory first valid index that specifies the earliest valid range bin for the leading edge (@sp-awi)
- Wrong Envisat ionospheric correction in CCI/C3S SGDR v3.0 Level-1 pre-processor leading to significant range noise. 
- Merged Warren99/AMSR2 snow climatology no longer causes an exception outside the October - April period


### Deprecated
- Class `pysiral.surface.SICCI2SurfaceType`. To be used in future: `pysiral.surface.ClassifierThresholdSurfaceType`

### Removed
- Deprecated output and processor settings files


## [0.9.5] - 2020-11-05

### Added
- Option to have multiple elevations from different retrackers to l2i output (range must be in l1p) 
- Option to reject CryoSat-2 l1b source file based on the status of orbit solution (coarse predicted orbits leads to freeboard elevation issues in CryoSat-2 baseline-D NRT data)

### Changed
- Updated the CCI+ v3.0-preview2 processor and output definitions

### Fixed
- Minor issue with detection first maximum in TFMRA (cython implementation)
- Degradation in CCI+ CryoSat-2 Antarctic data due to incorrect choice of geophysical range corrections
- Crash in Envisat SGDR input handler when constructing the lookup directory for SGDR data
- Crash of the southern hemisphere snow climatology happening on February 29. 

## [0.9.4] - 2020-10-19

### Fixed
- Some filters were not applied correctly in the sea-level anomaly estimation

## [0.9.3] - 2020-10-08

### Fixed
- Issues with Level-2 output format

## [0.9.2] - 2020-10-08

### Fixed
- Issues with Level-2 pre-processor

## [0.9.1] - 2020-10-08

### Fixed
- Added missing Level-2 pre-processor definitions for CryoSat-2 v2.3

## [0.9.0] - 2020-10-07

### Added
- Deblurring correction for ERS-1/2
- AWI CryoSat-2 v2.3 NRT and REP processor and output definitions
- CCI+ v3.0-preview2 processor and output definitions
- Snow density dependent geometric freeboard correction after Mallett et al., 2020
- dedicated `l2proc.procsteps` submodule
- tests for Level-2 processor definition files
- `core` module that includes basic pysiral functionality (introduced to avoid circular imports)
- SLA estimator that uses Gaussion Process Estimatation based on covariance and noise level of lead elevations (work in progress)

### Changed
- Format of this changelog now according to [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
- Merged W99/AMSR snow climatology is now daily instead of monthly to avoid unrealistic jumps at the beginning of a month
- Snow density parametrization in W99/AMSR snow climatology now follows Mallett et al. 2020
- Level-2 processor definition now contains a customizable list of algorithm steps (`procsteps`) instead of a fixed set to allow more customization
- Level-2 processor definition now contains a `metadata` group in level-2 processor definition files. Replaces older `id`, `mission`, `version_tag` and `hemisphere` entries
- Default output directory for level-2 data is defined as series of subfolders based on the `metadata` group in the processor file. 
- Logging system now based on loguru instead of logbook. New logging format.
- Split content of `surface_type` module between flag definitions and actual surface type classification processor steps
- Renamed `ssh` module to `sla` (more in line with L2 naming conventions)
- Made a distinction between dependancies (required packages) and installation reqirements (other github projects) in `setup.py` to avoid crash during `build_ext`
- Improving the way pysiral handles the different options for the location of the config files
- Performace improvements (approx. 30%) of Level-1 preprocessor by using `bottleneck.nanmean` instead of `np.nanmean`
- Added `fill_valid_sic_gaps` boolean flag to OSI-SAF sea ice type auxiliary classes to reduce data gaps near coasts
- Refactored code in `pysiral.sla`
- Dropped official python 3.5 support (github check does fail, but unrelated to pysiral code)

### Deprecated
- modules `iotools`, `roi`, `validator`  (will be removed shortly)

### Removed
- Errorhandler in pysiral.error
- Unused `bias` field in the L2 data container
- Deprecated Level-2 processor definitions (incompatible with the new Level-2 processor step definition)
- Deprecated output definitions
- `flag` module (functionality can be now found in `core`)
- `retracker.NoneRetracker` (no need to specify any longer)

### Fixed
- Some Level-2 processor definition files had incorrect names for ionospheric correction
- Issue with Level-3 array orientation (flip) when writing more than one output in the Level-3 processor
- Issue in the Level-2 processor that resulted in land flags being overwritten by the surface type classification algorithm
- Issue in CryoSat-2 baseline-D l1 pre-processor that resulted in all SARin waveforms being labeled as invalid due to a wrong interpretation of the measurement confidence flag
- Issue in the CryoSat-2 baseline-D l1 pre-processor that resulted in all files being ignored after a faulty netCDF file
- various minor bugs in settings files such as duplicate variable definitions/typos


## Version 0.8.0 (24. April 2020)

### Changes
- [python3] pysiral now supports python 3 (and python 3 only)
- [python3] replaced module `treedict` by `attrdict`
- [python3] replaced module `os` by `pathlib`
- [python3] various syntax changes with change to python 3
- [config] all configuration are no read once at the main import of the pysiral package and are accessible via the global variable pysiral.psrlcfg (replaces `pysiral.config.ConfigInfo`)
- [config] Various changes to the internal handling of auxiliary and altimeter platform metadata
- [config] Moved global variables in `pysiral.config` to config file
- [l1bpreproc] The old l1bpreproc Level-1 pre-processor has been removed 
- [test] Expanded tests of basic pysiral functionality (configuration management)

### Clean-up
- [pysiral] package wide removal of unused code bits
- [path] Removed `path` module
- [io_adapter] Removed `io_adapter` module
- [esa] Removed `esa` module
- [cryosat2] Removed modules associated with l1bpreproc
- [envisat] Removed modules associated with l1bpreproc
- [ers] Removed modules associated with l1bpreproc
- [icesat] Removed `icesat` module
- [legacy] Removed `legacy` module
- [maptools] Removed `maptools` module
- [orbit] Removed `orbit` module
- [l1bpreproc] Removed `l1bpreproc` module

### New Features
- [auxdata] Added functions to query all auxiliary data classes

### Bugfix
- [l3proc] flipping of y-coordinate of Level-3 output was handled badly for multiple outputs
 
 
## Version 0.7.3 (3. April 2020)

### Bugfix
- [config] Force remove time zones from data comparison when evaluationg requested time ranges

### Changes
- [auxdata] Added flexibility for the use of monthly or fixed thresholds in the SICCI surface type classifier class

### Settings
- [l2proc] Added ERS-1/2 v3.0-preview2 initial processor settings (surface type classification only at this point)
- [missions] Added ERS1 mission definition

### Technical
- [modules] replaced netCDF4.date2num with cftime.date2pynum throughout the package due to some weird inconsistencies of the return object type 


## Version 0.7.2.4 (13. February 2020)

### Bugfix
- [l1preproc] bugfix: local_machine_def_tag was not properly set to output handler

### Settings
- [l2proc] Added CCI+ output definitions
- [l2proc] Added l2 processor settings for AWI CryoSat-2 rep v2.2.1 (update of SIC ICDR)
- [l2proc] Added awi v2.2.1 l2p output definition
- [l3proc] Added awi v2.2.1 l3c output definition (normal and plus products)
- [l3proc] Update naming of cci level-3 processor and output definition files

## Version 0.7.2.3 (29. August 2019)

### Bugfix
- [envisat] final correction for range correction NaN issue in l1p files
- [sitype] corrected incorrect handling of time dimensions in the OSI-SAF CDR sea ice type handler 

### Settings
- [l2proc] Updated and corrected CCI+ settings for both platforms and hemispheres
- [output] Added reference for all public AWI v2.2 products

## Version 0.7.2.2 (25. August 2019)

### Bugfix
- [envisat] Additional fix for all-nan range correction variable in the original SGDR files.

## Version 0.7.2.1 (25. August 2019)

### Bugfix
- [envisat] Reinstate an exception for all-nan range correction variable in the original SGDR files.

## Version 0.7.2 (25. August 2019)

### New Features
- [l2proc] Allow retracker to add auxiliary variables
- [output] Variable attributes now can use template format

### Changes
- [auxdata] Modernized auxiliary data functions (sitype.OsiSafSITypeCDR, snow.ICDCSouthernClimatology)
- [l1preproc] variables in l1p/l1bdata now to not contain `scaling_factor` and `add_offset` to avoid issues with newer netCDF4 python module versions
- [l1preproc] added threshold for minimum size of l1 segments to be written to file
- [l2proc] cTFRMRA has now all features of SICCI2TfmraEnvisat (now obsolete)
- [l3proc] cTFMRA adds auxiliary variable `tfmra_threshold`

### Settings
- [l2proc] Added CCI+ Level-2 processor settings for Envisat (both hemispheres)
- [l2proc] Added CCI+ l2i data definition
- [l2proc] Added AWI l2p output definition

### Bugfix
- [l1bdata] An issue was resolved that lead to passing wrong datatypes to cythonized functions in cTFMRA
- [envisat] Issue with NaNs in Envisat l1p range correction output resolved.

## Version 0.7.1 (20. August 2019)

### New Features
- [l3proc] Compute temporal statistics of data in grid cell
- [l3proc] Grid classifiers for different surface_types

### Changes
- [l3proc] Major overhaul of Level-3 processing pipeline (requires new processor settings)
- [l3proc] Allow multiple output files

### Settings
- [l2proc] Added new Level-2 processor definition(s) for AWI CryoSat-2 v2.2
- [l2proc] Splitted Level-2 processor definition(s) for AWI CryoSat-2 v2.2 
- [l3proc] Added new Level-3 AWI CryoSat-2 v2.2 processor definitions for nrt and rep (different auxiliary data) 
- [l3proc] Added new Level-3 AWI CryoSat-2 v2.2 output definitions for regular and extended (plus) output
- [l3proc] Variable renaming and clean-up of Level-3 AWI CryoSat-2 v2.2 products

### Bugfixes
- [l1preproc] Improve stability of ERS Level-1 preprocessor for cases with incorrect global attributes in REAPER files

## Version 0.7.0 (7. August 2019)

### New Features
- [l1preproc] Added a new Level-1 preprocessor system (l1preproc) that will gradually replace the original one (l1bpreproc)
- [l1preproc] Add support for CryoSat-2 Baseline-D, Envisat SGDR v3.0, ERS-1/2 REAPER and Sentinel-3A/B for l1preproc
- [doc] Add support for automatic documentation building using SPHINX (-> https://pysiral.readthedocs.io/en/latest/)
- [config] Add a stop watch to track time

### Changes
- [l1bdata] Add antenna pitch, roll, heading to the time orbit group

### Settings
- [l1preproc] Added CryoSat-2 Baseline-D
- [l1preproc] Added Envisat SGDR v3.0
- [l1preproc] Added ERS-1/2 REAPER 
- [l1preproc] Added Sentinel-3A/B for l1preproc
- [l2proc] Added beta for CryoSat-2 v2.2

## Version 0.6.6 (14. November 2018)

### New Features
- [l2proc] Added direct transfer of l1p variables to the l2 object via processor configuration files
- [rio] Added rio auxiliary data type (rebase from `fmi` branch)
- [icechart] Added icechart auxiliary data type (rebase from `fmi` branch)

### Changes
- [l2proc] Auxiliary data handlers now have a unique id, allowing multiple handlers of the same type (e.g. icechart.canada, icechart.aari) to be registered. Old state was one handler per type.

### Settings
- [fmi] Added FMI l2i output (rebase from `fmi` branch)


## Version 0.6.5 (9. November 2018)

### New Features
- [config] directory for pysiral configuration can now be specified to allow running multiple instances with different configuration by the same user on the same machine (fully encapsulated installation)


## Version 0.6.4.5 (23. October 2018)

### Changes
- [config] set default value of overwrite protection to False
- [config] renamed period name `default_week` to `weekly`

## Version 0.6.4.4 (19. October 2018)

### Changes
- [setup] Added userhome pysiral update script to python package

## Version 0.6.4.3 (19. October 2018)

### Changes
- [l2proc] Rollback of CryoSat-2 AWI v2.1 surface type classification algorithm to the CCI one. RickerTC2014 was causing issues of SARin areas. 

## Version 0.6.4.2 (18. October 2018)

### Changes
- [l3proc] Added l2i pre-filtering to the AWI CryoSat-2 v2.1 Level-3 processor settings. To confusing to have radar freeboard over leads

## Version 0.6.4.1 (17. October 2018)

### Changes
- [l3proc] Allow masking of non-float Level-3 parameter (e.g. integer flags, which will be set to -1)

### Bugfixes
- [l3proc] Level-3 processor was not adapted for name change of l2 variable (`ice_density` -> `sea_ice_density`)
- [l3proc] Incorrect source variable naming for uncertainties of `radar_freeboard`, `freeboard` & `sea_ice_thickness` in AWI v2.1 l3 output definition 
- [snow] bugfix: missing message string caused crash of data handler (see #19)

## Version 0.6.4 (15. October 2018)

### New Features
- [mss] Added support for DTU18
- [snow] Added support for new (finalized) merged Warren/IUP AMSR2 snow climatology
- [region] Added new `region` auxiliary data type for Level-2 processor
- [l2proc] Added modified NSIDC region definition to Level-2 and higher level products
- [l3proc] Added label of the period (`monthly`, `weekly`) to the grid output subfolders
- [scripts] Added a script (`pysiral/bin/psrl_update_userhome_cfg.py`) that will update the pysiral configuration in the user home with the definition files of the pysiral package (excluding `local_machine_def.yaml`)

### Changes
- [auxdata] Major overhaul of auxiliary data ingestion engine (main part of this release). A custom list of handlers now replaces a static set of defined handlers (`mss`, `sic`, `sitype`, `snow`)
- [auxdata] All auxiliary data types are now moved the the sub-module: `pysiral.auxdata`
- [auxdata] Created a common interface for all auxiliary data classes. Adapted all existing auxiliary data handlers
- [general] Various documentation improvements and code cleanup

### Bugfixes
- [l2data] l2data.Level2Data._PARAMETER_CATALOG was a mutable class variable, and it could be unintentionaly changed for all instances 

### Settings
- [general] New format of Level-2 auxiliary definition
- [awi] Added AWI v2p1 Level-2 & Level-3 processor settings
- [awi] Added AWI v2p1 Level-2 & Level-3 output definition files

### Known Issues
- [settings] Not all Level-2 settings adapted to new auxiliary data engine


## Version 0.6.3 (26. August 2018)

### New Features
- [snow] Added support for Merged Warren99 / AMSR2 snow depth climatology (`snow.Warren99AMSR2Clim`)
- [auxdata] Added class for fast extraction of along-track data from a grid (`auxdata.GridTrackInterpol`)

### Changes
- [auxdata] Moved properties and methods from subclasses in the `auxdata.AuxdataBaseClass` superclass with the intention to reduce code for new auxiliary subclasses. 
- [sic] The class for OSI-401 products now has the option `fill_pole_hole` to avoid data loss around 88N

### Bugfixes
- [snow] Ice-type dependent scaling of Warren99 snow depth only valid for factor of 0.5
- [l2proc] l2proc crashes

### Settings
- [awi] Added initial version of the AWI v2p1 Level-2 processor settings 


## Version 0.6.2 (20. August 2018)

### New Features
- [github] Added bug report / pull request templates
- [github] Added contributing information 

### Changes

### Bugfixes
- [snow] Ice-type dependent scaling of Warren99 snow depth only valid for factor of 0.5
- [l2proc] l2proc crashes

### Settings 

## Version 0.6.1 (13. August 2018)

### New Features
- [setup] Allow pysiral to be installed using pip
- [config] Change of settings structure. Config files are categorized in `type` (proc, output, grid)" and `data_level` (l1p, l2i, l2p, l3, None for grid)
- [pysiral] Move to github

### Changes
- [config] pysiral configuration files are now in the user home directory in the `.pysiral-cfg` subfolder. The existence of the folder is checked upon calling pysiral and automatically created if it does not exist. 
- [config] Removing obsolete config files from the code repository
- [config] Improved error message for missing  `local_machine_def.yaml`
- [collection] Module `collection` removed from pysiral (moved to
 `pysiral-product-tools`)
- [catalog] Module `catalog` removed from pysiral (moved to `pysiral-product-tools`)
- [visualization] Inactive module `visualization` removed from pysiral
- [version] pysiral version tag moved to `pysiral/VERSION`

### Bugfixes

### Settings 

## Version 0.6.0 (1 August 2018)

### New Features
- [auxdata] Add the sic data records (cdr/icdr) for c3s
- [config] Added orbit inclination dict for level 3 status flag
- [catalog] Added a catalog module that inventarizes all pysiral product (L2P, L3C) in a given directory
- [collection] Added a collection module to create and manage time series of pysiral products
- [filter] Added simple smoothing filter for 2D Arrays with the option to preserve gaps
- [griddef] Added 12.5 km EASE2 grid for northern hemisphere
- [l2proc] allow multiple subfolders as run-tag (e.g. `-run-tag folder_1/folder_2/.../folder_n`)
- [l3proc] Allow to specifically set the flag values for the status flag in the l3 settings file
- [l3proc] Added additional masking method for the l2i stack (here for testing of ascending/descending grids)
- [l3proc] Add option to filter ascending/descending orbit parts before gridding
- [l3proc] Add additional L3 parameter: `negative_thickness_fraction`
- [l3proc] Added computation of level 3 uncertainties. (Error propagation for radar freeboard, freeboard, thickness in each cell based on stack and average values)
- [l3proc] Added property `data_record_type` to Level-3 processor to specify cdr or icdr output
- [product] Enable the writing of doi's into netcdf global attributes. doi must be supplied as a command line argument for l2preproc, l3proc
- [settings] added `time_dim_is_unlimited` grid option attribute
- [settings] Added C3S Envisat Level-2 algorithm settings
- [snow] Added evaluate method to the warren snow climatology (return snow parameters for a given set on lons, lats, and month)
- [visualization] Added ParameterScatterPlot (work in progress)

### Changes
- [l2data] renamed level-2 dimension from `n_records` to `time`
- [l3proc] Allow time dimension in l3 output to be UNLIMITED
- [l3proc] Improved error logging
- [mask] Added option to flip a mask to l3 mask
- [mask] Always return lon/lat variables for each mask
- [output] Allowed to use sanitized version of sensors names in filenames (e.g. RA-2 -> RA2)
- [output] Two new options in l3 output definition: `root.grid_options.flip_yc` (`True` or `False`) and `root.grid_options.time_dim_is_unlimited` (`True` or `False`)
- [settings] various improvements in netcdf attributes
- [settings] Increased threshold filter size for Envisat sea ice thickness to account for larger range noise
- [settings] Create level-3 settings for SICCI-2 CDR with other changes status flag values
- various code documentation improvements

### Bugfixes
- [cython] usage of full module name to prevent issues with compilers
- [l2proc] the freeboard filter was wrongly determined from the radar freeboard. Negative freeboards were over-filtered, since radar freeboard is always lower.
- [l3proc] Invalid check if `time_dim_is_unlimited` grid options attribute exists
- [l3proc] wrong config filter tag for orbit filter
- [maptools] pcolor grid (grid corner coordinates) not correctly computed
- [output] netcdf attribute time_bnds lacked the necessary time dimension
- [output] Umlaut in global attributes crashed the attribute template engine
- [output] template engine crashed if attribute is not a string
- [sentinel3a] output not written if `-hemisphere` option was set in the l1 pre-processor
- [settings] removed fixed value for orbit_inclination tag in l3 settings


## Version 0.5.0 (9 February 2018) 

### Bugfixes
- [auxdata] Longnames of sitype in `auxdata_def.yaml` were not recognized properly (incorrect key)
- [datahandler] Catch exception caused by feeding empty directory to datahandler
- [l1bpreproc] catch a rare exception when the l1b input file is just 1 data point
- [l2data] radar_mode information was lost during the Level-2 processing
- [l2proc] bugfix: Improper handling of pysiral-l2proc call for specific list of l1b input files
- [l2preproc] Catch exception if error during parsing a l2i file
- [l2preproc] Catch exception if l2i file(s) exist for a single day but none contains valid freeboard data
- [l2preproc] Fixed a bug that caused all uncertainties in l2p products to be zeros
- [l3proc] Catch exception if error during parsing a l2i file
- [l3proc] Fixed a bug in the generation of lat/lon grid coordinates that caused the grid to be incorrect by 50% of the grid size at the outer grid limits
- [sentinel3] Fixed several issues in the Sentinel-3A preprocessor
- [sitype] sitype api method did not return value, uncertainty, msg
- [general] Fixed a few of incorrect variable names

### New Features
- [catalog] Added catalog module for cataloging/querying pysiral product files.(Currently only working for l2p files, that will be extend to all data levels)
- [config] Added default_week (Monday - Sunday) as period in TimeRangeRequest
- [config] Added duration string in ISO format to TimeRangeRequest
- [config] Added monthly, weekly and daily iterators to TimeRangeRequest
- [config] Added isodate labels to TimeRangeRequest
- [datahandler] Added capability to find all l2i files for an arbitrary time_range
- [general] Implement timeliness attribute throughout the processing chain
- [general] pysiral-wide disabling of python warnings
- [general] keep internal record of official mission names (platform)
- [grid] new grid definition property: pyresample area definition
- [icesat] Added icesat module for converting GLAH13 data into a pysiral compliant data format
- [icesat] Added ICESat pre-processor
- [l1dbata] Added method to fill gaps in Level-1 data 
- [l1dbata] Added method to split l1bdata objects at discontinuities 
- [l2proc] Added class for snow freeboard class for converting radar freeboard to freeboard (pysiral.frb.SnowFreeboardAssumption)
- [l2proc] Added surface type classificator for ICESat based on (Khvorostovsky, K. and Rampal, P.: On retrieving sea ice freeboard from ...)
- [l2data] Added method to init l2data object from l2i netcdf file
- [l3proc] Added option to export grid description parameters as variable in netcdf output
- [l3proc] Added time bound property for Level-3 netcdf output
- [l3proc] Added default time and grid parameters to Level-3 netcdf output
- [l3proc] Allow flipping of Level-3 variable in netcdf output (to be consistent with osisaf cdr)
- [mask] new module handling masks
- [mask] Add mechanic to compute fixed masks for each l3 grid definition 
- [mask] Add l3 mask for Arctic area in which Warren99 is considered valid
- [mask] Add land/sea mask with is_land, is_sea & is_mixed flag
- [sentinel3] Add parser for EUMETSAT L2_WAT product
- [settings] New attribute template {uuid} to assign universally unique identifier (uuid) tracking id for l2/l3 products
- [settings] Attribute template {source_hemisphere} has now a select option to select text based on the hemisphere. (Usage: {source_mission_id:select;id1:Text id1;id2:Text id2;...}`)
- [settings] Attribute template {source_mission_id} has now a select option to select text based on the mission id. (Usage: {source_hemisphere:select;Text north;Text south}`)
- [sitype] Added sitype class for osi-saf sitype cdr products

### Changes
- [l2data] Renamed `timestamp` variable to `time` due to variable naming conventions (backward compability with older files exist)
- [l2data] Added version tag of the algorithm (required in all l2 settings files)
- [l3proc] Expanded the list of possible Level-3 parameters in preparation to status/quality flags
- [l3proc] Expanded the options for default Level-3 processing with external masks and post-processing
- [l3proc] Allowed to load external masks in Level-3 processor
- [l3proc] Allowed computation of quality indicator flags
- [l3proc] Allowed computation of status flags
- [l3proc] Updated Level-3 settings with improved variable initializuation (init value, dtype)
- [settings] Global attributes in l2/l3 products now follow the order prescribed in the output definition files (was alphabetic)
- [settings] Attribute template options are now separated by semikolon (was comma)
- [settings] Time coverage start/end in l2p files are now always defined by requested data period (full day) and not actually available data points
- [settings] Filenaming in settings files is now depends on data period (e.g. different filenaming for weekly or monthly l3c products)
- [general] Set the creation time for products once at the initialization of the python class (was multiple calls of datetime.now() during output generation)
- [griddef] Added grid name in grid definitions
- [general] Countless code documentation, style and logging improvements

### Settings
- [l2proc] Added Sentinel-3 l2proc setting files
- [l2proc] Added Level-2 algorithm settings recipe for C3S cdr v1
- [l3proc] Added default Level-3 settings for southern hemisphere without Warren99 is_valid check. Needs to be selected manually with the `-l3-settings` option. 
- [output] Added L3C output definition for C3S
- [output] Added l2p/l3c definition for AWI product in version 2.0
- [output] ccicdr output definitions (l2p/l3c) now follow CCI data standard requirements

## Version 0.4.7 (17 August 2017) [bugfix release]

### Bugfixes
- [settings] L3C filename was incorrectly written (request period start time instead of iteration period start time)


## Version 0.4.6 (16 August 2017)

### New Features
- [l2proc] Added input validation for l1bdata files (time coverage and versions)

### Settings
- [l3proc] Changed SICCI-2 L3C global attributes towards pysiral and Data Discovery standards


## Version 0.4.5 (24 July 2017) [bugfix release]

### Bugfixes
- [settings] Default uncertainty value for sitype south was missing


## Version 0.4.4 (24 July 2017) [bugfix release]

### Bugfixes
- [settings] Wrong SICCI-2 L3C version
- [l2proc] ICDC sitype class crashed on missing file


## Version 0.4.3 (17 July 2017) [bugfix release]

### Bugfixes
- [settings] Invalid default l2i settings solved


## Version 0.4.2 (7 July 2017)

### New Features
- [sitype] MYI concentration uncertainty support Added
- [l2 processing] Added l1b prefilter to the l2proc workflow
- [l2 processing] OsiSafSIC now allows auto switch between products (e.g. osisaf-409 until 2015-04-15 and osisaf-430 after)

### Bugfixes
- [l3 processing] #396: sea ice type uncertainty not filtered by freeboard mask in default l3 settings
- [output] #397: valid_fraction missing in SICCI-2 L3C output
- [output] #399: attribute tag `utcnow` not working

### Changes
- [sit] computation of sea ice density uncertainty updated
- [surface_type] Merged SICCI2Envisat and SICCI2CryoSat2 surface type classes

### Settings
- [l2 processing] Added Envisat backscatter drift correction as l1b prefilter
- [l2 processing] osisafcdr (merged osisaf-409 & osisaf-430) sic for all ccicdr settings


## Version 0.4.1 (3. July 2017)

### New Features
- [grid] New `grid` module for grid related tasks (work in progress)
- [l2 data] Added `sea_surface_height` property (ssh = mss + ssa)
- [auxdata] Added a `NoneHandler` (dummy handler that only returns NaN's)
- [l3 processing] Added masking of l2i data before gridding

### Bugfixes
- [l3 processing] Level-3 processor crash when no l2i data was available for given iteration

### Changes
- [l2 processing] Updated uncertainty description
- [l2 processing] Can now chose output handler (option `-l2-output` for `pysiral-l3proc.py`)
- [l3 processing] Updated `pysiral-l3proc` command line arguments (use `python pysiral-l3proc.py -help` for description)
- [l3 processing] Level-3 processor settings are now splitted in output definition and grid settings
- [l3 processing] Numerous changes to Level3Processing workflow
- [config] Update API of TimeRangeRequest

### Settings
- [auxdata] Added definition for OSI-SAF 430 SIC data
- [l2 settings] Added definition (l2 settings & outputdef) for SPICES NRT ssh product
- [outputdef] Minor updated in SICCI-2 output definition


## Version 0.4.0 (9. June 2017)

### New Features
- [visualization] Added mean filter for smooting in l3 maps (`create_l3s_map.py`: new option `-avfilt int:window_size`)
- [auxdata] Auxiliary data handler Added (default auxiliary handler needs type and id of data handler and gets all information from config file)
- [output] Level-2 output handler Added (reads new output definition files and automatically generates filenames and output folder for various output definitions)

### Bugfixes
- [maptools] incorrect computation of pcolor grid corner coordinates

### Changes
- [l2 processing] Option `-run-tag` has been made optional. If omitted, the run tag will be set to L2 settings id
- [l2 processing] Larger internal overhaul of the l2 processing approach (modularization of auxdata, l2 processing settings and output definitions)
- [l2 processing] new output definition files for generating output from the level-2 processor. Can be stacked. 
- [retracker] Added polyplane fit option of SICCI TFMRA retracker
- [uncertainties] computation of uncertainties improved
- [config] various changes to TimeRangeRequest
- [ssh] Added additional filters to SSASmoothedLinear


### Settings
- [ccicdr] Finalized version 1 of ccicdr settings files (Envisat & CryoSat, North & South)
- [general] separated l2 settings and output definitions
- [general] settings files can now be put in subfolders for better organization
- [general] settings in subfolders named `obsolete` will not be automatically recognized by the setting lookup mechanismen


## Version 0.3.0 (10. March 2017)

### New Features
- [l2 processing] Added SICCI-2 retracker class for Envisat (retracker.SICC2TfmraEnvisat)
- [l2 processing] Added SICCI-2 surface type class for Envisat and CryoSat-2 (surface_type.SICCI2Envisat, surface_type.SICCI2Cryosat2)
- [l2 processing] Added support for NASA-Team based MYI concentrations produced by Integrated Climate Data Center (ICDC) (-> sitype.ICDCNasaTeam)
- [l2 processing] Added support for Antarctic daily snow climatology by Integrated Climate Data Center (ICDC) (-> snow.ICDCSouthernClimatology)
- [l2 processing] Added error catching and reporting for algorithm classes. Run statistics, software environment parameters and a breakdown of encountered error is logged and written to `pysiral-l2proc-summary.txt` in the l2i export directory
- [l2 processing] Added capability to pass a list of input files using the `-l1b-files` argument (mutually exclusive with `-start` and `-stop`)
- [l2 processing] surface type classification algorithm can now access radar mode flag
- [l2 processing] retracker algorithm can now retrieve parameters from l1b object
- [l3 processing] Added surface type statistics (level-3) parameter (`n_total_waveforms`, `n_valid_waveforms`, `valid_fraction`, `lead_fraction`, `ice_fraction`)
- [l3 processing] Added sea ice concentration masking (grid parameter to nan if `sic < 5` or `sic == nan`)
- [l3 processing] Added netCDF variable attributes to L3s output
- [tools] Added customization options for the colorbar generator
- [scripts] Added bash script to download, extract and sort CryoSat-2 NRT files (requires username & password for ESA ftp)

### Bugfixes
- [config] Bug #280: month exclusion ignored for the month given by -start option in the l1bpreproc, l2proc and l3proc calls leading to incorrect start and stop dates of the processors

### Changes
- [l2 processing] reorganization of l2 processing scheme to allow computation of uncertainties
- [l2 processing] Improved error handling and reporting
- [l3 processing] Improved logging and documentation
- [l3 processing] Additional field `sea_ice_concentration_mask_targets` in level-3 output definition required
- [visualization] Minor changes to colormaps
- [visualization] Added sea ice concentration to map background
- [visualization] Added capability to create Antarctic maps

### Settings
- [level-2] Added SICCI-2 prototype 002 settings for Arctic and Antarctic
- [level-2] Added SICCI-2 prototype 003 settings for Arctic and Antarctic


## Version 0.2.0 (19. September 2016)

### New Features
- [l1b pre-processing] time range for can be specified for individual days
- [l1b pre-processing] pre-processor can be set to specific hemispheres
- [l1b pre-processing] obsolete l1bdata files can be automatically removed (`--remove-old` option to pysiral-l1bpreproc)
- [l1b pre-processing] new classifier parameter for CryoSat-2, Envisat and ERS
- [l2 processing] time range for can be specified for individual days
- [l2 processing] obsolete l2 output files can be automatically removed (`--no-overwrite-protection` & `--remove-old` options to pysiral-l2proc)
- [l2 processing] l2i output files log the settings of the level2 processor in the global attributes
- [l2 processing] l2i output variables with proper netCDF attributes
- [retracker] cythonized implementation of the TFRMA (3 times faster than pure python)
- [tools] `l1bdata_report.py`: creates a pdf summary of l1bdata files
- [tools] `zipl2month.py` creates monthly archives of l2i orbit files

### Bugfixes

- [l1b pre-processing] TAI to UTC timestamp conversion for CryoSat-2 incorrect
- [l1b pre-processing] Wrong values for open ocean percentage when subsetting and merging l1b segments
- [l1b pre-processing] pysiral version tag was not properly passed to l1bdata netcdf files

### Changes

- [l1b pre-processing] moved hard coded settings from python parser to mission config file (Envisat, ERS, Sentinel-3)
- [l1b pre-processing] All geophysical corrections are now piped into l1bdata files instead of selection in the mission configuration. Requires changes in the l2 settings files.
- [l2 processing] Internal handling of processor settings


