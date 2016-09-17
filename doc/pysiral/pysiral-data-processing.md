# Pysiral Data Processing

## Level-1 Pre-Processing

The main purpose of the L1 pre-processor is to create Arctic and Antarctic subsets in a unified data format (netCDF) from the mission specific input files. The pre-processing main routine is placed in `bin/pysiral-l1bpreproc.py` of the pysiral repository. 

### Usage

    python pysiral-l1bpreproc.py [options]

The output will be written in the `l1bdata` folder in `local_machine_def.yaml` for the specific mission, data version, hemisphere and month. 

### Options

Argument | Required | Parameter | Description
--- | --- | --- | --- | ---
**-h** | optional | *None* | Display manual
**-mission** *mission_id* | yes |string (e.g. `-mission cryosat2`) |  list of pysiral recognized mission names (see `missions` in `config\mission_def.yaml`)
**-start** *start_date* | yes | integer list | year, month [, day] definining the start of the data period. If day is ommittted, always the first day of the month will be used
**-stop** *stop_date* | yes | integer list | year, month [, day] definining the end of the data period. If day is ommittted, always the last day of the month will be used
**-hemisphere** *hemisphere_id* | optional | string | select the hemisphere for the generation of l1bdata files , either `north`, `south` or `global` (default)
**-exclude-month** *month_list*  | optional | integer list | list of month (from 01 - 12) that will be skipped)
**-input-version** *version_id* | optional | string  | version name of the mission specific data (default: `default`). Any additional version must have a corresponding entry in `local_machine_def.yaml` for the particular mission
**--remove-old** | optional | *None* | delete all l1bdata files in output directory before pre-processing. Triggers manual confirmation prompt
**--no-critical-prompt** | optional | *None* | skip the manual confirmation prompt for the `--remove-old` option


### Cookbook

1) Pre-processing CryoSat-2 data for on Arctic winter season (2015-10-01 till 2016-04-30)

     python pysiral-l1bpreproc.py -mission cryosat2 -start 2015 10 -stop 2016 04 -hemisphere north

2) Pre-processing ERS-2 data for 1 week (north and south)

     python pysiral-l1bpreproc.py -mission ers2 -start 2003 03 01 -stop 2003 03 07

3) Batch pre-processing CryoSat-2 winter data over several years

     python pysiral-l1bpreproc.py cryosat2 -start 2010 10 -stop 2016 04 -hemisphere north -exclude-month 05 06 07 08 09

4) Pre-Processing CryoSat-2 baseline-B data

     python pysiral-l1bpreproc.py cryosat2 -start 2014 03 -stop 2014 04 -hemisphere north -input-version baseline_b

Needs the corresponding entry besides `default` in `local_machine_def.yaml`, e.g.:

    l1b_repository:
    
        # CryoSat-2
        cryosat2:
            default: # default is baseline-c
                l1bdata: 'E:\altim\data\altimetry\cryosat2\baseline-c\l1bdata'
                sar: 'E:\altim\data\altimetry\cryosat2\baseline-c\SIR_SAR_L1'
                sin: 'E:\altim\data\altimetry\cryosat2\baseline-c\SIR_SIN_L1'     
            baseline_b:
                l1bdata: 'E:\altim\data\altimetry\cryosat2\baseline-b\l1bdata'
                sar: 'E:\altim\data\altimetry\cryosat2\baseline-b\SIR_SAR_L1'
                sin: 'E:\altim\data\altimetry\cryosat2\baseline-b\SIR_SIN_L1'


___



## Level-2 Processing

The level-2 processor computes geophysical parameters from along-track (pre-processed) along-track l1b orbits. The output is saved in netCDF files that corresponds in number of records to the intput l1b netCDF input files. The main script is `bin/pysiral-l2proc.py` of the pysiral repository. 

### Usage

    python pysiral-l2proc.py [options]

### Options

Argument | Required | Parameter | Description
--- | --- | --- | --- | ---
**-h** | optional | *None* | Display manual
**-l2-settings** *settings_definition* | yes | string | level-2 processor settings yaml file (either a the basename of pre-defined level-2 settings (`in settings/l2`, or a full file path). The settings file contains additional information, e.g. mission id, hemisphere, choice of auxiliary data sets ... See below for a more detailed documentation of the level-2 processor settings. 
**-run-tag** *id* | yes | string | An id of the specific processor run that is used to generate the output directory structure
**-start** *start_date* | yes | integer list | year, month [, day] definining the start of the data period. If day is ommittted, always the first day of the month will be used
**-stop** *stop_date* | yes | integer list | year, month [, day] definining the end of the data period. If day is ommittted, always the last day of the month will be used
**-exclude-month** *month_list*  | optional | integer list | list of month (from 01 - 12) that will be skipped)
**-input-version** *version_id* | optional | string  | version name of the mission specific data (default: `default`). Any additional version must have a corresponding entry in `local_machine_def.yaml` for the particular mission
**--no-overwrite-protection** | optional | *None* | Flag that prevents the generation of a directory with the current time in order to produce a unique output directory. If this flag is set and the run tag has been used before, output files from different processor settings or version will end up in the same output directory
**--remove-old** | optional | *None* | delete all l2 files for all output types in their specific directory before processing. Triggers manual confirmation prompt. Has no effect if `--no-overwrite-protection` is not set
**--no-critical-prompt** | optional | *None* | skip the manual confirmation prompt for the `--remove-old` option

### Level-2 processor settings

The mandatory settings file needs to be in YAML format and contain all necessary information on the 

1. mission
- region of interest
- auxiliary data sets (mean sea surface, sea ice concentration, sea ice type, snow depth)
- geophysical range corrections
- surface type classification algorithm
- retracker for different surface types
- sea surface anomaly algorithm
- freeboard retrieval algorithm
- freeboard to thickness conversion
- filters for freeboard and thickness (individual points)
- validation methods (used to classify entire orbit segment as valid/invalid)
- postprocessing methods
- generation of output files (can be multiple files per orbit segment)

For each of this items a specific python class can be selected and initialized with parameter options from the l2 settings file, thus enabling a fully customizable sea ice radar altimeter processor. The general structure of the settings file is defined as:

```python
# Level 2 processor settings are mission specific
mission:
    id:  # mission Id

# Regions Settings (for file selection and potential subsetting
roi:
    pyclass:         # name of region selection class (roi module)
    hemisphere: 
    options:
        name: value
        ...

# Settings of Level-2 orbit processing
level2:

    # Sources of ancillary datasets
    # (the tag "name" links to the corresponding tag in config/auxdata.yaml)
    auxdata:
        
        mss:
            name:     # mss data source id (-> config/auxdata_def.yaml) 
            options:
                name: value
                ... 
        sic:
            name:     # sic data source id (-> config/auxdata_def.yaml)  
            options:
                name: value
                ... 
        sitype:
            name:     # sitype data source id (-> config/auxdata_def.yaml)  
            options:
                name: value
                ... 
        snow: 
            name:     # snow data source id (-> config/auxdata_def.yaml)  
            options:
                name: value
                ... 
   
    # geophysical corrections applied to the l1b range window
    corrections:
        -  # list of range corrections 
        
    # Surface type classification algorithm
    surface_type: 
        pyclass:    # name of surface type classificator class (surface_type module)
        options:    # options for different surface types
            ocean:
                name: value
                ...
            lead: 
                name: value
                ...
            sea_ice: 
                name: value
                ...
                
    # Retracking algorithm dependent on surface type
    retracker: 
        ocean: 
            pyclass:    # name of ocean retracker class (retracker module)
            options:
                name: value
                ...
        lead: 
            pyclass:    # name of lead retracker class (retracker module) 
            options:
                name: value
                ...
        sea_ice: 
            pyclass:    # name of sea ice retracker class (retracker module) 
            options:
                name: value
                ... 
                
    # Algorithm for instantaneos sea surface height (mss + ssa) and radar freeboard 
    ssa: 
        pyclass:    # name of instanteneous sea surface height estimater class (mss module)
        options:      
            name: value
            ...  
            
    # Algorithm for converting radar freeboard into freeboard
    frb: 
        pyclass:    # name of freeboard (not radar freeboard) estimator class (frb module)
        options: 
            name: value
            ...  
            
    # Algorithm for getting sea ice thickness from other l2 parameters
    sit: 
        pyclass:    # name of freeboard to thickness converter class (sit module)
        options: 
            name: value
            ...
            
    # List of filter at different stages of l2 processing 
    # (can be list of filters)
    filter: 
    
        # Filters after freeboard computation
        freeboard: 
            frb_valid_range: 
                pyclass:    # name of freeboard filter class (filter module)
                options:
                    name: value
                    ...   
                    
        # Filters after thickness computation
        thickness: 
            sit_valid_range: 
                pyclass:    # name of thickness filter class (filter module)
                options: 
                    name: value
                    ... 
                    
    # Tests if l1b orbit file is valid
    validator: 
    
        # Tests (Can be several) after surface type classification
        surface_type: 
            n_leads: 
                pyclass:    # name of lead number validator class (validator module) 
                options:
                    name: value
                    ...  
            
    # Post-Processing (tbd if needed)            
    post_processing: null
    
    # Definition of output files of the l2 orbit datasets
    # (can be several if needed)
    output: 
        l2i:                  # example output id
            pyclass:          # L2 output class (output module)
            path: null        # Will be added during processing
            options: 
                subfolders:   # subfolder structure after main_product_dir/run_tag/[unique_run_id]/subfolder[1]/subfolder[2]/
                    - year 
                    - month
                parameter:    # list of l2 parameter names 
                    - ...          
```

The key mission settings files are under revision control in the subfolders `settings/l2` in the pysiral repository.

### Processor output

The type(s) of output files need to be defined in the settings file (`level2.output.l2_output_id.XXX`) including their subfolder structure. The full path of the l2 files is then based on these information, the options of the main processor call and the main product folder defined in in `local_machine_def.yaml`.
This can be summarized as: 

    product_repository/run_tag/[unique_dir]/output_id/subfolders[1]/subfolders[2]...

key| defined by
--- | --- | 
**`product_repository`** | `key product_repository` in `local_machine_def.yaml`
**`run_tag`** | from l2 processor argument `-run-tag`
**`unique_dir`** | automatic generated directory name (YYYMMDDTHHMiSS). Will not be ommitted if `-no-overwrite-protection` is set
**`subfolders`** | defined in the l2 settings file, typically the year and the month of the input data files

### Auxiliary data set definitions

The Level-2 processor requires additional datasets or parametrisations for the mean sea surface, sea ice concentration, type and snow depth & density. There are two locations were these datasets need to be defined

1. `config/auxdata_def.yaml`: default options and name of the corresponding pysiral handler class. 
2. `local_machine_def.yaml`: local directory name for the different auxiliary data sets. 

The link between both configuration files is made by the id of the auxiliary data set (e.g. osisaf for sic). The content of the l2 settings file only need this id, though default option can be overriden if necessary.  


### Cookbook

1) Processing CryoSat-2 data for on Arctic winter season (2015-10-01 till 2016-04-30) with the default CryoSat-2 AWI configuration. (Note: the option `-l2-settings cryosat2_awi_north` refers the settings file `settings/l2/cryosat2_awi_north.yaml`)

     python pysiral-l2proc.py -l2-settings cryosat2_awi_north -run-tag some_run -start 2015 10 -stop 2016 04 

2) Batch processing CryoSat-2 data for several Arctic winters with the default CryoSat-2 AWI configuration. 

     python pysiral-l2proc.py -l2-settings cryosat2_awi_north -run-tag some_long__run -start 2011 10 -stop 2016 04 -exclude_month 5 6 7 8 9

3) Processing from a specific l1b data version

     python pysiral-l2proc.py -l2-settings cryosat2_awi_north -run-tag old_stuff -input-version baseline_b -start 2012 10 -stop 2013 04 

4) Processing a single day

     python pysiral-l2proc.py -l2-settings cryosat2_awi_north -run-tag some_short_run -start 2016 03 01 -stop 2016 03 01

5) Do not write output data in a unique directory 

     python pysiral-l2proc.py -l2-settings cryosat2_awi_north -run-tag some_run -start 2015 10 -stop 2016 04 --no-overwrite-proctection

5) Do another run with same settings as in 5), but remove the previous data first (without the safety prompt)

     python pysiral-l2proc.py -l2-settings cryosat2_awi_north -run-tag some_run -start 2015 10 -stop 2016 04 --no-overwrite-proctection --remove-old --no-critical-prompt