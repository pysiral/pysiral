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
**--remove-old** | optional | *None* | erase all files in output directory before pre-processing **(not yet implemented)**

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