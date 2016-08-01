# PysiralVisualization

## Visualization of Level-1 pre-processed orbit segments

The content of the pre-processed l1bdata netCDF files can be summarized with `tools\l1bdata-report.py`. This tools creates a pdf report including a table listing of the metadata and graphical visualization of the time-orbit, waveform, surface-type and classifier data groups. 

### Usage

    python l1bdata_report.py [options] l1bdata_netcdf_file(s)

### l1bdata input files

The input files argument can either be a link to one speficic file or a simple search pattern, e.g.:

    python l1bdata_report.py E:\altim\data\altimetry\cryosat2\baseline-c\l1bdata\north\2015\04\l1bdata_v00_north_cryosat2_026390_20150401T003811_20150401T003915.nc

or

    python l1bdata_report.py E:\altim\data\altimetry\cryosat2\baseline-c\l1bdata\north\2015\04\*.nc

### Options

The options can be used to modify the plot style, add labels and also to create difference plots. 

Argument | Required | Parameter | Description
--- | --- | --- | --- | ---
**-h** | optional | *None* | Display manual
**-o** *output_directory* | optional | str | valid target directory for pdf reports. If omitted, a `report` subdirectory of the input directory will be created 




*****




## Visualization of Level-3 gridded data

Visualization of L3 data products (and differences) is handled by `tools\create_l3s_map.py`. It is based on the pysiral.visualization module. The definition of the colorbar (colormap, limits) is taken from the parameter configuration files in `config\parameter_def.yaml`. 

### Usage

    python create_l3s_map.py [options] source_file(s)

### Source file

The source files is either a full path to L3s netCDF file or a search pattern recognized by the python glob module, e.g.:
   
    python create_l3s_map.py [options] E:\altim\products\altimetry\pysiral\some_result\l3s\*.nc

If a search pattern is used, a map will be produced for all netCDF files in this folder.

### Options

The options can be used to modify the plot style, add labels and also to create difference plots. 

Argument | Required | Parameter | Description
--- | --- | --- | --- | ---
**-h** | optional | *None* | Display manual
**-parameter** *parameter list* | yes |string list (e.g. '-parameter freeboard sea_ice_thickness') |  list of pysiral recognized l3 parameters separated by blanks 
**-o -output** *folder* | optional | valid folder | destination folder for plots (default: `l3_grid_folder\visuals`)
**-type** *presentation type* | optional | string | plot style, either `presentation` (default) or `paper`
**-diff** *filename*  | optional | valid filename (abolute path) | will produce a difference plot with *filename* as reference grid
**-annotation** *label string in * | optional | quoted string (e.g. `-annotation "DTU15 - DTU13"`)  | additional label that will be included in output filename and the map 
**--cs2awi** | optional | *None* | source (and diff) file originates from cs2awi IDL processor
**--cb** | optional (default) | *None* |  draw the colorbar 
**--no-cb** | optional | *None* | do not draw the colorbar

### Cookbook

1) Create a map of freeboard, thickness and sea surface height anomaly of a given l3s grid

     python create_l3s_map.py -parameter freeboard sea_ice_thickness sea_surface_anomaly E:\some_product\l3s\2016\l3s_v00_cryosat2_EASE2_25km_20160401_20160430.nc

Here the output will be written to `E:\altim\products\altimetry\pysiral\some_result\l3s\2016\visuals\presentation`

2) Batch script to produce a series of freeboard difference plots in paper style with custom output folder and annotations:

     @ECHO OFF

     SET pyargs=-W ignore
     SET bin=D:\workspace\pysiral\tools\create_l3s_map.py
     SET gridDir=E:\altim\product\altimetry\cs2awi\
     SET refGrid=E:\altim\product\altimetry\cs2awi\arctic-baselineC-ucl13\grid\cs2awi_nh_201504.nc
     SET output=X:\paper\2016-MSS-Intercomparison\figures\mss-freeboards-diff
     SET args=-parameter sea_ice_freeboard -type paper -output %output% --cs2awi

     python %pyargs% %bin% %args% -annotation "UCL13 - No-MSS" -diff %refGrid% %gridDir%arctic-baselineC-nomss\grid\cs2awi_nh_201504.nc
     python %pyargs% %bin% %args% -annotation "UCL13 - EGM2008" -diff %refGrid% %gridDir%arctic-baselineC-egm2008\grid\cs2awi_nh_201504.nc
     python %pyargs% %bin% %args% -annotation "UCL13 - DTU10" -diff %refGrid% %gridDir%arctic-baselineC-dtu10\grid\cs2awi_nh_201504.nc
     python %pyargs% %bin% %args% -annotation "UCL13 - DTU13" -diff %refGrid% %gridDir%arctic-baselineC-dtu13\grid\cs2awi_nh_201504.nc
     python %pyargs% %bin% %args% -annotation "UCL13 - DTU15" -diff %refGrid% %gridDir%arctic-baselineC-dtu15\grid\cs2awi_nh_201504.nc
