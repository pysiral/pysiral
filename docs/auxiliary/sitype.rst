Sea Ice Type / Classification (sitype)
======================================

Several :term:`Sea Ice Type` products are available in pysiral
for ingestion in the :term:`Level-2 Processor`.


Global Sea Ice Type (Multi-Sensor) - (OSI-403)
------------------------------------------------

The OSI-403 sea ice type data product is an operational
product available with a :term:`timeliness` of only 5 hours 
and thus suitable for near real-time applications.

.. note:: 
    This dataset is not consistently generated for the full time coverage
    and several versions (`-a`, `-b`, `-d`, `-d`, ...) exist for different
    time periods. The use of this product should be reserved to 
    near-real-time applications and long-term pysiral generated 
    data records should be based on (interim) climate data records. 


Configuration Parameters
^^^^^^^^^^^^^^^^^^^^^^^^

The :term:`Auxiliary Dataset ID` is ``sitype:osi-403`` since pysiral version 0.12 and the configuration 
in the auxiliary data definition of the Level-2 processor definition file is:

.. code-block:: yaml

    - sitype:
        name: osi-403
        options:
            exception_on_error: False
            fill_valid_sic_gaps: True

**Options**

- ``fill_valid_sic_gaps``: Flag to determine whether gap-filling should be used. Setting the value to ``True`` ensures that there is valid sea ice type information is available everywhere there is sea ice indicated by the sea ice concentration information. This may be necessary in near coastal regions. The gap filling is based on nearest-neighbour interpolation within a certain range. If this range is exceeded, the `sea ice type` value will be set to `0.5` (ambiguous).


.. note:: 
    The dataset id was renamed in pysiral version 0.12 from ``sitype:osisaf-operational`` to ``sitype:osi-403``.


Data Variables
^^^^^^^^^^^^^^

The following variables are returned from this dataset:

- ``sea_ice_type``: Sea ice type as multi-year ice area fraction
- ``sea_ice_type_uncertainty``: Uncertainty estimate of sea ice type as multi-year ice area fraction



Data Source and Storage
^^^^^^^^^^^^^^^^^^^^^^^

The data access is specified on the `OSI-SAF website <https://osi-saf.eumetsat.int/products/osi-403-d>`_.
The expected grid and data format is the 10km polarstereographic grid 
in netCDF format. Data files from both hemispheres need to be organized in the following data structure: 

.. code-block::

    <source directory as defined in local_machine_def.yaml>
    └── <year>
        └── <month>
            ├── ice_type_nh_polstere-100_multi_<year><month><day>1200.nc
            ├── ice_type_sh_polstere-100_multi_<year><month><day>1200.nc
            ├── ...


C3S Sea Ice Type (Interim) Climate Data Record (C3S)
----------------------------------------------------

The C3S sea ice type data product [SITYPE-C3S-REF]_ is an (interim) climate data record that provides daily gridded data
on sea ice type from 1978 to present. It is based on satellite observations and is suitable for climate 
research and long-term monitoring of sea ice conditions. 

The climate data record (:term:`CDR`) and interim climate data record (:term:`iCDR`) are technically 
two separate datasets, but are treated as one in pysiral to allow seamless use of the `iCDR` after
the end of the `CDR` generation period. 

.. [SITYPE-C3S-REF] Copernicus Climate Change Service, Climate Data Store, (2020): Sea ice edge and type daily gridded data from 1978 to present derived from satellite observations. Copernicus Climate Change Service (C3S) Climate Data Store (CDS). DOI: 10.24381/cds.29c46d83 

.. note:: 

    The availability of this dataset is limited to the northern hemisphere until version 3.0. 
    Southern hemisphere data is added in v4.0. 


Configuration Parameters
^^^^^^^^^^^^^^^^^^^^^^^^

The :term:`Auxiliary Dataset ID` is ``sitype:c3s`` and the configuration 
in the auxiliary data definition of the Level-2 processor definition file is:

.. code-block:: yaml

    - sitype:
        name: c3s
        options:
            version: v3p0
            fill_valid_sic_gaps: True
            exception_on_error: False

**Options**

- ``version``: The version of the dataset to be used. This should be set to the desired version string, e.g. `v3p0`.
- ``fill_valid_sic_gaps``: Flag to determine whether gap-filling should be used. Setting the value to `True` ensures that there is valid sea ice type information is available everywhere there is sea ice indicated by the sea ice concentration information. This may be necessary in near coastal regions. The gap filling is based on nearest-neighbour interpolation within a certain range. If this range is exceeded, the `sea ice type` value will be set to `0.5` (ambiguous).
- ``exception_on_error``: Flag, when set to `True`, will raise an exception if the auxiliary dataset is not available
  for the requested time period. If set to `False`, the processor will continue without the auxiliary dataset and
  all output is set to `NaN`. An error will be logged in this case.


Data Variables
^^^^^^^^^^^^^^

The following variables are returned from this dataset:

- ``sea_ice_type``: Sea ice type as multi-year ice area fraction
- ``sea_ice_type_uncertainty``: Uncertainty estimate of sea ice type as multi-year ice area fraction


Data Source and Storage
^^^^^^^^^^^^^^^^^^^^^^^

Both `CDR` and `iCDR` data files need to be organized in a common directory structure to allow 
pysiral automatic file lookup depending on the target date. 

.. code-block::

    <source directory as defined in local_machine_def.yaml>
    └── <data record type: cdr|icdr>
        └── <version: v{major}p{minor}>      # e.g v3p0
            └── <year>
                └── <month>
                    ├── ice_type_<hemisphere>_ease2-250_<cdr|icdr>-<version>_<year><month><day>1200.nc
                    ├── ...

There are two options of downloading the datatset: 

1. Climate Data Store: https://cds.climate.copernicus.eu/datasets/satellite-sea-ice-edge-type
2. Thredds Server at MET Norway: https://thredds.met.no/thredds/c3s/c3s.html (Select 
   `Sea Ice Type Climate Data Record (CDR)` and `Sea Ice Type Interim Climate Data Record (iCDR)` 
   of the target version).
