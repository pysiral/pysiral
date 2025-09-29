Sea Ice Concentration (sic)
===========================

Several :term:`Sea Ice Concentration` products are available in pysiral
for ingestion in the :term:`Level-2 Processor`. The code is implemented in :py:mod:`pysiral.auxdata.sic`.

.. _SIC_OSI401:

Global Sea Ice Concentration (SSMIS) - (OSI-401)
------------------------------------------------

The OSI-401 sea ice concentration data product is an operational
product available with a :term:`timeliness` of only 5 hours and thus
suitable for near real-time applications.

.. versionchanged:: 0.11

    The dataset id was renamed from ``sic:osisaf-operational`` to ``sic:osi-401``.

.. tip:: 
    This dataset is not consistently generated for the full time coverage
    and several versions (`-a`, `-b`, `-d`, `-d`, ...) exist for different
    time periods. The use of this product should be reserved to 
    near-real-time applications and long-term pysiral generated 
    data records should be based on (interim) climate data records
    sea ice concentration. 

.. warning:: 
    Due to the end of data distribution of SSMIS data in September 2026, 
    it is recommended to use the AMSR2-based sea ice concentration products, 
    such as ``sic:osi-408`` or ``sic:osi-cdr-amsr2``.


Configuration Parameters
^^^^^^^^^^^^^^^^^^^^^^^^

The :term:`Auxiliary Dataset ID` is ``sic:osi-401`` since pysiral version 0.12 and the configuration 
in the auxiliary data definition of the Level-2 processor definition file is:

.. code-block:: yaml

    - sic:
        name: osi-401
        options:
            exception_on_error: False

**Options**

To be passed to :py:class:`pysiral.auxdata.sic.OsiSafSIC`.

- ``exception_on_error``: Flag, when set to `True`, will raise an exception if the auxiliary dataset is not available
  for the requested time period. If set to `False`, the processor will continue without the auxiliary dataset and
  all output is set to `NaN`. An error will be logged in this case.



Data Variables
^^^^^^^^^^^^^^

The following variables are returned from this dataset:

- ``sea_ice_concentration``: Sea ice concentration in percent
- ``distance_to_ocean``: Distance to the nearest ocean grid cell in meters. This variable is defined as the distance to the nearest grid cell with `sea_ice_concentration < 15%`
- ``distance_to_low_ice_concentration``: Similar to `distance_to_ocean``, but distance to the nearest grid cell with `sea_ice_concentration < 70%`.



Data Source and Storage
^^^^^^^^^^^^^^^^^^^^^^^

The data access is specified on the `OSI-SAF website <https://osi-saf.eumetsat.int/products/osi-401-d>`_.
The expected grid and data format is the 10km polarstereographic grid 
in netCDF format. Data files from both hemispheres need to be organized in the following data structure: 

.. code-block::

    <source directory as defined in local_machine_def.yaml>
    └── <year>
        └── <month>
            ├── ice_conc_nh_polstere-100_multi_<year><month><day>1200.nc
            ├── ice_conc_sh_polstere-100_multi_<year><month><day>1200.nc
            ├── ...


Global Sea Ice Concentration (AMRS2), (OSI-408)
-----------------------------------------------


The OSI-408 sea ice concentration data product is an operational
product similar to :ref:`SIC_OSI401`, but based on AMSR2 input data. 
Datafiles are read with the same python class and handling and 
extracted variables are similar to the OSI-401 product.

.. tip:: 
    This dataset is not consistently generated for the full time coverage
    and several versions (`-a`, `-b`, `-d`, `-d`, ...) exist for different
    time periods. The use of this product should be reserved to 
    near-real-time applications and long-term pysiral generated 
    data records should be based on (interim) climate data records
    sea ice concentration. 

.. versionadded:: 0.11

    This auxiliary dataset was added in pysiral version 0.11 


Configuration Parameters
^^^^^^^^^^^^^^^^^^^^^^^^

The :term:`Auxiliary Dataset ID` is ``sic:osi-408`` since pysiral version 0.12 and the configuration 
in the auxiliary data definition of the Level-2 processor definition file is:

.. code-block:: yaml

    - sic:
        name: osi-408
        options:
            exception_on_error: False
            fill_pole_hole:
                pole_hole_lat_threshold: 87.0
                pole_hole_fill_value: 100.

The main difference to `sic:osi-401` is that the AMSR2 data contains a small pole hole that can be filled
with a static value. 

**Options**

- ``exception_on_error``: Flag, when set to `True`, will raise an exception if the auxiliary dataset is not available
  for the requested time period. If set to `False`, the processor will continue without the auxiliary dataset and
  all output is set to `NaN`. An error will be logged in this case.
- ``fill_pole_hole``: Options to fill the pole hole in the AMSR2 data.
  
    - ``pole_hole_lat_threshold``: Latitude threshold in decimal degrees for the pole hole filling (default: `87.0`).
    - ``pole_hole_fill_value``: Sea ice concentration fill value for the pole hole (default: `100.`).


Data Variables
^^^^^^^^^^^^^^

The following variables are returned from this dataset:

- ``sea_ice_concentration``: Sea ice concentration in percent
- ``distance_to_ocean``: Distance to the nearest ocean grid cell in meters. This variable is defined as the distance to the nearest grid cell with `sea_ice_concentration < 15%`
- ``distance_to_low_ice_concentration``: Similar to `distance_to_ocean``, but distance to the nearest grid cell with `sea_ice_concentration < 70%`.



Data Source and Storage
^^^^^^^^^^^^^^^^^^^^^^^

The data access is specified on the `OSI-SAF website <https://osi-saf.eumetsat.int/products/osi-401-d>`_.
The expected grid and data format is the 10km polarstereographic grid 
in netCDF format. Data files from both hemispheres need to be organized in the following data structure: 

.. code-block::

    <source directory as defined in local_machine_def.yaml>
    └── <year>
        └── <month>
            ├── ice_conc_nh_polstere-100_multi_<year><month><day>1200.nc
            ├── ice_conc_sh_polstere-100_multi_<year><month><day>1200.nc
            ├── ...


Global Sea Ice Concentration (i)CDR, (OSI-430/450)
--------------------------------------------------

.. admonition:: ToDo
  
    To be documented

Global Sea Ice Concentration (i)CDR (AMSR2), (OSI-438/458)
----------------------------------------------------------

.. admonition:: ToDo
  
    To be documented


C3S Sea Ice Concentration (Interim) Climate Data Record (C3S)
-------------------------------------------------------------

The C3S sea ice concentration data product [SIC-C3S-REF]_ is an (interim) climate data record that provides daily gridded data
on sea ice concentration from 1978 to present. It is based on satellite observations and is suitable for climate
research and long-term monitoring of sea ice conditions.

The climate data record (:term:`CDR`) and interim climate data record (:term:`iCDR`) are technically 
two separate datasets, but are treated as one in pysiral to allow seamless use of the `iCDR` after
the end of the `CDR` generation period. 

.. [SIC-C3S-REF] Copernicus Climate Change Service (C3S) (2020): Sea ice concentration daily gridded data from 1978 to present derived from satellite observations. Copernicus Climate Change Service (C3S) Climate Data Store (CDS). DOI: 10.24381/cds.3cd8b812 



Configuration Parameters
^^^^^^^^^^^^^^^^^^^^^^^^

The :term:`Auxiliary Dataset ID` is ``sic:c3s`` and the configuration 
in the auxiliary data definition of the Level-2 processor definition file is:

.. code-block:: yaml

    - sic:
        name: c3s
        options:
            version: v3p0
            exception_on_error: False

**Options**

- ``version``: The version of the dataset to be used. This should be set to the desired version string, e.g. `v3p0`.
- ``exception_on_error``: Flag, when set to `True`, will raise an exception if the auxiliary dataset is not available
  for the requested time period. If set to `False`, the processor will continue without the auxiliary dataset and
  all output is set to `NaN`. An error will be logged in this case.


Data Variables
^^^^^^^^^^^^^^

The following variables are returned from this dataset:

- ``sea_ice_concentration``: Sea ice concentration in percent
- ``distance_to_ocean``: Distance to the nearest ocean grid cell in meters. This variable is defined as the distance to the nearest grid cell with `sea_ice_concentration < 15%`
- ``distance_to_low_ice_concentration``: Similar to `distance_to_ocean``, but distance to the nearest grid cell with `sea_ice_concentration < 70%`.



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
                    ├── ice_conc_<hemisphere>_ease2-250_<cdr|icdr>-<version>_<year><month><day>1200.nc
                    ├── ...

There are two options of downloading the datatset: 

1. Climate Data Store: https://cds.climate.copernicus.eu/datasets/sea-ice-concentration
2. Thredds Server at MET Norway: https://thredds.met.no/thredds/c3s/c3s.html (Select 
   `Sea Ice Concentration Climate Data Record (CDR)` and `Sea Ice Concentration Interim Climate Data Record (iCDR)` 
   of the target version).
