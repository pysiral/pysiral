Sea Ice Concentration (sic)
===========================

Several :term:`Sea Ice Concentration` products are available in pysiral
for ingestion in the :term:`Level-2 Processor`.


.. _SIC_OSI401:

Global Sea Ice Concentration (SSMIS) - (OSI-401)
------------------------------------------------

The OSI-401 sea ice concentration data product is an operational
product available with a :term:`timeliness` of only 5 hours and thus
suitable for near real-time applications.

.. note:: 
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

.. note:: 
    The dataset id was renamed in pysiral version 0.12 from ``sic:osisaf-operational`` to ``sic:osi-401``.


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

.. note:: 
    This dataset is not consistently generated for the full time coverage
    and several versions (`-a`, `-b`, `-d`, `-d`, ...) exist for different
    time periods. The use of this product should be reserved to 
    near-real-time applications and long-term pysiral generated 
    data records should be based on (interim) climate data records
    sea ice concentration. 


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


Global Sea Ice Concentration (i)CDR (AMSR2), (OSI-438/458)
----------------------------------------------------------
