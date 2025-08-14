Sea Ice Concentration
=====================

Several :term:`Sea Ice Concentration` products are available in pysiral
for ingestion in the :term:`Level-2 Processor`.


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

The auxiliary dataset id is ``sic:osi-401`` (for pysiral versions >= 0.12) and the configuration 
in the auxiliary data definition of the Level-2 processor definition file is:

.. code-block:: yaml

    - sic:
        name: osi-401
        options:
            exception_on_error: False

.. note:: 
    The dataset id was renamed in pysiral version 0.12 from ``sic:osisaf-operational`` to ``sic:osi-401``.


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



Global Sea Ice Concentration (i)CDR, (OSI-430/450)
--------------------------------------------------


Global Sea Ice Concentration (i)CDR (AMSR2), (OSI-438/458)
----------------------------------------------------------
