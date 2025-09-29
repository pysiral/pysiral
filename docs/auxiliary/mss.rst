Mean Sea Surface (mss)
======================

Mean sea surface auxiliary data files are mostly global static files 
that contain the mean sea surface height relative to a certain ellipsoid. 

This information is required for deriving freeboard acurately in the
Level-2 processor.


DTU21 Mean Sea Surface - (DTU21)
--------------------------------

Configuration
^^^^^^^^^^^^^

The :term:`Auxiliary Dataset ID` is ``mss:dt21`` since pysiral version 0.12 and the configuration 
in the auxiliary data definition of the Level-2 processor definition file is:

.. code-block:: yaml

    - mss:
        name: dtu21
        options:
            latitude_range: [45.0, 90.0]

**Options**

- ``latitude_range``: The latitude range to be used for subsetting the global mean sea surface 
  data before mapping the content to the radar altimeter trajectory. 
  The latitude range is defined as `[<latitude_min>, <latitude_max>]` and the values in 
  decimal degrees need to be adapted to the target :term:`Hemisphere`.


Data Variables
^^^^^^^^^^^^^^

The following variables are returned from this dataset:

- ``mean_sea_surface``: Mean sea surface height in meters relative to the WGS84 ellipsoid


Data Source and Storage
^^^^^^^^^^^^^^^^^^^^^^^^

The DTU21 mean sea surface data is sourced from the DTU Space [DTU21-MSS]_, Technical University of Denmark. The required file is the global 1 minute grid in netCDF format and referenced to the TOPEX ellipsoid.

.. code-block::

    <source directory as defined in local_machine_def.yaml>
    └── DTU21MSS_1min.mss.nc


.. [DTU21-MSS] Andersen, Ole Baltazar (2022). DTU21 Mean Sea Surface. Technical University of Denmark. Dataset. https://doi.org/10.11583/DTU.19383221.v2
